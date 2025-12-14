# Updated PPI analysis app with clustering coefficient, KEGG validation,
# cancer-specific PPI export, improved visualizations, and validation sheet.
# Place this file alongside your original resources and run as before.

import os
import gzip
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import webbrowser
import json
from collections import Counter
import threading
from functools import wraps
import requests
import time

try:
    import community.community_louvain as community_louvain
except Exception:
    community_louvain = None

try:
    from Bio import Entrez
except Exception:
    Entrez = None


def read_maybe_gzip(path, sep=None, **kwargs):
    open_fn = gzip.open if str(path).endswith('.gz') else open

    if sep is None:
        sep = r'\s+'

    try:
        with open_fn(path, 'rt', encoding='utf-8') as fh:
            df = pd.read_csv(fh, sep=sep, **kwargs)
        return df
    except Exception as e:
        raise e


def strip_taxon_prefix(s):
    """Strip numeric taxon prefix like '9606.' if present."""
    if isinstance(s, str) and '.' in s:
        return s.split('.', 1)[1]
    return s


def load_string_links(links_path, score_threshold=0, sample_fraction=None):
    """Load STRING links file and coerce into standard columns."""
    path = links_path
    open_fn = gzip.open if str(path).endswith('.gz') else open
    with open_fn(path, 'rt', encoding='utf-8') as fh:
        sample = fh.read(8192)
    delim = None
    if '\t' in sample:
        delim = '\t'
    elif ',' in sample:
        delim = ','
    else:
        delim = r'\s+'
    try:
        if delim == r'\s+':
            df = pd.read_csv(links_path, sep=delim, engine='python', compression='infer', header=0)
        else:
            df = pd.read_csv(links_path, sep=delim, compression='infer', header=0)
    except Exception:
        df = pd.read_csv(links_path, compression='infer', engine='python', header=0)

    # Normalize columns
    if not any(k in df.columns for k in ['protein1', 'protein2', 'protein1', 'protein2']):
        if df.shape[1] >= 3:
            df = df.iloc[:, 0:3]
            df.columns = ['protein1', 'protein2', 'combined_score']

    score_col = None
    for n in ['combined_score', 'score']:
        if n in df.columns:
            score_col = n
            break
    if score_col is None:
        df['combined_score'] = 0
        score_col = 'combined_score'

    links = df[['protein1', 'protein2', score_col]].copy()
    links.columns = ['protein1', 'protein2', 'combined_score']

    # Strip taxon prefix if present (e.g., "9606.ENSP0000...")
    links['p1_short'] = links['protein1'].astype(str).map(strip_taxon_prefix)
    links['p2_short'] = links['protein2'].astype(str).map(strip_taxon_prefix)
    # Filter by score threshold if provided (STRING score is typically 0-1000)
    if score_threshold and 'combined_score' in links.columns:
        links = links[links['combined_score'].astype(float) >= float(score_threshold)]
    # Optionally sample
    if sample_fraction and 0 < sample_fraction < 1.0:
        links = links.sample(frac=sample_fraction, random_state=42)

    return links


def load_string_info(path):
    df = read_maybe_gzip(path, sep='\t', header=0, dtype=str)
    cols = [c.strip() for c in df.columns]
    df.columns = cols
    id_col = None
    gene_col = None
    candidates_id = ['protein_external_id', 'protein_id', 'protein', 'string_protein_id', 'string_id']
    candidates_gene = ['preferred_name', 'gene', 'gene_name', 'preferred_name']
    for c in candidates_id:
        if c in df.columns:
            id_col = c
            break
    for c in candidates_gene:
        if c in df.columns:
            gene_col = c
            break
    if id_col is None:
        id_col = df.columns[0]
    if gene_col is None:
        gene_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
    df = df[[id_col, gene_col]].drop_duplicates()
    df.columns = ['string_id', 'gene']
    df['short_id'] = df['string_id'].map(strip_taxon_prefix)
    df['gene_upper'] = df['gene'].str.upper()
    return df


def load_oncokb_tsv(path):
    df = read_maybe_gzip(path, sep='\t', header=0, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    if 'HUGO_SYMBOL' in df.columns:
        col = 'HUGO_SYMBOL'
    elif 'Gene' in df.columns:
        col = 'Gene'
    else:
        col = df.columns[0]
    genes = df[col].dropna().astype(str).str.strip().str.upper().unique()
    return set(genes), col

def load_cosmic_tsv(path):
    df = read_maybe_gzip(path, sep='\t', header=0, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    for c in ['Gene Symbol', 'GENE_SYMBOL', 'Gene']:
        if c in df.columns:
            return set(df[c].dropna().str.upper().unique())
    return set()

def load_tcga_gene_list(path):
    df = read_maybe_gzip(path, sep='\t', header=0, dtype=str)
    return set(df.iloc[:, 0].dropna().str.upper().unique())

def load_biogrid(path):
    """Load BioGRID tab2 if provided (optional)."""
    df = read_maybe_gzip(path, sep='\t', header=0, dtype=str)
    possible_pairs = [('Official Symbol Interactor A', 'Official Symbol Interactor B'),
                      ('Official Symbol A', 'Official Symbol B'),
                      ('Interactor A', 'Interactor B'),
                      ('Gene A', 'Gene B')]
    a_col = b_col = None
    for a, b in possible_pairs:
        if a in df.columns and b in df.columns:
            a_col, b_col = a, b
            break
    if a_col is None:
        a_col, b_col = df.columns[0], df.columns[1]
    return df[[a_col, b_col]].dropna().rename(columns={a_col: 'A', b_col: 'B'})

def safe_save_df(df, path):
    try:
        df.to_csv(path, index=False)
        return True
    except Exception as e:
        print("Failed to save:", e)
        return False

# Threading helper
def run_in_thread(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        self = args[0] if args else None

        def target():
            try:
                if hasattr(self, 'set_busy'):
                    self.set_busy(True)
                fn(*args, **kwargs)
            finally:
                if hasattr(self, 'set_busy'):
                    self.set_busy(False)

        t = threading.Thread(target=target, daemon=True)
        t.start()
        return t

    return wrapper

class PPINetworkAnalysis:
    def __init__(self):
        self.string_links = None
        self.string_info = None
        self.oncokb = set()
        self.oncokb_column = None
        self.biogrid = None
        self.G = None
        self.G_sub = None
        self.metrics_df = None
        self.community = None
        self.params = {
            'score_threshold': 700,
            'include_biogrid': False,
            'use_weighted_degree': True,
            'louvain_resolution': 1.0,
            'sample_fraction': 0.01
        }
        self.cosmic = set()
        self.tcga = set()

    def load_cosmic(self, path):
        self.cosmic = load_cosmic_tsv(path)

    def load_tcga(self, path):
        self.tcga = load_tcga_gene_list(path)

    def load_string(self, links_path, info_path=None):
        self.string_links = load_string_links(links_path, self.params.get('score_threshold', 0),
                                             self.params.get('sample_fraction'))
        if info_path is not None:
            print("Loading protein info mapping...")
            self.string_info = load_string_info(info_path)
        else:
            self.string_info = None

    def load_oncokb(self, oncokb_path):
        s, col = load_oncokb_tsv(oncokb_path)
        self.oncokb = set([g.upper() for g in s])
        self.oncokb_column = col

    def detect_label_propagation(self, graph=None):
        if graph is None:
            graph = self.G_sub if self.G_sub is not None else self.G

        communities = nx.algorithms.community.label_propagation_communities(graph)
        partition = {}
        for i, comm in enumerate(communities):
            for node in comm:
                partition[node] = i

        nx.set_node_attributes(graph, partition, 'label_propagation')
        return partition

    def detect_girvan_newman(self, graph=None, level=1):
        if graph is None:
            graph = self.G_sub if self.G_sub is not None else self.G

        gen = nx.algorithms.community.girvan_newman(graph)
        communities = next(gen)
        for _ in range(level - 1):
            communities = next(gen)

        partition = {}
        for i, comm in enumerate(communities):
            for node in comm:
                partition[node] = i

        nx.set_node_attributes(graph, partition, 'girvan_newman')
        return partition
    
    def detect_cosmic_seeded_modules(self, graph=None, hop=1):
        if graph is None:
            graph = self.G_sub if self.G_sub is not None else self.G

        modules = {}
        cid = 0

        for seed in self.cosmic:
            if seed not in graph:
                continue

            nodes = {seed}
            frontier = {seed}
            for _ in range(hop):
                next_frontier = set()
                for n in frontier:
                    next_frontier.update(graph.neighbors(n))
                nodes.update(next_frontier)
                frontier = next_frontier

            for n in nodes:
                if n not in modules:
                    modules[n] = cid
            cid += 1

        nx.set_node_attributes(graph, modules, 'cosmic_module')
        return modules
    
    def annotate_tcga_nodes(self, graph=None):
        if graph is None:
            graph = self.G_sub if self.G_sub is not None else self.G

        attrs = {n: {'tcga_gene': (n in self.tcga)} for n in graph.nodes()}
        nx.set_node_attributes(graph, attrs)

    def construct_graph(self):
        if self.string_links is None:
            raise ValueError("STRING links not loaded.")

        print("Building network graph")
        df = self.string_links.copy()

        # Map to gene names if info available
        if self.string_info is not None:
            gene_map = self.string_info[['short_id', 'gene_upper']].drop_duplicates()
            df = df.merge(gene_map, left_on='p1_short', right_on='short_id', how='left')
            df = df.merge(gene_map, left_on='p2_short', right_on='short_id', how='left', suffixes=('_1', '_2'))
            df['g1'] = df['gene_upper_1'].fillna(df['p1_short']).str.upper()
            df['g2'] = df['gene_upper_2'].fillna(df['p2_short']).str.upper()
        else:
            df['g1'] = df['p1_short'].str.upper()
            df['g2'] = df['p2_short'].str.upper()

        df = df[df['g1'] != df['g2']]

        df['weight'] = pd.to_numeric(df['combined_score'], errors='coerce').fillna(0.0) / 1000.0

        print(f"Creating graph with {len(df)} edges...")
        G = nx.from_pandas_edgelist(df, source='g1', target='g2', edge_attr='weight')
        if self.params.get('include_biogrid') and self.biogrid is not None:
            for _, row in self.biogrid.iterrows():
                a = str(row['A']).upper()
                b = str(row['B']).upper()
                if a == b:
                    continue
                if not G.has_edge(a, b):
                    G.add_edge(a, b, weight=0.001)
        self.G = G
        return G

    @run_in_thread
    def compute_metrics(self):
        if self.G is None:
            raise ValueError("Graph not constructed yet.")
        G = self.G
        deg = dict(G.degree())
        wdeg = dict(G.degree(weight='weight'))
        try:
            pr = nx.pagerank(G, weight='weight')
        except Exception:
            pr = {n: 0.0 for n in G.nodes()}
        df = pd.DataFrame({
            'node': list(G.nodes()),
            'degree': [deg[n] for n in G.nodes()],
            'wdegree': [wdeg[n] for n in G.nodes()],
            'pagerank': [pr[n] for n in G.nodes()]
        })
        df['is_oncokb'] = df['node'].str.upper().isin(self.oncokb)
        n_nodes = G.number_of_nodes()
        # Compute betweenness & closeness only if not too large
        if n_nodes <= 2000:
            print("Computing betweenness & closeness (may take a minute)...")
            bet = nx.betweenness_centrality(G, weight='weight', normalized=True)
            clo = nx.closeness_centrality(G)
            df['betweenness'] = df['node'].map(bet)
            df['closeness'] = df['node'].map(clo)
        else:
            print(f"Graph has {n_nodes} nodes — skipping betweenness/closeness (use smaller subgraph).")
            df['betweenness'] = 0.0
            df['closeness'] = 0.0

        # Clustering coefficient (added per requirement)
        try:
            cc = nx.clustering(G)
            df['clustering'] = df['node'].map(cc)
        except Exception:
            df['clustering'] = 0.0

        for col in ['degree', 'wdegree', 'pagerank', 'betweenness', 'closeness', 'clustering']:
            if col in df.columns:
                df[f'{col}_rank'] = df[col].rank(method='min', ascending=False)

        self.metrics_df = df.sort_values('degree', ascending=False).reset_index(drop=True)

        if self.metrics_df is not None:
            self.metrics_df['is_cosmic'] = self.metrics_df['node'].isin(self.cosmic)
            self.metrics_df['is_tcga'] = self.metrics_df['node'].isin(self.tcga)
            
        return self.metrics_df

    def extract_cancer_subgraph(self, hop=1, min_degree=0):
        """
        Build a cancer-centered subgraph:
        - nodes: OncoKB-mapped proteins present in G, plus their neighbors up to 'hop'
        - prunes nodes below min_degree (in the subgraph)
        """
        if self.G is None:
            raise ValueError("Graph not constructed yet.")
        G = self.G
        onc_nodes = [n for n in G.nodes() if n.upper() in self.oncokb]

        print(f"Found {len(onc_nodes)} cancer genes in the network")
        if len(onc_nodes) == 0:
            print("WARNING: No cancer genes found in network. This might happen with very small samples.")
            print("Try increasing sample_fraction or score_threshold.")

        nodes = set(onc_nodes)
        if hop >= 1:
            for n in list(nodes):
                nodes.update(G.neighbors(n))
        if hop > 1:
            frontier = set(nodes)
            for _ in range(hop - 1):
                new_frontier = set()
                for n in frontier:
                    new_frontier.update(G.neighbors(n))
                nodes.update(new_frontier)
                frontier = new_frontier
        sub = G.subgraph(nodes).copy()

        if min_degree > 0:
            remove = [n for n, d in dict(sub.degree()).items() if d < min_degree and n.upper() not in self.oncokb]
            sub.remove_nodes_from(remove)

        cancer_in_sub = sum(1 for n in sub.nodes() if n.upper() in self.oncokb)
        print(f"Cancer subgraph: {sub.number_of_nodes()} total nodes ({cancer_in_sub} cancer genes)")

        self.G_sub = sub
        return sub

    #@run_in_thread
    def detect_communities(self, graph=None, resolution=1.0):
        if graph is None:
            graph = self.G_sub if self.G_sub is not None else self.G
        if community_louvain is None:
            raise RuntimeError("python-louvain (community) package is not installed.")
        partition = community_louvain.best_partition(graph, weight='weight', resolution=resolution)
        nx.set_node_attributes(graph, partition, 'louvain_community')
        self.community = partition
        return partition

    def export_gexf(self, graph=None, path='ppi_export.gexf'):
        if graph is None:
            graph = self.G_sub if self.G_sub is not None else self.G
        nx.write_gexf(graph, path)
        return path

    def generate_pyvis(self, graph=None, path='ppi_interactive.html', show_in_browser=False, highlight_oncokb=True,
                       max_nodes=500):
        if graph is None:
            graph = self.G_sub if self.G_sub is not None else self.G

        n_nodes = graph.number_of_nodes()
        if n_nodes > max_nodes:
            print(f"Graph has {n_nodes} nodes. Reducing to top {max_nodes} by degree for visualization...")
            degrees = dict(graph.degree())
            top_nodes = sorted(degrees.keys(), key=lambda x: degrees[x], reverse=True)[:max_nodes]
            graph = graph.subgraph(top_nodes).copy()
            print(f"Visualizing {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")

        net = Network(height='900px', width='100%', notebook=False, bgcolor='#222222', font_color='white')
        net.barnes_hut()

        for n in graph.nodes():
            degree = graph.degree(n)
            title = f"{n}<br>degree={degree}"
            is_onco = (n.upper() in self.oncokb) if highlight_oncokb else False
            color = '#ff4444' if is_onco else '#4488ff'
            size = 10 + min(degree * 2, 40)
            net.add_node(n, label=(n if is_onco else ''), title=title, color=color, size=size)

        for u, v, data in graph.edges(data=True):
            w = data.get('weight', 1.0)
            net.add_edge(u, v, value=float(w))

        net.save_graph(path)
        if show_in_browser:
            webbrowser.open(path)
        return path

    def top_hubs(self, topn=50, by='degree'):
        if self.metrics_df is None:
            raise ValueError("Metrics not computed yet.")
        if by not in self.metrics_df.columns:
            by = 'degree'
        return self.metrics_df.sort_values(by, ascending=False).head(topn)

    def query_pubmed(self, gene_symbol, email, retmax=5, sleep=0.4):
        """Query PubMed for a specific gene symbol. Returns list of (pmid, title, url).
           Requires Biopython Entrez and internet access. User must provide email.
        """
        if Entrez is None:
            raise RuntimeError("Biopython is not installed; PubMed queries unavailable.")
        Entrez.email = email
        term = f"{gene_symbol}[Title/Abstract] AND cancer[Title/Abstract]"
        handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
        record = Entrez.read(handle)
        ids = record.get('IdList', [])
        results = []
        if not ids:
            return results
        # small sleep to be polite to NCBI
        time.sleep(sleep)
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        recs = Entrez.read(handle)
        for article in recs.get('PubmedArticle', []):
            try:
                pmid = article['MedlineCitation']['PMID']
                art = article['MedlineCitation']['Article']
                title = art.get('ArticleTitle', '')
                url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                results.append((str(pmid), title, url))
            except Exception:
                continue
        return results

    def kegg_gene_pathways(self, gene_symbol):
        """
        Query KEGG REST to find pathways for a human gene symbol.
        This function uses KEGG REST service and will fail gracefully if no internet.
        Returns a list of KEGG pathway entries (strings) or empty list.
        """
        try:
            # Step 1: find human gene entry in KEGG by name
            # The KEGG 'find' endpoint can search gene names. We'll search "hsa:GENE" as well as plain GENE.
            # Try searching "hsa:GENE" mapping; often KEGG uses ids like 'hsa:7157' for TP53.
            # Try to find 'hsa:GENE' directly:
            resp = requests.get(f"http://rest.kegg.jp/find/genes/{gene_symbol}", timeout=10)
            if resp.status_code != 200:
                return []
            lines = resp.text.strip().splitlines()
            # Filter lines for human (hsa:)
            hsa_lines = [ln for ln in lines if '\thsa:' in ln or ln.startswith('hsa:') or 'Homo sapiens' in ln]
            entry_ids = []
            for ln in hsa_lines:
                parts = ln.split('\t')
                if parts:
                    # first token can be "hsa:####" or "gene_id\thsa:####"
                    token = parts[0].strip()
                    if token.startswith('hsa:'):
                        entry_ids.append(token)
                    else:
                        # try to find hsa: in rest
                        for p in parts:
                            if 'hsa:' in p:
                                entry_ids.append(p.strip().split()[0])
                                break
            # fallback: parse first line tokens for hsa:
            if not entry_ids and lines:
                for ln in lines:
                    if 'hsa:' in ln:
                        # find token with hsa:
                        for t in ln.split():
                            if t.startswith('hsa:'):
                                entry_ids.append(t)
                                break
            pathways = set()
            for eid in entry_ids:
                # Link to pathways
                r2 = requests.get(f"http://rest.kegg.jp/link/pathway/{eid}", timeout=10)
                if r2.status_code != 200:
                    continue
                for l in r2.text.strip().splitlines():
                    parts = l.split('\t')
                    if len(parts) == 2:
                        # second part is e.g., path:hsa05200
                        pathways.add(parts[1].strip())
            # Convert pathway ids to human-readable names
            pathway_names = []
            for pid in pathways:
                r3 = requests.get(f"http://rest.kegg.jp/get/{pid}", timeout=10)
                if r3.status_code != 200:
                    continue
                # parse NAME field
                for ln in r3.text.splitlines():
                    if ln.startswith("NAME"):
                        pname = ln.replace("NAME", "").strip()
                        pathway_names.append(f"{pid} -- {pname}")
                        break
            return pathway_names
        except Exception:
            return []

    def save_cancer_specific_ppi(self, outpath):
        """Export a CSV of edges (and simple metadata) that include OncoKB genes."""
        if self.string_links is None:
            raise ValueError("STRING links not loaded.")
        df = self.string_links.copy()
        # Map to gene symbols if possible
        if self.string_info is not None:
            gene_map = self.string_info[['short_id', 'gene_upper']].drop_duplicates()
            df = df.merge(gene_map, left_on='p1_short', right_on='short_id', how='left')
            df = df.merge(gene_map, left_on='p2_short', right_on='short_id', how='left', suffixes=('_1', '_2'))
            df['g1'] = df['gene_upper_1'].fillna(df['p1_short']).str.upper()
            df['g2'] = df['gene_upper_2'].fillna(df['p2_short']).str.upper()
        else:
            df['g1'] = df['p1_short'].str.upper()
            df['g2'] = df['p2_short'].str.upper()

        # keep edges where at least one end is OncoKB
        onc_set = set([g.upper() for g in self.oncokb])
        mask = df['g1'].str.upper().isin(onc_set) | df['g2'].str.upper().isin(onc_set)
        cancer_df = df[mask].copy()
        # Select columns to write
        out_df = cancer_df[['g1', 'g2', 'combined_score']].rename(columns={'g1': 'protein_a', 'g2': 'protein_b'})
        out_df.to_csv(outpath, index=False)
        return outpath

    def save_validation_results(self, rows, outpath):
        """Save a list of validation dicts to CSV."""
        df = pd.DataFrame(rows)
        df.to_csv(outpath, index=False)
        return outpath

class PPIApp:
    def __init__(self, root):
        self.root = root
        root.title("Cancer PPI Analysis (STRING + OncoKB) - Extended")
        root.geometry("900x740")
        self.core = PPINetworkAnalysis()
        self._build_widgets()

    def set_busy(self, busy: bool):
        try:
            if hasattr(self, 'run_buttons') and self.run_buttons:
                state = 'disabled' if busy else 'normal'
                for b in self.run_buttons:
                    try:
                        b.config(state=state)
                    except Exception:
                        pass
            if hasattr(self, 'results_text') and self.results_text:
                try:
                    if busy:
                        self.results_text.insert('end', '\n== BUSY ==\n')
                        self.results_text.see('end')
                    else:
                        self.results_text.insert('end', '\n== DONE ==\n')
                        self.results_text.see('end')
                except Exception:
                    pass
        except Exception:
            pass

    def _build_widgets(self):
        frm_top = ttk.Frame(self.root, padding=8)
        frm_top.pack(fill='x')

        ttk.Label(frm_top, text="STRING Links (file):").grid(row=0, column=0, sticky='w')
        self.string_entry = ttk.Entry(frm_top, width=70)
        self.string_entry.grid(row=0, column=1, sticky='w')
        ttk.Button(frm_top, text="Browse", command=self.browse_string).grid(row=0, column=2)

        ttk.Label(frm_top, text="STRING Info (protein.info) (optional):").grid(row=1, column=0, sticky='w')
        self.info_entry = ttk.Entry(frm_top, width=70)
        self.info_entry.grid(row=1, column=1, sticky='w')
        ttk.Button(frm_top, text="Browse", command=self.browse_info).grid(row=1, column=2)

        ttk.Label(frm_top, text="OncoKB TSV:").grid(row=2, column=0, sticky='w')
        self.onco_entry = ttk.Entry(frm_top, width=70)
        self.onco_entry.grid(row=2, column=1, sticky='w')
        ttk.Button(frm_top, text="Browse", command=self.browse_onco).grid(row=2, column=2)

        ttk.Label(frm_top, text="COSMIC TSV (optional):").grid(row=4, column=0, sticky='w')
        self.cosmic_entry = ttk.Entry(frm_top, width=70)
        self.cosmic_entry.grid(row=4, column=1, sticky='w')
        ttk.Button(frm_top, text="Browse", command=self.browse_cosmic).grid(row=4, column=2)

        ttk.Label(frm_top, text="TCGA Gene List TSV (optional):").grid(row=5, column=0, sticky='w')
        self.tcga_entry = ttk.Entry(frm_top, width=70)
        self.tcga_entry.grid(row=5, column=1, sticky='w')
        ttk.Button(frm_top, text="Browse", command=self.browse_tcga).grid(row=5, column=2)

        ttk.Label(frm_top, text="BioGRID (optional):").grid(row=3, column=0, sticky='w')
        self.bgrid_entry = ttk.Entry(frm_top, width=70)
        self.bgrid_entry.grid(row=3, column=1, sticky='w')
        ttk.Button(frm_top, text="Browse", command=self.browse_biogrid).grid(row=3, column=2)

        frm_params = ttk.LabelFrame(self.root, text="Parameters", padding=8)
        frm_params.pack(fill='x', padx=8, pady=6)

        ttk.Label(frm_params, text="STRING score threshold (0-1000):").grid(row=0, column=0, sticky='w')
        self.score_var = tk.IntVar(value=700)
        ttk.Entry(frm_params, textvariable=self.score_var, width=8).grid(row=0, column=1, sticky='w')

        ttk.Label(frm_params, text="Cancer subnetwork hop (neighbors):").grid(row=0, column=2, sticky='w')
        self.hop_var = tk.IntVar(value=1)
        ttk.Entry(frm_params, textvariable=self.hop_var, width=6).grid(row=0, column=3, sticky='w')

        self.include_bg_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(frm_params, text="Include BioGRID edges", variable=self.include_bg_var).grid(row=0, column=4,
                                                                                                  sticky='w')

        ttk.Label(frm_params, text="Louvain resolution:").grid(row=1, column=0, sticky='w')
        self.louvain_var = tk.DoubleVar(value=1.0)
        ttk.Entry(frm_params, textvariable=self.louvain_var, width=8).grid(row=1, column=1, sticky='w')

        ttk.Label(frm_params, text="Sample fraction (0.01=1%, 1.0=all):").grid(row=1, column=2, sticky='w')
        self.sample_var = tk.DoubleVar(value=0.01)
        ttk.Entry(frm_params, textvariable=self.sample_var, width=8).grid(row=1, column=3, sticky='w')

        frm_actions = ttk.Frame(self.root, padding=8)
        frm_actions.pack(fill='x')

        btns = []
        btns.append(ttk.Button(frm_actions, text="Label Propagation", command=self.core.detect_label_propagation))
        btns[-1].grid(row=2, column=0, padx=6)
        btns.append(ttk.Button(frm_actions, text="Girvan–Newman", command=self.core.detect_girvan_newman))
        btns[-1].grid(row=2, column=1, padx=6)
        btns.append(ttk.Button(frm_actions, text="COSMIC Modules", command=self.core.detect_cosmic_seeded_modules))
        btns[-1].grid(row=2, column=2, padx=6)
        btns.append(ttk.Button(frm_actions, text="Load files & Construct Graph", command=self.load_and_construct))
        btns[-1].grid(row=0, column=0, padx=6)
        btns.append(ttk.Button(frm_actions, text="Compute Metrics", command=self.compute_metrics))
        btns[-1].grid(row=0, column=1, padx=6)
        btns.append(ttk.Button(frm_actions, text="Detect Communities (Louvain)", command=self.detect_communities))
        btns[-1].grid(row=0, column=2, padx=6)
        btns.append(ttk.Button(frm_actions, text="Build Cancer Subnetwork", command=self.build_cancer_subgraph))
        btns[-1].grid(row=0, column=3, padx=6)
        btns.append(ttk.Button(frm_actions, text="Export Cancer PPI CSV", command=self.export_cancer_ppi))
        btns[-1].grid(row=0, column=4, padx=6)
        btns.append(ttk.Button(frm_actions, text="Validate Hubs (PubMed+KEGG)", command=self.validate_hubs))
        btns[-1].grid(row=0, column=5, padx=6)
        btns.append(ttk.Button(frm_actions, text="Global Plot (Matplotlib)", command=self.plot_global_network))
        btns[-1].grid(row=1, column=0, padx=6, pady=8)
        btns.append(ttk.Button(frm_actions, text="Degree Highlight Plot", command=self.plot_degree_highlight))
        btns[-1].grid(row=1, column=1, padx=6, pady=8)
        btns.append(ttk.Button(frm_actions, text="Community-Colored Plot", command=self.plot_community_colored))
        btns[-1].grid(row=1, column=2, padx=6, pady=8)
        btns.append(ttk.Button(frm_actions, text="Interactive PyVis (HTML)", command=self.gen_interactive_and_open))
        btns[-1].grid(row=1, column=3, padx=6, pady=8)
        btns.append(ttk.Button(frm_actions, text="Export GEXF for Cytoscape", command=self.export_gexf))
        btns[-1].grid(row=1, column=4, padx=6, pady=8)
        btns.append(ttk.Button(frm_actions, text="Save Metrics CSV", command=self.save_metrics))
        btns[-1].grid(row=1, column=5, padx=6, pady=8)

        self.run_buttons = btns

        frm_bottom = ttk.LabelFrame(self.root, text="Top Hubs (Results)", padding=8)
        frm_bottom.pack(fill='both', expand=True, padx=8, pady=6)

        self.results_text = tk.Text(frm_bottom, wrap='none', height=18)
        self.results_text.pack(side='left', fill='both', expand=True)
        scrollbar = ttk.Scrollbar(frm_bottom, orient='vertical', command=self.results_text.yview)
        scrollbar.pack(side='right', fill='y')
        self.results_text.config(yscrollcommand=scrollbar.set)

        frm_pub = ttk.LabelFrame(self.root, text="PubMed Evidence (optional)", padding=8)
        frm_pub.pack(fill='x', padx=8, pady=6)
        ttk.Label(frm_pub, text="Email (for NCBI Entrez):").grid(row=0, column=0, sticky='w')
        self.email_var = tk.StringVar(value="")
        ttk.Entry(frm_pub, textvariable=self.email_var, width=40).grid(row=0, column=1, sticky='w')
        ttk.Label(frm_pub, text="Gene symbol (single):").grid(row=0, column=2, sticky='w')
        self.pub_gene_var = tk.StringVar(value="")
        ttk.Entry(frm_pub, textvariable=self.pub_gene_var, width=14).grid(row=0, column=3, sticky='w')
        ttk.Button(frm_pub, text="Query PubMed for Gene", command=self.query_pubmed_for_gene).grid(row=0, column=4, padx=6)

    def browse_cosmic(self):
        p = filedialog.askopenfilename(title="Select COSMIC TSV", filetypes=[("TSV files", "*.tsv;*.txt;*.gz"), ("All files", "*.*")])
        if p:
            self.cosmic_entry.delete(0, tk.END)
            self.cosmic_entry.insert(0, p)

    def browse_tcga(self):
        p = filedialog.askopenfilename(title="Select TCGA gene list TSV", filetypes=[("TSV files", "*.tsv;*.txt;*.gz"), ("All files", "*.*")])
        if p:
            self.tcga_entry.delete(0, tk.END)
            self.tcga_entry.insert(0, p)

    def browse_string(self):
        p = filedialog.askopenfilename(title="Select STRING links file", filetypes=[("All files", "*.*")])
        if p:
            self.string_entry.delete(0, tk.END)
            self.string_entry.insert(0, p)

    def browse_info(self):
        p = filedialog.askopenfilename(title="Select STRING protein.info", filetypes=[("All files", "*.*")])
        if p:
            self.info_entry.delete(0, tk.END)
            self.info_entry.insert(0, p)

    def browse_onco(self):
        p = filedialog.askopenfilename(title="Select OncoKB TSV", filetypes=[("TSV files", "*.tsv;*.txt;*.gz"),
                                                                               ("All files", "*.*")])
        if p:
            self.onco_entry.delete(0, tk.END)
            self.onco_entry.insert(0, p)

    def browse_biogrid(self):
        p = filedialog.askopenfilename(title="Select BioGRID file (optional)", filetypes=[("All files", "*.*")])
        if p:
            self.bgrid_entry.delete(0, tk.END)
            self.bgrid_entry.insert(0, p)

    @run_in_thread
    def load_and_construct(self):
        try:
            links_p = self.string_entry.get().strip()
            info_p = self.info_entry.get().strip() or None
            onco_p = self.onco_entry.get().strip()
            bgrid_p = self.bgrid_entry.get().strip() or None
            if not links_p or not onco_p:
                messagebox.showerror("Error", "Please provide at least STRING links file and OncoKB TSV.")
                return

            self.core.params['score_threshold'] = int(self.score_var.get())
            self.core.params['louvain_resolution'] = float(self.louvain_var.get())
            sample_val = float(self.sample_var.get())
            self.core.params['sample_fraction'] = sample_val if 0 < sample_val < 1.0 else None

            print(f"\n{'=' * 60}")
            print(f"Loading with settings:")
            print(f"  Score threshold: {self.core.params['score_threshold']}")
            print(f"  Sample fraction: {sample_val * 100:.1f}%")
            print(f"{'=' * 60}\n")

            self.core.load_string(links_p, info_p)
            self.core.load_oncokb(onco_p)

            cosmic_p = self.cosmic_entry.get().strip() if hasattr(self, 'cosmic_entry') else None
            tcga_p = self.tcga_entry.get().strip() if hasattr(self, 'tcga_entry') else None

            if cosmic_p:
                self.core.load_cosmic(cosmic_p)

            if tcga_p:
                self.core.load_tcga(tcga_p)

            if bgrid_p:
                self.core.load_biogrid(bgrid_p)
                self.core.params['include_biogrid'] = True

            self.core.construct_graph()
            messagebox.showinfo("Done",
                                f"Graph constructed: {self.core.G.number_of_nodes()} nodes, {self.core.G.number_of_edges()} edges")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def compute_metrics(self):
        try:
            self.core.params['score_threshold'] = int(self.score_var.get())
            # compute_metrics runs in background thread; it returns a Thread object
            thread = self.core.compute_metrics()
            # We can poll for completion if desired; here we just inform user to wait in UI
            messagebox.showinfo("Computing", "Metrics computation started in background. Use the UI when done.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def detect_communities(self):
        try:
            res = float(self.louvain_var.get())
            # call background louvain
            thread = threading.Thread(target=lambda: self._detect_communities_worker(res), daemon=True)
            thread.start()
            messagebox.showinfo("Started", "Community detection started in background.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _detect_communities_worker(self, res):
        try:
            partition = self.core.detect_communities(resolution=res)
            # partition is returned from thread function, but our decorator wraps to thread; attempt to call directly if available
            # To ensure UI feedback, compute counts if partition available
            if isinstance(partition, dict):
                counts = Counter(partition.values())
                s = "Detected communities (id: size)\n"
                for k, v in sorted(counts.items(), key=lambda x: -x[1])[:20]:
                    s += f"{k}: {v}\n"
            else:
                s = "Community detection dispatched."
            messagebox.showinfo("Communities", s)
        except Exception as e:
            messagebox.showerror("Error (Community)", str(e))

    @run_in_thread
    def build_cancer_subgraph(self):
        try:
            hop = int(self.hop_var.get())
            sub = self.core.extract_cancer_subgraph(hop=hop, min_degree=0)
            messagebox.showinfo("Done", f"Cancer subgraph built: {sub.number_of_nodes()} nodes, {sub.number_of_edges()} edges")
            deg = dict(sub.degree())
            pr = nx.pagerank(sub, weight='weight') if sub.number_of_nodes() > 0 else {}
            rows = []
            for n in sub.nodes():
                rows.append((n, deg.get(n, 0), pr.get(n, 0.0), n.upper() in self.core.oncokb))
            df = pd.DataFrame(rows, columns=['node', 'degree', 'pagerank', 'is_oncokb']).sort_values('degree', ascending=False)
            txt = ["Top 30 nodes in cancer subgraph\n"]
            for i, r in df.head(30).iterrows():
                onc = " (OncoKB)" if r['is_oncokb'] else ""
                txt.append(f"{r['node']}    degree={r['degree']}    pagerank={r['pagerank']:.4f}{onc}\n")
            self.results_text.delete('1.0', tk.END)
            self.results_text.insert(tk.END, ''.join(txt))
        except Exception as e:
            messagebox.showerror("Error", str(e))

    @run_in_thread
    def plot_static_subgraph(self):
        try:
            sub = self.core.G_sub if self.core.G_sub is not None else self.core.G
            if sub is None or sub.number_of_nodes() == 0:
                messagebox.showerror("Error", "No subgraph available. Build cancer subgraph first or construct graph.")
                return
            n_nodes = sub.number_of_nodes()
            if n_nodes > 300:
                degs = dict(sub.degree())
                top_nodes = sorted(degs.keys(), key=lambda k: degs[k], reverse=True)[:300]
                viz = sub.subgraph(top_nodes).copy()
                print(f"Visualizing top 300 of {n_nodes} nodes")
            else:
                viz = sub
            print("Computing layout (this may take a moment)...")
            plt.figure(figsize=(12, 10))
            pos = nx.spring_layout(viz, k=0.15, iterations=50, seed=42)
            node_colors = ['red' if (n.upper() in self.core.oncokb) else 'skyblue' for n in viz.nodes()]
            node_sizes = [50 + 10 * viz.degree(n) for n in viz.nodes()]
            nx.draw_networkx_nodes(viz, pos, node_size=node_sizes, node_color=node_colors, alpha=0.9)
            nx.draw_networkx_edges(viz, pos, alpha=0.3, width=0.5)
            labels = {n: n for n in viz.nodes() if (n.upper() in self.core.oncokb) or (viz.degree(n) >= 10)}
            nx.draw_networkx_labels(viz, pos, labels=labels, font_size=8)
            plt.title("Cancer-centered PPI subnetwork (static)")
            plt.axis('off')
            plt.tight_layout()
            plt.show()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    @run_in_thread
    def gen_interactive_and_open(self):
        try:
            sub = self.core.G_sub if self.core.G_sub is not None else self.core.G
            if sub is None or sub.number_of_nodes() == 0:
                messagebox.showerror("Error", "No subgraph available to visualize. Build cancer subgraph or construct graph first.")
                return
            path = filedialog.asksaveasfilename(defaultextension=".html", initialfile="ppi_interactive.html",
                                                filetypes=[("HTML", "*.html")])
            if not path:
                return
            self.core.generate_pyvis(graph=sub, path=path, show_in_browser=True, highlight_oncokb=True)
            messagebox.showinfo("Done", f"Interactive network saved to {path} and opened in browser.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def export_gexf(self):
        try:
            path = filedialog.asksaveasfilename(defaultextension=".gexf", initialfile="ppi_export.gexf",
                                                filetypes=[("GEXF", "*.gexf")])
            if not path:
                return
            sub = self.core.G_sub if self.core.G_sub is not None else self.core.G
            self.core.export_gexf(graph=sub, path=path)
            messagebox.showinfo("Done", f"Exported graph to {path}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def save_metrics(self):
        try:
            if self.core.metrics_df is None:
                messagebox.showerror("Error", "No metrics computed. Run 'Compute Metrics' first.")
                return
            path = filedialog.asksaveasfilename(defaultextension=".csv", initialfile="ppi_metrics.csv",
                                                filetypes=[("CSV", "*.csv")])
            if not path:
                return
            self.core.metrics_df.to_csv(path, index=False)
            messagebox.showinfo("Done", f"Saved metrics to {path}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def query_pubmed_for_gene(self):
        try:
            gene = self.pub_gene_var.get().strip()
            email = self.email_var.get().strip()
            if not gene or not email:
                messagebox.showerror("Error", "Please provide both gene symbol and email for Entrez queries.")
                return

            def worker():
                try:
                    results = self.core.query_pubmed(gene_symbol=gene, email=email, retmax=8)
                    if not results:
                        out = f"No PubMed results found for {gene} (query limited to cancer-related titles/abstracts).\n"
                    else:
                        out = f"PubMed results for {gene} (top {len(results)}):\n"
                        for pmid, title, url in results:
                            out += f"{pmid}: {title}\n{url}\n\n"
                    self.results_text.delete('1.0', tk.END)
                    self.results_text.insert(tk.END, out)
                except Exception as ex:
                    messagebox.showerror("Error (PubMed)", str(ex))

            thread = threading.Thread(target=worker, daemon=True)
            thread.start()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def export_cancer_ppi(self):
        try:
            if self.core.string_links is None:
                messagebox.showerror("Error", "STRING links not loaded. Load and construct first.")
                return
            path = filedialog.asksaveasfilename(defaultextension=".csv", initialfile="cancer_ppi.csv",
                                                filetypes=[("CSV", "*.csv")])
            if not path:
                return
            out = self.core.save_cancer_specific_ppi(path)
            messagebox.showinfo("Done", f"Exported cancer-specific PPI to {out}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def validate_hubs(self):
        """Gather top hubs and validate each with PubMed and KEGG; export a validation CSV."""
        try:
            if self.core.metrics_df is None:
                messagebox.showerror("Error", "No metrics available. Run 'Compute Metrics' first.")
                return
            # Ask for basic params via simple dialogs
            topn = 50
            topndlg = filedialog.asksaveasfilename(defaultextension=".csv", initialfile="validation_sheet.csv",
                                                   filetypes=[("CSV", "*.csv")])
            if not topndlg:
                return
            save_path = topndlg
            email = self.email_var.get().strip()
            if not email and Entrez is not None:
                # prompt user to input email if Biopython available and no email provided
                messagebox.showinfo("Note", "Provide an email in the PubMed box if you want PubMed validation.")
            # Gather top hubs
            hubs = self.core.top_hubs(topn=topn, by='degree')
            genes = hubs['node'].tolist()

            # Run validation in background
            thread = threading.Thread(target=lambda: self._validate_hubs_worker(genes, email, save_path), daemon=True)
            thread.start()
            messagebox.showinfo("Started", f"Validation running in background; results will be saved to {save_path} when done.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _validate_hubs_worker(self, genes, email, save_path):
        rows = []
        for g in genes:
            item = {'gene': g, 'is_oncokb': (g.upper() in self.core.oncokb)}
            # PubMed
            pubmed_hits = []
            try:
                if Entrez is not None and email:
                    res = self.core.query_pubmed(gene_symbol=g, email=email, retmax=5)
                    pubmed_hits = [{'pmid': pmid, 'title': title, 'url': url} for pmid, title, url in res]
                else:
                    pubmed_hits = []
            except Exception:
                pubmed_hits = []
            item['pubmed_count'] = len(pubmed_hits)
            # store sample pmids/titles as JSON string (safe for CSV)
            item['pubmed_hits'] = json.dumps(pubmed_hits)

            # KEGG
            kegg_paths = []
            try:
                kegg_paths = self.core.kegg_gene_pathways(g)
            except Exception:
                kegg_paths = []
            item['kegg_pathway_count'] = len(kegg_paths)
            item['kegg_pathways'] = json.dumps(kegg_paths)

            rows.append(item)
            # be polite to remote APIs
            time.sleep(0.2)
        try:
            self.core.save_validation_results(rows, save_path)
            messagebox.showinfo("Validation done", f"Validation completed and saved to {save_path}")
        except Exception as e:
            messagebox.showerror("Error saving validation", str(e))

    @run_in_thread
    def plot_global_network(self):
        """Plot a global network (or top-N nodes) with degree-based sizes and save PNG."""
        try:
            G = self.core.G
            if G is None:
                messagebox.showerror("Error", "No global graph constructed yet.")
                return
            n_nodes = G.number_of_nodes()
            if n_nodes > 800:
                # reduce to top 800 by degree for visual clarity
                degs = dict(G.degree())
                top_nodes = sorted(degs.keys(), key=lambda k: degs[k], reverse=True)[:800]
                viz = G.subgraph(top_nodes).copy()
                print(f"Visualizing top 800 of {n_nodes} nodes")
            else:
                viz = G
            print("Computing layout (global)...")
            plt.figure(figsize=(14, 10))
            pos = nx.spring_layout(viz, k=0.12, iterations=50, seed=42)
            degrees = dict(viz.degree())
            node_sizes = [30 + degrees[n] * 3 for n in viz.nodes()]
            nx.draw_networkx_nodes(viz, pos, node_size=node_sizes, node_color='skyblue', alpha=0.8)
            nx.draw_networkx_edges(viz, pos, alpha=0.2, width=0.4)
            # label only top degree nodes
            top_labels = {n: n for n in sorted(degrees.keys(), key=lambda x: degrees[x], reverse=True)[:40]}
            nx.draw_networkx_labels(viz, pos, labels=top_labels, font_size=8)
            plt.title("Global PPI Network (top nodes labeled)")
            plt.axis('off')
            plt.tight_layout()
            # save file
            p = filedialog.asksaveasfilename(defaultextension=".png", initialfile="global_network.png",
                                             filetypes=[("PNG", "*.png")])
            if p:
                plt.savefig(p, dpi=200)
                messagebox.showinfo("Saved", f"Global network plot saved to {p}")
            plt.show()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    @run_in_thread
    def plot_degree_highlight(self):
        """Plot network with nodes colored by degree bin (low/medium/high)."""
        try:
            G = self.core.G_sub if self.core.G_sub is not None else self.core.G
            if G is None:
                messagebox.showerror("Error", "Construct graph or build subgraph first.")
                return
            degrees = dict(G.degree())
            # degree bins quantiles
            deg_vals = list(degrees.values())
            if not deg_vals:
                messagebox.showerror("Error", "Graph empty.")
                return
            q1 = pd.Series(deg_vals).quantile(0.33)
            q2 = pd.Series(deg_vals).quantile(0.66)
            colors = []
            for n in G.nodes():
                d = degrees.get(n, 0)
                if d <= q1:
                    colors.append('lightgrey')
                elif d <= q2:
                    colors.append('orange')
                else:
                    colors.append('red')
            sizes = [20 + degrees.get(n, 0) * 4 for n in G.nodes()]
            pos = nx.spring_layout(G, k=0.15, iterations=50, seed=42)
            plt.figure(figsize=(12, 10))
            nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=sizes, alpha=0.9)
            nx.draw_networkx_edges(G, pos, alpha=0.2)
            # label high-degree nodes
            labels = {n: n for n, dd in degrees.items() if dd >= q2}
            nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
            plt.title("Network colored by degree (low / med / high)")
            plt.axis('off')
            plt.tight_layout()
            p = filedialog.asksaveasfilename(defaultextension=".png", initialfile="degree_highlight.png",
                                             filetypes=[("PNG", "*.png")])
            if p:
                plt.savefig(p, dpi=200)
                messagebox.showinfo("Saved", f"Degree-highlight plot saved to {p}")
            plt.show()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    @run_in_thread
    def plot_community_colored(self):
        """Plot graph with nodes colored by detected community (Louvain)."""
        try:
            G = self.core.G_sub if self.core.G_sub is not None else self.core.G
            if G is None:
                messagebox.showerror("Error", "Construct graph or build subgraph first.")
                return
            # Ensure communities exist
            if not self.core.community:
                # attempt to detect communities quickly
                try:
                    self.core.detect_communities(graph=G, resolution=self.core.params.get('louvain_resolution', 1.0))
                except Exception:
                    pass
            part = None
            if isinstance(self.core.community, dict):
                part = self.core.community
            else:
                # attempt to read node attribute
                part = nx.get_node_attributes(G, 'louvain_community')
            if not part:
                messagebox.showerror("Error", "No community partition available. Run 'Detect Communities' first.")
                return
            # Map community id to color
            comms = list(set(part.values()))
            cmap = plt.get_cmap('tab20')
            color_map = {c: cmap(i % 20) for i, c in enumerate(sorted(comms))}
            node_colors = [color_map.get(part.get(n, -1), (0.7, 0.7, 0.7)) for n in G.nodes()]
            pos = nx.spring_layout(G, k=0.12, iterations=50, seed=42)
            plt.figure(figsize=(12, 10))
            nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=40, alpha=0.9)
            nx.draw_networkx_edges(G, pos, alpha=0.15)
            # label large communities' top nodes
            degree = dict(G.degree())
            labels = {}
            for c in comms:
                members = [n for n in G.nodes() if part.get(n) == c]
                if not members:
                    continue
                # top node by degree
                top = sorted(members, key=lambda x: degree.get(x, 0), reverse=True)[:2]
                for t in top:
                    labels[t] = t
            nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
            plt.title("Community-colored network (Louvain)")
            plt.axis('off')
            plt.tight_layout()
            p = filedialog.asksaveasfilename(defaultextension=".png", initialfile="community_colored.png",
                                             filetypes=[("PNG", "*.png")])
            if p:
                plt.savefig(p, dpi=200)
                messagebox.showinfo("Saved", f"Community-colored plot saved to {p}")
            plt.show()
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == '__main__':
    root = tk.Tk()
    app = PPIApp(root)
    root.mainloop()