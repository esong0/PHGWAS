# Author: Euijun Song
# Python >= 3.8

import os
import networkx as nx
import numpy as np
from ripser import ripser
from persim import plot_diagrams
import pickle
import matplotlib.pyplot as plt

class PHGWAS:
    def __init__(self, file_PPI="interactome.pkl", file_GWAS="gwashit_t2dm_LCC.tsv", use_LCC=True):
        with open(file_PPI, 'rb') as f:
            self.G = pickle.load(f)
        
        seeds_, pvals_, mlogpvals_ = [], [], []
        for line in open(file_GWAS, 'r'):
            if line[0] == '#':
                continue
            line_data = line.strip().split('\t')
            if line_data[0] in self.G.nodes:
                seeds_.append(line_data[0]) # Gene
                pvals_.append(float(line_data[2])) # P-value
                mlogpvals_.append(float(line_data[3])) # -log(P-value)
        print("# of imported seeds in PPI:", len(seeds_))
        
        self.seeds, self.pvals, self.mlogpvals = zip(*sorted(zip(seeds_, pvals_, mlogpvals_), key=lambda x: x[1]))
        self.seeds, self.pvals, self.mlogpvals = list(self.seeds), list(self.pvals), list(self.mlogpvals)
        
        self.G = self.G.subgraph(self.seeds)
        if use_LCC:
            self.G = self.G.subgraph(max(nx.connected_components(self.G), key=len)).copy()
            print("Observable disease module:", len(self.G.nodes))
    
    def homology_compute(self, max_dim=2, thr_type="pval"):
        A = np.full((len(self.seeds), len(self.seeds)), np.inf)
        np.fill_diagonal(A, 0)
        gene_to_idx = {g:i for i, g in enumerate(self.seeds)}

        # assuming self.pvals is sorted.
        for i, gene in enumerate(self.seeds):
            if thr_type == "pval":
                delta = self.pvals[i]
            elif thr_type == "logpval":
                delta = -self.mlogpvals[i]
            else:
                raise Exception("Threshold type error.")
            A[i][i] = delta
            for gene_adj in self.G[gene]:
                j = gene_to_idx[gene_adj]
                A[i][j] = delta
                A[j][i] = delta
        
        perst_dgms = ripser(A, distance_matrix=True, maxdim=max_dim)
        diagrams = perst_dgms['dgms']
        print("# of PD:", [len(dgm) for dgm in diagrams])
        for i, dgm in enumerate(diagrams):
            print("------ %d-th persistent homology ------" % i)
            print(dgm)
        return diagrams
    
    def persistence_diagram(self, diagrams, file_plot="persistence_diagram.pdf"):
        non_empty_pd = []
        for i, dgm in enumerate(diagrams):
            if len(dgm) > 0:
                non_empty_pd.append(i)
        
        plt.figure(figsize=(5,4))
        plot_diagrams(diagrams, plot_only=non_empty_pd, colormap='seaborn')
        plt.xlabel('Birth')
        plt.ylabel('Death')
        plt.tight_layout()
        plt.savefig(file_plot)
        plt.show()
    
    def h1_cycle_compute(self):
        h1_cycles = [c for c in nx.minimum_cycle_basis(self.G) if len(c)>3]
        print("# of H1 cycles:", len(h1_cycles))
        h1_cycle_union = set()
        for c in h1_cycles:
            h1_cycle_union = h1_cycle_union | set(c)
        print("Union of H1 cycles:", sorted(h1_cycle_union))
        print("Size:", len(h1_cycle_union))
        return h1_cycle_union
