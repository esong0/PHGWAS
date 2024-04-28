# Author: Euijun Song
# Python >= 3.8


import PHGWAS as phg
import gseapy as gp
from gseapy import barplot, dotplot
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.weight'] = "bold"


if __name__ == "__main__":
    ph_gwas = phg.PHGWAS()
    diagrams = ph_gwas.homology_compute()
    h1_dz_module = ph_gwas.h1_cycle_compute()
    h1_dz_module = sorted(list(h1_dz_module))

    enr = gp.enrichr(gene_list=h1_dz_module,
                    gene_sets=['KEGG_2021_Human'],
                    organism='human',
                    outdir='functional_analysis_result',
                    )
    ax = dotplot(enr.results,
                column="Adjusted P-value",
                group='Gene_set',
                top_term=10,
                figsize=(3,5),
                cmap=plt.cm.viridis,
                ofname='functional_analysis_result/plot_KEGG.pdf'
                )

    enr = gp.enrichr(gene_list=h1_dz_module,
                    gene_sets=['ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X'],
                    organism='human',
                    outdir='functional_analysis_result',
                    )
    ax = dotplot(enr.results,
                column="Adjusted P-value",
                group='Gene_set',
                top_term=10,
                figsize=(3,5),
                cmap=plt.cm.viridis,
                ofname='functional_analysis_result/plot_TF.pdf'
                )

    enr = gp.enrichr(gene_list=h1_dz_module,
                    gene_sets=['miRTarBase_2017'],
                    organism='human',
                    outdir='functional_analysis_result',
                    )
    ax = dotplot(enr.results,
                column="Adjusted P-value",
                group='Gene_set',
                top_term=10,
                figsize=(3,3),
                cmap=plt.cm.viridis,
                ofname='functional_analysis_result/plot_miRNA.pdf'
                )
