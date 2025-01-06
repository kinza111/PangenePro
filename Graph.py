import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from collections import defaultdict, Counter
import os

def load_orthologs(file_path):
    data = pd.read_csv(file_path, sep='\t', header=None, names=['Genome1_Gene', 'Genome2_Gene', 'Similarity'])
    data[['Genome1', 'GeneID1']] = data['Genome1_Gene'].str.split('|', expand=True)
    data[['Genome2', 'GeneID2']] = data['Genome2_Gene'].str.split('|', expand=True)
    return data

def classify_genes(data):
    gene_presence = defaultdict(set)
    for _, row in data.iterrows():
        gene_presence[row['Genome1']].add(row['GeneID1'])
        gene_presence[row['Genome2']].add(row['GeneID2'])

    core_genes = set.intersection(*gene_presence.values())
    accessory_genes = (set.union(*gene_presence.values()) - core_genes) - set.union(*(gene_presence[g] - core_genes for g in gene_presence))
    unique_genes = {genome: genes - set.union(*(gene_presence[g] for g in gene_presence if g != genome)) for genome, genes in gene_presence.items()}

    return core_genes, accessory_genes, unique_genes, gene_presence

def generate_output_table(core_genes, accessory_genes, unique_genes, gene_presence, output_file):
    rows = []
    for gene in core_genes:
        rows.append({"Gene": gene, "Category": "Core", "Genomes": list(gene_presence.keys())})
    
    for gene in accessory_genes:
        present_genomes = [genome for genome, genes in gene_presence.items() if gene in genes]
        rows.append({"Gene": gene, "Category": "Accessory", "Genomes": present_genomes})
    
    for genome, genes in unique_genes.items():
        for gene in genes:
            rows.append({"Gene": gene, "Category": "Unique", "Genomes": [genome]})
    
    output_df = pd.DataFrame(rows)
    output_df.to_csv(output_file, index=False)
    print(f"Pangenes sets Table saved to {output_file}.")

def save_bar_plot(core_genes, accessory_genes, unique_genes):
    plt.figure(figsize=(8, 6))
    plot_data = {
        'Core Genes': len(core_genes),
        'Accessory Genes': len(accessory_genes),
        'Unique Genes': sum(len(genes) for genes in unique_genes.values())
    }
    plt.bar(plot_data.keys(), plot_data.values(), color=['#1f77b4', '#ff7f0e', '#2ca02c'])
    plt.title("Core, Accessory, and Unique Genes")
    plt.ylabel("Gene Count")
    plt.savefig("bar_plot.svg")
    plt.close()
    print("Bar plot saved as bar_plot.svg.")

def save_upset_plot(gene_presence):
    gene_sets = []
    for gene, presences in pd.DataFrame([{genome: gene in gene_presence[genome] for genome in gene_presence}
                                        for gene in set.union(*gene_presence.values())]).iterrows():
        gene_sets.append(tuple(presences.values))

    gene_set_counts = Counter(gene_sets)
    genomes = list(gene_presence.keys())

    fig, (ax_bar, ax_matrix) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(14, 8))
    fig.subplots_adjust(hspace=0.05)

    x_labels = []
    bar_heights = []
    for presence_comb, count in gene_set_counts.items():
        x_labels.append(" & ".join([genomes[i] for i, present in enumerate(presence_comb) if present]))
        bar_heights.append(count)
    ax_bar.bar(range(len(bar_heights)), bar_heights, color='#4c72b0')
    ax_bar.set_xticks([])
    ax_bar.set_ylabel("Gene Count")

    for i, presence_comb in enumerate(gene_set_counts):
        for j, is_present in enumerate(presence_comb):
            color = '#4c72b0' if is_present else '#d9d9d9'
            ax_matrix.plot(i, len(genomes) - j - 1, 'o', color=color)

    ax_matrix.set_yticks(range(len(genomes)))
    ax_matrix.set_yticklabels(genomes, fontsize=8)
    ax_matrix.set_xticks(range(len(bar_heights)))
    ax_matrix.set_xticklabels(x_labels, rotation=90, fontsize=8)
    ax_matrix.set_ylabel("Genomes")

    plt.suptitle("Custom UpSet-style Plot of Gene Overlaps")
    plt.savefig("custom_upset_plot.svg")
    plt.close()
    print("Custom UpSet-style plot saved as custom_upset_plot.svg.")

def save_venn_diagram(core_genes, gene_presence):
    if len(gene_presence) == 3:
        genome_keys = list(gene_presence.keys())
        plt.figure(figsize=(8, 8))
        venn3(
            subsets=(
                len(gene_presence[genome_keys[0]] - gene_presence[genome_keys[1]] - gene_presence[genome_keys[2]]),
                len(gene_presence[genome_keys[1]] - gene_presence[genome_keys[0]] - gene_presence[genome_keys[2]]),
                len(gene_presence[genome_keys[0]].intersection(gene_presence[genome_keys[1]]) - gene_presence[genome_keys[2]]),
                len(gene_presence[genome_keys[2]] - gene_presence[genome_keys[0]] - gene_presence[genome_keys[1]]),
                len(gene_presence[genome_keys[0]].intersection(gene_presence[genome_keys[2]]) - gene_presence[genome_keys[1]]),
                len(gene_presence[genome_keys[1]].intersection(gene_presence[genome_keys[2]]) - gene_presence[genome_keys[0]]),
                len(core_genes)
            ),
            set_labels=genome_keys
        )
        plt.title("Venn Diagram of Gene Overlaps")
        plt.savefig("venn_diagram.svg")
        plt.close()
        print("Venn diagram saved as venn_diagram.svg.")

def main():
    orthologs_file = "/output/orthagogue/Orthologs.abc"
    output_table = "gene_sets_summary.csv"

    data = load_orthologs(orthologs_file)
    core_genes, accessory_genes, unique_genes, gene_presence = classify_genes(data)

    generate_output_table(core_genes, accessory_genes, unique_genes, gene_presence, output_table)
    save_bar_plot(core_genes, accessory_genes, unique_genes)
    save_upset_plot(gene_presence)
    save_venn_diagram(core_genes, gene_presence)

# Execute the main function
if __name__ == "__main__":
    main()
