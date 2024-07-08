import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import matplotlib.pyplot as plt

# Define four-fold degenerate codons
four_fold_degenerate_codons = {
    'GGT', 'GGC', 'GGA', 'GGG',  # Gly
    'CCT', 'CCC', 'CCA', 'CCG',  # Pro
    'ACT', 'ACC', 'ACA', 'ACG',  # Thr
    'GCT', 'GCC', 'GCA', 'GCG',  # Ala
    'TCT', 'TCC', 'TCA', 'TCG',  # Ser
    'CGT', 'CGC', 'CGA', 'CGG',  # Arg
    'GTT', 'GTC', 'GTA', 'GTG'   # Val
}

# Read homolog gene pairs from TSV file
def read_homologs(tsv_file):
    return pd.read_csv(tsv_file, sep='\t')

# Extract gene sequences from FASTA file
def extract_sequences(fasta_file, gene_list):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in gene_list:
            sequences[record.id] = record
    return sequences

# Perform multiple sequence alignment using Muscle
def run_muscle(input_fasta, output_fasta):
    muscle_executable = "./muscle"  # Assuming muscle is in the same directory
    subprocess.run([muscle_executable, "-in", input_fasta, "-out", output_fasta])

# Extract four-fold degenerate sites
def extract_four_fold_sites(alignment):
    sites = []
    for i in range(0, alignment.get_alignment_length(), 3):
        codon1 = str(alignment[0, i:i+3].seq)
        codon2 = str(alignment[1, i:i+3].seq)
        if codon1[:2] == codon2[:2] and codon1 in four_fold_degenerate_codons and codon2 in four_fold_degenerate_codons:
            sites.append((codon1[2], codon2[2]))
    return sites

# Calculate transversions
def calculate_transversions(sites):
    transversions = 0
    for base1, base2 in sites:
        if base1 != base2 and (base1 in "AG" and base2 in "CT" or base1 in "CT" and base2 in "AG"):
            transversions += 1
    return transversions

# Main function
def main(tsv_file, fasta_file, output_file):
    homologs = read_homologs(tsv_file)
    all_sequences = extract_sequences(fasta_file, set(homologs["Gene1"]).union(set(homologs["Gene2"])))
    results = []

    for index, row in homologs.iterrows():
        gene1 = row["Gene1"]
        gene2 = row["Gene2"]

        # Create temporary FASTA file
        with open("temp.fasta", "w") as temp_fasta:
            SeqIO.write([all_sequences[gene1], all_sequences[gene2]], temp_fasta, "fasta")

        # Run Muscle for alignment
        run_muscle("temp.fasta", "aligned.fasta")

        # Read alignment result
        alignment = AlignIO.read("aligned.fasta", "fasta")

        # Extract four-fold degenerate sites
        four_fold_sites = extract_four_fold_sites(alignment)

        # Calculate transversions and 4dTV value
        transversions = calculate_transversions(four_fold_sites)
        total_sites = len(four_fold_sites)
        if total_sites > 0:
            four_dtv = transversions / total_sites
            results.append((gene1, gene2, four_dtv))
            print(f"{gene1} - {gene2} 4dTV: {four_dtv}")
        else:
            results.append((gene1, gene2, None))
            print(f"{gene1} - {gene2} 4dTV: No four-fold degenerate sites found")

    # Save results to TSV file
    save_results_to_tsv(results, output_file)

    # Visualize 4dTV values
    visualize_results(results)

# Save results to TSV file
def save_results_to_tsv(results, output_file):
    df = pd.DataFrame(results, columns=["Gene1", "Gene2", "4dTV"])
    df.to_csv(output_file, sep='\t', index=False)

# Visualize results
def visualize_results(results):
    fig, ax = plt.subplots()
    genes = [f"{gene1}-{gene2}" for gene1, gene2, _ in results if _ is not None]
    values = [value for _, _, value in results if value is not None]

    ax.barh(genes, values)
    ax.set_xlabel('4dTV')
    ax.set_title('4dTV Values for Gene Pairs')
    plt.show()

if __name__ == "__main__":
    tsv_file = input("PATH to homologs.tsv: ")
    fasta_file = input("PATH to sequences.fasta: ")
    output_file = input("PATH to output_file.tsv: ")
    main(tsv_file, fasta_file, output_file)
