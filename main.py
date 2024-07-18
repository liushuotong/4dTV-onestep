import subprocess
from io import StringIO
from Bio import SeqIO, AlignIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import os
import uuid

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

def run_muscle(input_sequences):
    # define unique id
    unique_id = str(uuid.uuid4())
    input_file = f"input_sequences_{unique_id}.fasta"
    output_file = f"aligned_sequences_{unique_id}.fasta"

    # Write input sequences to temporary file
    with open(input_file, "w") as f:
        f.write(input_sequences)

    # Run muscle
    command = ["muscle.exe", "-in", input_file, "-out", output_file]

    # Suppress output
    with open(os.devnull, 'w') as devnull:
        subprocess.run(command, stdout=devnull, stderr=devnull)

    # Read aligned sequences
    with open(output_file, "r") as f:
        aligned_sequences = f.read()
    
    # Remove temporary files
    os.remove(input_file)
    os.remove(output_file)
    
    return aligned_sequences

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

# Process each gene pair
def process_gene_pair(args):
    gene1, gene2, all_sequences = args

    # Create in-memory FASTA file content
    temp_fasta_content = StringIO()
    SeqIO.write([all_sequences[gene1], all_sequences[gene2]], temp_fasta_content, "fasta")
    temp_fasta_content.seek(0)

    # Run Muscle for alignment
    aligned_fasta_content = run_muscle(temp_fasta_content.read())
    
    # 读取对齐结果
    aligned_fasta_content = StringIO(aligned_fasta_content)
    try:
        alignment = AlignIO.read(aligned_fasta_content, "fasta")
    except ValueError:
        print(f"Error reading alignment for {gene1} and {gene2}")
        return (gene1, gene2, None)

    # Extract four-fold degenerate sites
    four_fold_sites = extract_four_fold_sites(alignment)

    # Calculate transversions and 4dTV value
    transversions = calculate_transversions(four_fold_sites)
    total_sites = len(four_fold_sites)
    if total_sites > 0:
        four_dtv = transversions / total_sites
        return (gene1, gene2, four_dtv)
    else:
        return (gene1, gene2, None)

# Main function
def main(tsv_file, fasta_file, output_file):
    homologs = read_homologs(tsv_file)
    all_sequences = extract_sequences(fasta_file, set(homologs["Gene1"]).union(set(homologs["Gene2"])))
    results = []

    # Prepare arguments for multiprocessing
    args = [(row["Gene1"], row["Gene2"], all_sequences) for _, row in homologs.iterrows()]

    with Pool(cpu_count()) as pool:
        for result in tqdm(pool.imap_unordered(process_gene_pair, args), total=len(args)):
            results.append(result)
    
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
    # Extracting 4dTV values
    values = [value for _, _, value in results if value is not None]
    
    # Check if there are any values to plot
    if len(values) == 0:
        print("No 4dTV values to visualize.")
        return
    
    # Plot the KDE plot using seaborn
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.kdeplot(values, bw_adjust=0.5, fill=True)
    
    # Set plot labels and title
    plt.xlabel('4dTV')
    plt.ylabel('Density')
    plt.title('Distribution of 4dTV Values for Gene Pairs')
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    tsv_file = input("PATH to homologs.tsv: ")
    fasta_file = input("PATH to sequences.fasta: ")
    output_file = input("PATH to output_file.tsv: ")
    main(tsv_file, fasta_file, output_file)
