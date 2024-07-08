# 4dTV-onestep Calculation Script

This script calculates the four-fold degenerate transversion (4dTV) values for homologous gene pairs. The 4dTV value is a measure of the transversions at four-fold degenerate sites, providing insights into the evolutionary distances between gene pairs.

## Prerequisites

Ensure you have the following installed:
- Python 3.x
- Biopython
- pandas
- matplotlib
- Muscle (place the executable in the same directory as this script)

## Usage

1. **Prepare your data:**
   - A TSV file containing homologous gene pairs with columns "Gene1" and "Gene2".
   - A FASTA file containing the sequences of the genes listed in the TSV file.

2. **Run the script:**

   ```bash
   python main.py
   ```

   You will be prompted to provide the paths to your TSV file, FASTA file, and the desired output file.

## Script Workflow

1. **Read homolog gene pairs:** The script reads the gene pairs from the TSV file.
2. **Extract gene sequences:** It extracts the corresponding gene sequences from the FASTA file.
3. **Multiple sequence alignment:** The script performs multiple sequence alignment using Muscle.
4. **Extract four-fold degenerate sites:** It identifies the four-fold degenerate sites from the aligned sequences.
5. **Calculate transversions:** The script calculates the number of transversions at these sites and computes the 4dTV value.
6. **Save and visualize results:** Finally, it saves the results to a TSV file and generates a bar plot of the 4dTV values.

## Example

### Input

**homologs.tsv**

```
Gene1  Gene2
geneA1 geneA2
geneB1 geneB2
...
```

**sequences.fasta**

```
>geneA1
ATGCGT...
>geneA2
ATGCGT...
>geneB1
ATGCGT...
>geneB2
ATGCGT...
...
```

### Output

**output_file.tsv**

```
Gene1  Gene2  4dTV
geneA1 geneA2 0.123
geneB1 geneB2 0.456
...
```

### Visualization

A bar plot showing the 4dTV values for each gene pair.

## Explanation of 4dTV

- **Four-fold degenerate codons:** Codons where any change in the third position does not alter the amino acid.
- **Transversions:** Substitutions of a purine for a pyrimidine or vice versa (e.g., A ↔ T, G ↔ C).
- **4dTV value:** The ratio of transversions at four-fold degenerate sites to the total number of four-fold degenerate sites.

## Troubleshooting

- Ensure Muscle is executable and in the same directory as the script.
- Verify the FASTA and TSV files are correctly formatted.
- Ensure all genes in the TSV file are present in the FASTA file.

## License

This script is released under the MIT License.

## Contact

For any questions or issues, please contact [liushuotong0218@gmail.com].
