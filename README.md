# RG-RGG-proteins
Data and code used for bioinformatical analysis of RG(G) proteins, presented in Figure 1A and Figure 2A of the main text and Supplementary Figure 5. 

## Occurrence of RG and RGG repeats in the human IDR-ome

In the Fig 1A and Fig 2B repository, the find_rg_repeats.py is used to generate the occurrence_matrix_RG.txt and occurrence_matrix_RGG.txt files in the output repository. The script needs to be called from the Fig 1A and Fig 2A directory as:
```
find_rg_repeats.py /path/to/fasta/file motif
```
where motif is either RG or RGG (find_rg_repeats.py /path/to/fasta/file RG    or    find_rg_repeats.py /path/to/fasta/file RGG).

The input .fasta file used to generate the output data is given in the same repository (UP000005640_9606_SPOTD_MIN_30AA.fasta).

The output files are in the form of a matrix where each row corresponds to a different IDR and columns correspond to the number of times each of the RG.{0,5}RG (+.{0,5}RG)*11 were found in each of the IDRs. In other words, for each IDR we are counting the number of times two to twelve consecutive RG.{0,5} are found in their sequence. 

The data was further analyzed in Excel.

## Frequency of each amino acid in RG regions

In the SuppFig 5, the RG_molecular_grammar.py is used to generate the X_composition.txt files in the output repository. The analysis was done on three .fasta files: 
- UP000005640_9606_SPOTD_MIN_30AA.fasta  ~  human IDRome
- stress_granules_seq.fasta  ~  sequences of proteins known to localize in stress granules (taken from [Youn et al. 2018, MolCell][https://pubmed.ncbi.nlm.nih.gov/29395067/])
- TNPO1_binding_seq.fasta  ~  sequences of proteins know to bind to TNPO1 (taken from [Mackmull et al. 2017, MolSystBiol][https://pubmed.ncbi.nlm.nih.gov/29254951/] and [Kimura et al. 2018, eLife][https://elifesciences.org/articles/21184])

The RG_molecular_grammar.py script needs to be called from the SuppFig 5 directory separtely for each .fasta file:
```
RG_molecular_grammar.py /path/to/fasta/file
```
as:
```
RG_molecular_grammar.py UP000005640_9606_SPOTD_MIN_30AA.fasta
RG_molecular_grammar.py Stress_granules_seq.fasta
RG_molecular_grammar.py TNPO1_binding_seq.fasta
```

The script counts the number of times each amino acid appears in RG regions of all of the sequences given in a .fasta file. The RG region is defined as everything between the first arginine of the current RG motif and glycine of the last RG in the motif, plus 10 amino acids before and after region. The script iteratively searches for RG motifs of 3, 4 and 5 consecutive RG repeats: RG.{0,5}RG.{0,5}RG, RG.{0,5}RG.{0,5}RG.{0,5}RG.{0,5} and RG.{0,5}RG.{0,5}RG.{0,5}RG.{0,5}RG. 

For each .fasta file we give as input, we get a RG_3_x_composition.txt, RG_4_x_composition.txt and RG_5_x_composition.txt files where a global frequency (number of times each amino acid occurs in the current RG region across all IDRs from the .fasta file divided by the total number of amino acids in RG regions). These frequencies are then compared to the amino acid composition in human IDRome based on the null hypothesis, which are given in IDRome_AA_COMPOSITION.txt file. 

The final step of the analysis was done in Excel. The expected frequencies from the null hypothesis are subtracted from the calculated frequencies in the RG regions of different sets of IDRs (IDRome, stress granules localizing, TNPO1 binding) and the log odds analysis was performed on the difference.


[https://pubmed.ncbi.nlm.nih.gov/29395067/]: https://pubmed.ncbi.nlm.nih.gov/29395067/
[https://pubmed.ncbi.nlm.nih.gov/29254951/]: https://pubmed.ncbi.nlm.nih.gov/29254951/
[https://elifesciences.org/articles/21184]: https://elifesciences.org/articles/21184