# Molecular mimicry identification

Scripts for the recreation of the pipeline described in [Identification of potential molecular mimicry in pathogen-host interactions](https://doi.org/10.1101/2023.06.14.544818). All commands and parameters used are detailed in [pipeline_commands.sh](https://github.com/Kayleerich/molecularmimicry/blob/main/pipeline_commands.sh). Supplemental data files can be found in [supp_data](https://github.com/Kayleerich/molecularmimicry/tree/main/supp_data), including comprehensive Gene Ontology enrichment results, the full images referred to in Figure 6 and TSV files containing candidate mimicry sequences/sequence coordinates.

The pipeline requires *Python 3* with `pandas` and `numpy` as well as pre-installation of [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [Phobius](https://phobius.sbc.su.se/), [bedtools](https://github.com/arq5x/bedtools2), and [POPScomp](https://github.com/Fraternalilab/POPScomp). 

Ensure that sequences in FASTA files have unique names. Sequence names **should not contain spaces, `\` or `:`**. Uniprot identifiers work well (for example: `>sp|P04004|VTNC_HUMAN Vitronectin` to `>P04004`), and are easily renamed using `bioawk`: 

	bioawk -c fastx '{split($name, a, "|"); print ">"a[2]"\n"$seq}' sequences.fasta
