## rename FASTA sequences before starting to only include unique uniprot ID or similar alphanumeric name
## BLASTP full-length parasite proteins against control proteins
makeblastdb -in CONTROLS.fasta -dbtype prot -parse_seqids
blastp -db CONTROLS.fasta -query PARASITE.fasta -num_threads 10 -outfmt 6 -out CONTROL_hits.out -evalue 1e-10 
cut -f 1 CONTROL_hits.out | sort -u > CONTROL_hits.list

## identify and mask signal peptides
perl phobius.pl -short PARASITE.fasta > PARASITE_phobius.out 2> PARASITE_phobius.err
python3 make_bed.py -i PARASITE_phobius.out -t phobius
bedtools maskfasta -fi PARASITE.fasta -bed PARASITE_phobius.bed -fo PARASITE_phobius_masked.fasta -mc X

## remove full-length parasite proteins with significant hit to a control sequence
grep '>' PARASITE.fasta | sed 's/>//' | sort -u > PARASITE.list
comm -23 PARASITE.list CONTROL_hits.list > PARASITE_specific_proteins.list
seqtk subseq PARASITE_phobius_masked.fasta PARASITE_specific_proteins.list | fold -w 60 > PARASITE_specific_proteins.fasta

## make k-mers from parasite-specific proteins
python3 make_kmers.py -i PARASITE_specific_proteins.fasta -k 14 -n 0

## BLASTP parasite k-mers against control proteins
blastp -db CONTROLS.fasta -query PARASITE_kmers.fasta -num_threads 10 -outfmt 6 -out CONTROL_hit_kmers.out -evalue 1 -comp_based_stats 0 -word_size 2 -ungapped

## filter parasite-vs-controls k-mer BLASTP results, remove k-mers with high identity to control
python3 filter_blast_control_file.py -i CONTROL_hit_kmers.out | sort -u > CONTROL_hit_kmers.list
grep '>' PARASITE_kmers.fasta | sed 's/>//' | sort -u > PARASITE_kmers.list
comm -23 PARASITE_kmers.list CONTROL_hit_kmers.list > PARASITE_specific_kmers.list
seqtk subseq PARASITE_kmers.fasta PARASITE_specific_kmers.list | fold -w 60 > PARASITE_specific_kmers.fasta

## BLASTP parasite-specific k-mers against host proteins
makeblastdb -in HOST.fasta -dbtype prot -parse_seqids
blastp -db HOST.fasta -query PARASITE_specific_kmers.fasta -num_threads 10 -outfmt 6 -out HOST_hit_kmers.out -evalue 1 -comp_based_stats 0 -word_size 2 -ungapped

## filter parasite-vs-host k-mer BLASTP results
python3 filter_host.py -i HOST_hit_kmers.out > PARASITE_candidate_kmers.out

## count number of hits per identity and k-mers for each (optional)
python3 count_filtered_kmers.py -i PARASITE_candidate_kmers.out > PARSITE_highID_counts.txt

## calculate QSASA for every protein structure, combine output from individual chains 
ls *.pdb.gz > structures_file_names.list
while IFS= read -r file; do
    name=$( echo ${file} | sed 's/^AF-//' | sed 's/.pdb.gz//' )
    ~/POPScomp/POPSC/src/pops --pdb {wdir}/${file}  --zipped --residueOut --popsOut pops_${name}.out 
done < structure_file_names.list

## merge overlapping regions and filter mimicry candidates by mean residue solvent accessibility
python3 region_qsasa_avgs.py -i PARASITE_candidate_kmers.out -p PARASITE -s HOST -k 14 -d path/to/structure/parent/directory -q 0.75

## filter mimicry candidates by LCR content
segmasker -in HOST.fasta -infmt fasta -parse_seqids -out HOST_segmask_intervals.txt
python make_bed.py -i HOST_segmask_intervals.txt -t segmasker
cut -f6,9,10 PARASITE_paired_qsasa_avg075.tsv | sort -k1,1 -k2,2n > HOST_filtered_regions.bed
bedtools intersect -a HOST_filtered_regions.bed -b HOST_segmask_intervals.bed -wo -f 0.5 -v > HOST_LCR_filtered.bed
cut -f 1,2,3,4,5,6,9,10 PARASITE_paired_qsasa_avg075.tsv | grep -f HOST_LCR_filtered.bed | sort -u > LCR_filtered_mimicry_candidates.tsv
