# *C. glabrata* MEC isolate analysis

## Isolates
There are 98 *C. glabrata* MEC isolates analyzed for this project (see
"Cglabrata_sequencing_paths.txt"). All MEC isolates are single colony
subcultures from blood cultures collected from UMN hospitals
between  2019 and 2021. These were collected under IRB ID STUDY00006473.

## Metadata
Clinical metadata is available in a Box spreadsheet and the UMN data shelter
for study members. De-identified isolate data is available in a Selmecki lab
REDCap  database (Selmecki Candida Isolate Data).

## DNA extraction and Sequencing
DNA extraction was performed with the standard Selmecki lab phenol protocol
(bead beating setting may have varied) in multiple batches between 2021 and
2022.
Illumina sequencing for all isolates was performed at SeqCenter (formerly MiGS)
on a NextSeq 2000 (the 400 Mb sequencing package, paired end 151 bp reads).
Insert sizes range from ~100 bp to > 400bp and median coverage ranges from 34x
to 135x.

## Illumina sequencing
All fastq files are uploaded to MSI (see "Cglabrata_sequencing_paths.txt" for
file locations). All files have been submitted to SRA (NCBI BioProject
PRJNA1069094). See the REDCap database or NCBI submission portal for biosample
and SRA accession numbers.

## Trimming, alignment, QC
Samples were trimmed with BBDuk (BBMap v38.94) and aligned to the CBS138
s05-m03-r02 reference genome with bwa mem v0.7.17, default parameters, followed
by sorting and duplicate marking with samtools v1.10 (see "array_bbduk_align_sort.sh").
Basic stats were generated with samtools flagstats. Trimmed reads and
alignments were QC'd with FASTQC v0.11.9  and Qualimap v2.2.2-dev, and results
compiled with MultiQC v1.13 (see "basic_qc.sh").

## FreeBayes variant calling, annotation, filtering
Variants were called with freebayes v1.1.0, with parameters of -p 1 and -C 10,
using each chromosome as a separate region (see "array_alb_freebayes.sh").
Chromosome concatenation, quality filtering and annotation were performed with
BCFtools v1.17 and SnpEff v5.0e (see "concat_filter_annotate.sh"). The final
set of filtered SNPs and indels for all 98 MEC isolates is
align_variants_cgd_cbs138/bwa/Cglabrata_MEC_bwa_filtered_annotated.vcf.gz.

## Clustering
Variant sites with biallelic SNPs and no missing genotypes were used to
generate a matrix for clustering ("concat_filter_annotate.sh"). Multiple
correspondence was performed using the FactoMiner package v2.9 (see
"cglabrata_clustering.sh" and "cglabrata_mca.R"). R data was saved as
align_variants_cgd_cbs138/bwa/clustering/Cglabrata_SNP_mca.rda for
reuse/additional plotting.

## Multilocus sequence typing
"Consensus sequences" of the 6 MLST loci were generated for all 98 MEC
isolates (see "array_mlst_subsetting.sh"). The PubMLST database was queried via API to
assign allele IDs and sequence types and to upload to REDCap (see "mlst_api.sh"
and "pubmlst_rest.py"). New alleles were submitted to PubMLST for assignment of
an ID.

## Genome visualization
A local version of YMAP's genome-wide CNV plotting was performed using
Deeptools for GC bias correction and samtools for depth and SNP calling. Allele
frequencies and read depth were checked to confirm haploid status. See
"candida_gc.sh", "candida_ymap.sh", "berman_count_snps_v5.py" and
"genome_vis.R". All plots are uploaded to the Selmecki google drive in
2023_Clinical_MEC_isolates.

## PacBio sequencing
Five isolates (MEC001, MEC175, MEC207, MEC258, MEC339) were submitted to
SeqCenter for HMW DNA extraction and multiplexed PacBio sequencing (PacBio
Sequel IIe for 15hr collection). Reads are in
/home/selmecki/shared/disaster_recovery/Sequencing_Runs/AnnaSelmecki230621 and
have been uploaded to PRJNA1069094. Reads were QC'd with NanoPlot 1.41.6.

## Long read assembly


## Phylogenetic tree building
### MEC isolates only:


### Expanded *C. glabrata* data set:


## Copy number analysis


## Phenotypic data


## Karyotyping (CHEF gels)

