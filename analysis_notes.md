# *C. glabrata* MEC isolate analysis

## MEC001 assembly
### 2023-12-18
Self-alignment of long reads to manually corrected fasta still shows uneven
coverage at most chromosome ends. RDNA arrays on L and M are likely to be a lot
more than just 4 arrays as in the CBS138 long read assembly. Looking at the
canu assembly as well to see if I can find a better version of the RDNA and
flanking genes. Without bubbles there are 28 contigs and the extras all seem to
be slightly different versions of the RDNA repeats. Maybe it's reflecting real
variation in the arrangement of the genes if they get shuffled during
replication? Also, there do appear to be ORFs for TNR(2 - low complexity repeat
gene families and strain differences are going to be a problem in this species)
and EPA13 and 14 (notably a and b are identical in the CBS138 long read
assembly, was this clear in the paper?) but there are mnay stops in the RDNAs.
I don't know how many introns may be present or what's unique to the different
subunits. Stops also present in the reference genome and in SC5314 A21, so this
is not just misassembly problems. The subunits all blast between 98 and 100%
identity to CBS138 so should annotate correctly when that step is reached.

## MEC001 assembly
### 2023-12-15
Reviewing MEC001 assembly. Telomeres are present on all but chr G of assembly
(and right arms of L and M but present on the collapsed rDNA contig). Will use
blast to identify rDNA and ITS structure/spacing and manually stitch 4 repeats
together to apply to clean-cut of right side of both L and M. Will use
identical manual correction for each because EPA14a and EPA14b both appear to
be identical in CBS138 (old and long read?).

## MEC001 assembly
### 2023-12-11
Blast individual ribosomal genes in contig 8 as well as flanking genes for both
ChrL and ChrM. It looks like contig 8 is the rDNA array for both of them,
collapsed (flanking genes which are parts of repeat-rich gene families in both
chrs both blast with identity in upper 80s to low 90s).

## TE discussion with Ursula and Anna
### 2023-12-06
Ursula has analyzed the *C. glabrata* reference genome (CGD s05m03r02) and the
PacBio fastq data for the 5 MEC isolates. There appears to be only a single TE
(family?) in all isolates - TCG3, which is a Ty element - and it may actually
be active. There are 6 TCG3 copies in CBS138 and variable numbers in the MEC
isolates.

## Centromere MSA
### 2023-11-07
IGV views and MAFFT MSAs don't really agree - can see much more variation in
IGV - BUT this is coming from soft-clipping in BWA-MEM alignments. To avoid
soft-clipping, will use BBMAP for global alignment of reads. BBMAP global
aligment of short reads and minimap2 alignment of long reads (MEC001) show no
major deletions or even really any "complex" variants, just a somewhat limited
number of SNPs. Will re-run centromere subsetting, consensus sequence
generation and MSA using the BBMAP BAM files.

## Centromere MSA
### 2023-11-03
Generated new VCF file for all MEC samples, with only simple quality filtering
(mapping quality and supporting reads); subset to centromere regions only
(vcf_cens.sh). Remade consensus sequences (array_cens_subsetting.sh,
concat_cens.sh) and re-ran MAFFT (mafft_cglabrata.sh). Used Table S1 from Dujon
et al 2004 (Genome evolution in yeast, original C. glabrata CBS138 assembly)
for CDE information. CDEII region has greater variation than CDEI or the reported
essential regions of CDEIII, but still remains AT-rich.

## Centromere MSA
### 2023-11-02
*C. glabrata* has point centomeres similar to *S. cerivisiae* with 3 distinct
regions - CDEI, II and III. Since megablast doesn't work in each full centomere
region as a whole, compared each centromere sequence to check the known
regions.
Generated consensus sequences from all MEC isolates and CBS138 with bcftools, and
concatenated each centromere to a multi-sample fasta. Ran MAFFT v7.475 with
recommended settings from Ursula and viewed output in Jalview 2.11.2.7. Several
deletions appear to be missing compared to IGV view. "Master" VCF file
filtering is too strict for centromere regions (complex variants).

## De novo assembly
### 2023-09-29
Built local blast database of MEC001 flye assembly for centromere search.
Generated consensus sequences from CGD reference genome + VCF file and blasted
against assembly. Changed which centromeres were found. Determined that
megablast will not identify all centromeres between our strains and CBS138
because identity drops to ~80% in some cases. Blastn identifies all
centromeres.

## De novo assembly
### 2023-09-28
Built local database of CGD CBS138 assembly for blasting the 27k contig
(contig_8 in original flye assembly) to identify areas of similarity. Aligns to
~8kb segment of ChrE and ~17kb segment of ChrM. Using minimap2 to align the
flye assembly to CBS138 and generating dot plot and genome coverage plot using
pafR shows alignment to ChrM.
Re-ran flye on the full fastq sequencing set, but set parameters to do the
initial disjointig with ~60x coverage downsampled from 400x, as flye will use
only the longest reads to accomplish this. Again, like the filtlong subsetting,
the final output is same number of contigs, with similar lengths, but poorer
graph paths in the longer contigs - again missing terminal end nodes and having
a path through multiple contigs instead of a single. So far the full coverage
gives best first-pass output. As of version 2.9, flye does not scaffold. So
another option is to try with the scaffolding parameter on, as well as trying
to patch or merge using ragtag and the CGD reference genome. Still need to pull
the read IDs corresponding to centromeres in CBS138, in order to understand
what's happening to them in assembly.

## De novo assembly
### 2023-09-27
Re-ran flye assembly with filtlong-subset fastq. Same number of final contigs,
but the actual repeat graphs have deteriorated. The longest contigs in the
subset assembly don't all have 2 terminal end nodes and the graph path includes
multiple contigs, rather than a single. The 27kb fragment that appears to be
some kind of repeat is still present.
Aligned pacbio reads to the CGD CBS138 reference genome to compare to Illumina
data. Overall agreement, and so far still haven't identified where the 27kb
fragment corresponds to the ref genome - is it an adhesin with increased copy
number? Is it part of ChrE vs ChrM?

## De novo assembly
### 2023-09-26
Downloaded HiFiAdapterFilt, shell script, for checking raw reads for
contamination. No adapters or contamination identified in MEC001. Install
filtlong to subset MEC001 fastq to a minimum of 6kb read length and downsample
to 1 million bases  (from 5 million) in hopes of improving flye assembly.

## De novo assembly
### 2023-09-23
Latest long read CBS138 ref genome for *C. glabrata*:
Downloaded fasta and GFF from CGD, indexed and re-aligned MEC isolates (in
2022_Cglabrata/align_variants_cgd_cbs138).

## De novo assembly
### 2023-09-22
Re-ran ragout on MEC001 using flye and hifiasm assemblies as draft references
for canu assembly. July error occurs because no contigs or scaffolds are placed
- everything ends up in the unplaced contigs fasta.
Will work directly on flye assemblies (comparing raw and pilon-polished) -
starting with building local blast database of the assemblies to identify where
the centromeres are.

## Fix clustering
### 2023-08-25
Clustering scripts were filtering for indels instead of SNPs (altered script
template in early 2022 when comparing SNPs and indels, never changed it back,
Christopher pointed out error), so re-ran to subset genotypes tables to only
SNPs, and repeat MCA and plotting.

## De novo assembly
### 2023-07-13
Getting ragout failure without final scaffold output, log files
suggest a missing intermediate scaffold input.

## De novo assembly
### 2023-07-12
Hifiasm raw assemblies look bad. Parameter for haploid/homozygous organisms
doesn't seem to work as expected?
For ragout input, try flye, masurca, canu and 1 or 2 hybrid assemblies? Or just
install all tools and run the Gabaldon tool interactively.

## De novo assembly
### 2023-07-11
Flye assembly looks best so far in comparing unpolished results with between 14
adn 17 contigs per assembly, with all but 1 having a circular mitochondrial
scaffold. Need to try polishing with Illumina. Multiple sources seem to do ~3
rounds of polishing. Need to incorporate contaminant detection and removal
prior to annotation.
Consider using ragout with multiple draft assemblies used as input for each
tool's raw assembly target? Like the Gabaldon lab's approach (LongHam...).
Consider using same approaches for polishing and cleaning up *C. lusitaniae*
assembly.

## MLST typing
### 2023-06-26
Get locus names for all *Candida* species with MLST typing schemes. Build
instrument to include ST, allele IDs of exact matches or best match data
(allele ID, % identity, number of mismatches) and isolate sequences used to
generate allele information. Finish python script to include REDCap API call to
upload dict in JSON format (pubmlst_rest.py) and then run array script for both
*C. albicans* and *C. glabrata* to process all MEC isolates (mlst_api.sh).

## MLST typing
### 2023-06-23
Write a python script to perform single-isolate query of *C. glabrata* sequence
types. Create a dict from API response values including ST and locus exact
matches. Loop over missing loci to find best matches.

## MLST typing
### 2023-06-22
Generate locus fasta sequences for *C. glabrata*: upload CBS138 reference fasta
to Pubmlst.org single sequence query, download the results as a fasta sequence
to get locus coordinates in the genome. Get fasta regions using samtools faidx
with -r option and regions txt file (format is chr:from-to, one per line),
piped to bcftools (local module for version 1.17) consensus, with -I option for
IUPAC symbols and -s sample for individual sample genotype. Compared bcftools
consensus output and IGV consensus sequence output, and results from sequence
query at Pubmlst. Everything matches.


## Genome visualization (local ymap-ish)
### 2023-05-19
Remove mosdepth results (3/27/23) in favor of samtools + local R script. Move
gc-corrected bams to global scratch. Run local genome visualization (scripts:
"candida_ymap.sh" and "genome_vis.R"). Modify genome_vis.R to include SNP
density (by summing positions with AF >= 0.9 and != ref position over the
specified window used for copy number) instead of LOH; using the two-color
scale that works for LOH does not work for SNP density. Even excluding NA
values, ggplot geom_density does not recognize SNP sum as a continuous value.
Therefore, removed any kind of SNP plotting from all *C. glabrata* isolates and
re-ran as copy number only.

## Clustering
### 2023-05-17
Review clustering results and zoom in to specific clusters using
coord_cartesian, including in dims 1-2 and dims 2-3. Also colored samples by
"demographic" clusters. No obvious groupings in the demographic clusters.
Scanned one SNP-based cluster in IGV with an outlier (MEC 96, 358, 264, 265,
281 and outlier 364). Need to assign groupings by MLST  typing and review
clusters again. Approximately 7 clusters of isolates by SNPs - does this just track
with MLST/clade? or are clades dispersed among clusters? What are typical
SNP/kb values within and between MLST groups?

## Clustering
### 2023-05-02
Run cglabrata_mca.R to plot SNP-based clustering on genotypes table (unfixed
SNPs relative to CBS138). Need to write a function to add samples of interest
for zooming and labeling.

## Copy number by read depth
### 2023-03-27
Perform GC-correction and mosdepth (2000bp window) on all bam files for use in
copy-number plotting (script: "gc_correct_mosdepth.sh").

## MEC 328 series review
### 2023-03-10
CHEF gel of MEC328-MEC335 series shows 3 different karyotypes. Reviewed
spreadsheets in box to confirm MRN is same for all samples. Subset VCF file to
only series samples, filtered to remove reference-only and fixed SNP sites
(need to doublecheck this bcftools syntax). Visualized in IGV. All isolates
appear clonal. *PDR1* ortholog has fixed SNPs, plus one unique missense SNP per
karyotype.

## Variant annotation and filtering
### 2022-12-05
Built snpEff database with standard nuclear codon table.
Concatenated chr vcfs and filtered (see script for parameters).
Annotated filtered vcf, further filtering with snpsift to only predicted high
and moderate impact variants.

## Variant calling
### 2022-11-16
Chr 12 timed out (12 hours).
Restarted with 22 hours wall time.

## Variant calling
### 2022-11-15
"Population" variant calling: performed variant calling for all samples with freebayes, as
array for each chromosome.
Freebayes parameters: -C 10, -F 0.9, -p 1

## QC
### 2022-11-12
Basic QC: ran fastqc, qualimap bamqc and multiqc for all samples.
All bbduk trimmed fastq, fastqc, bamqc files kept on global scratch.
Bam files, samtools flagstats and multiqc output kept in this project
directory.

## Alignment to CBS138 ASM254v2 reference genome
### 2022-11-11
Alignment array failed for five samples (out of memory): MEC136-MEC140.
All are part of a series from one patient (full series is MEC136-MEC141).
MEC141 aligned successfully, so ran centrifuge classifier on MEC and MEC141
using ncbi nt index.
Both classified as *C. glabrata*, with ~97% reads for MEC141 and ~95% reads for
MEC.
Re-ran bwa mem alignment to CBS138 for samples MEC136-140 with increased memory
allocation.

## Pre-processing, alignment to CBS138 ASM254v2 reference genome
### 2022-11-09
Generated sample sheet of space-delimited sample ID, path/to/fastq/R1 and
path/to/fastq/R2, named Cglab_sequencing_paths.txt.
Performed adapter and quality trimming with bbduk.
Aligned with bwa-mem to reference genome CBS138.
