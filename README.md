README: Trait-specific chromatin architectures channel pleiotropic gene expression toward sexually dimorphic development in horned beetles.

Abstract:

This readme file was generated on 2026-01-10 by Erica M. Nadolski

GENERAL INFORMATION
Title of Dataset: Trait-specific chromatin architectures channel pleiotropic gene expression toward sexually dimorphic development in horned beetles

1st Author/Corresponding Author Information

Name: Erica M. Nadolski
ORCID: 0000-0003-3314-5305
Institution: Indiana University Bloomington
Address: 915 E 3rd St, Bloomington, IN 47405
Institutional Email: emnadols@iu.edu
Alternative Email: erica.nadolski@gmail.com

2nd Author/Principal Investigator Information

Name: Armin P. Moczek
ORCID: 0000-0002-3478-9949
Institution: Indiana University Bloomington
Address: 915 E 3rd St, Bloomington, IN 47405
Email: armin@iu.edu

Dates of data collection: September-November 2022 (ATAC-seq data); November 2023 (RNA-seq data); May - October 2025 (RNAinterference data)

Information about funding sources that supported the collection of the data: NSF grant IOS-xxx

DATA & FILE OVERVIEW

File List:
ATACseq_read_processing_pipeline.txt
Ot_allpeaks.bed
Ot_peak_counts.txt
Ot_salmon_counts.txt
Ot_sample_info_ATAC.txt
Ot_sample_info_RNA.txt
Otau2_prot_annotations.txt
Otau_sex_diff.R
RNAseq_read_processing_pipeline.txt
geneID2GO.txt
salmon-dsx-quant.tar.gz

Relationship between files:
Otau_sex_diff.R to be used with Ot_peak_counts.txt, Ot_sample_info_ATAC.txt, Ot_sample_info_RNA.txt, Otau2_prot_annotations.txt, geneID2GO.txt, and salmon-dsx-quant.tar.gz

ATACseq_read_processing_pipeline.txt to be used with 
RNAseq_read_processing_pipeline.txt

Additional related data collected that was not included in the current data package: N/A

METHODOLOGICAL INFORMATION

Description of methods used for collection/generation of data: See manuscript.

Methods for processing the data: See manuscript.

Software-specific information needed to interpret the data:
fastqc (v0.11.9)
trimmomatic (v0.36)
salmon (v1.10.1)
python (v3.10.5)
bowtie2 (v2.5.1)
macs (v2.2.9.1)
samtools (v1.17)
bedtools (v2.31.0)
R (v4.2.2)
RStudio (v2022.12.0+353)
BiocManager (v1.30.22)
tximport (v1.34.0)
DESeq2 (v1.38.3)
edgeR (v3.40.2)
RColorBrewer (v1.1-3)
tidyverse (v2.0.0)
ggplot2 (v3.5.2)
pheatmap (v1.0.12)
gplots(v3.2.0)
topGo(v2.58.0)

Standards and calibration information, if appropriate: N/A

Environmental/experimental conditions: See manuscript.

Describe any quality-assurance procedures performed on the data: N/A

People involved with sample collection, processing, analysis and/or submission:

Phillip L. Davidson - sample rearing and collection
Yongsoo Choi - sample rearing
Miranda Towse - sample rearing

DATA-SPECIFIC INFORMATION FOR: Ot_sample_info_RNA.txt

Number of variables: 4

Number of cases/rows: 60

Variable List:
Sample - sample identification code.
Sex - sex of pupa.
Trait - body region origin of dissected epithelial tissue.
Replicate - identification number of dissected pupa.

Missing data codes: N/A

Specialized formats or other abbreviations used: N/A

DATA-SPECIFIC INFORMATION FOR: Ot_sample_info_ATAC.txt

Number of variables: 5

Number of cases/rows: 50

Variable List:
Sample - sample identification code.
Sex - sex of pupa.
Trait - body region origin of dissected epithelial tissue.
Replicate - unique identification number of dissected pupa.
Batch - identification number of sequencing run.

Missing data codes: N/A

Specialized formats or other abbreviations used: N/A

DATA-SPECIFIC INFORMATION FOR: Ot_peak_counts.txt

Number of variables: 54

Number of cases/rows: 175885

Variable List:
chr - scaffold of peak location. 
start - first base pair of peak location. 
end - last base pair of peak location. 
peak - unique identification number of peak.
variables 5:54 -  matrix of peak read counts per sample

Missing data codes: N/A

Specialized formats or other abbreviations used: N/A

DATA-SPECIFIC INFORMATION FOR: Otau2_prot_annotations.txt

Number of variables: 2

Number of cases/rows: 16057

Variable List:
OT_ID - unique gene ID within Otau3.0 genome
Gene.Ontology.IDs - list of all gene ontology terms mapping to gene.

Missing data codes: N/A

Specialized formats or other abbreviations used: N/A

DATA-SPECIFIC INFORMATION FOR: geneID2GO.txt

Number of variables: 5

Number of cases/rows: 16734

Variable List:
gene - unique gene ID within Otau3.0 genome.
OT2_ID - unique gene ID within Otau2.0 genome.
pident - percent shared identity of amino acids in protein.
eval - BLAST expect value, estimated number of random alignments with a particular score or better that could be found by chance in a given database search.
OT2_description - gene annotation from Otau2.0 genome.

Missing data codes: N/A

Specialized formats or other abbreviations used: N/A

DATA-SPECIFIC INFORMATION FOR: xxx.txt

Missing data codes: Missing data represented by NAs and blank cells.

Specialized formats or other abbreviations used: N/A
