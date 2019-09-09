# Material and methods used for the QC RNA-seq project
---

This repository gathers the scripts and references to analyse the data provide by the [SpliceLauncher pipeline](https://github.com/raphaelleman/SpliceLauncher,"tittle"). The RNA-seq experiment was performed on the behalf of the [ENIGMA consortium](https://enigmaconsortium.org/,"tittle").

![ENIGMA](https://github.com/raphaelleman/QC_RNA-seq_analysis/blob/master/refData/ENIGMA-logo.png)

To use these scripts you will need the dependancies :
+ [DEXseq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html,"tittle")
+ [HTSeq count](https://htseq.readthedocs.io/en/release_0.11.1/,"tittle")
+ [SpliceLauncher](https://github.com/raphaelleman/SpliceLauncher,"tittle") + its dependencies
+ R software with the librairies:
    * ade4
    * stats
    * VennDiagram
    * gplots

to downlod this repository:
```bash
git clone https://github.com/raphaelleman/QC_RNA-seq_analysis.git
```

The participating laboratories are :

| Acronymes | Laboratories, Country | RNA-seq approach |
|:------------:|:--------:|:------------:|
| Logan | Department of Pathology, University of Otago, Christchurch, New Zealand | PCR enrichment RNA-seq |
| Petra | Laboratoř onkogenetiky, Ústav bio chemie a experimentální onkologie, 1. LF UK v Praze | targeted RNA-seq |
| Kvist | Division of Oncology and Pathology, Department of Clinical Sciences, Lund University, Lund, Sweden | targeted RNA-seq |
| Douglas | Human Development and Health Academic Unit, Faculty of Medicine, University of Southampton, Southampton General Hospital, Southampton, UK. | whole transcriptome RNA-seq |
| Vang |     Department of Clinical Genetics, Odense University Hospital, Odense, Denmark. | whole transcriptome RNA-seq |
| Rien | Department of Clinical Genetics, Maastricht University Medical Centre, Maastricht, The Netherlands. | targeted RNAseq |
| Rachid | Ambry genetics, USA | Clone-seq |
| Caen | Centre François Baclesse, France | targeted RNA-seq |
| Inge | ClinicalCancerResearch Center,AalborgUniversityHospital,Aalborg,Denmark | whole transcriptome RNA-seq |


## Read count for gene and exon

To compare the coverage across the different sites, we propose a short pipeline that uses HTSeq count and DEXseq. We provide also in refData, examples of reference GTF file used by these tools for the *BRCA1* and *BRCA2* genes (hg19). To perform this read count, the pipeline starts from the aligned BAM files. It uses samtools to sort the BAM files, then convert them in SAM files. From these SAM files, the pipeline counts the read number of each region defined in the GTF file, examples ar in [refData](https://github.com/raphaelleman/QC_RNA-seq_analysis/blob/master/refData, "tittle"). The pipeline performs the read count differently if the data is single-end or paired-end (option: **-p**).

Usage example:

```bash
bash ./scripts/countReadGeneExon.sh -B /path/to/bamfolder \
    -O /path/to/output/ \
    -p \
    --samtools /path/to/samtools \
    --htseq /path/to/htseq \
    --dexseq /path/to/dexseq script \
    --refGene /path/to/geneCoordonate.gtf \
    --refExon /path/to/exonCoordonate.gtf
```

## Descritpive analysis from read count

We propose two type of unsupervised analyses to compare the different laboratory data: qualitative and quantitative. The qualityative analysis is illustrated by a Venn diagramm. Please note that over 5 ensembles can't be plotted. For the quantitative analyses we used the PCA. Before these analysis, we need to build a pooled matrix of read count, the junction expression is representated by the geometric average of each loboratory.

1 **Merge matrix count**

Use the script *mergeMatrixCount.r*, define the folder were are saved the matrix count (.txt) and the name of output file. By default the script taka only into account the *BRCA1/2* transcripts

Usage example:

```bash
Rscript ./scripts/mergeMatrixCount.r -I /path/to/matrices directory -O /path/to/output file
```

2 **Venn diagramm**

Here we provide the script *getVennDiagram.r* to generate the Venn diagram from the common matrix count. The number of common junction across laboratories independently of number of ensembles (=laboratories). However the diagramm could be draw only if there is less than 5 ensembles (option: **--plot**). The user can also define a list of juncyion of interest in the format *chrXX_start_end_strand*, *e.g.* the known alternative splicing of *BRCA1/2* in the file *listCoordAltEvent.txt* (option: **-L**).

```bash
Rscript getVennDiagram.r -I /path/to/pooled matrix \
        -L /path/to/list of junction file \
        -O /path/to/output file.pdf \
        --plot
```

3 **PCA**

Here we provide the script *getACP.r* to generate the PCA analysis from the common matrix count. The user can also define a list of juncyion of interest in the format *chrXX_start_end_strand*, *e.g.* the known alternative splicing of *BRCA1/2* in the file *listCoordAltEvent.txt* (option: **-L**).

```bash
Rscript getACP.r -I /path/to/pooled matrix \
        -L /path/to/list of junction file \
        -O /path/to/output file.pdf
```


## Comparison of IVS10-2 samples

For the moment, we don't provide a script to merge speciquelly the IVS10-2 samples, due to the heterogenous data. However we propose the script *analyzeIVS10-2data.r* to get the illustration of the expression correlation and we made the distinction between NMD inhibitor treated samples and untreated samples.

```bash
Rscript analyzeIVS10-2data.r \
        -I1 /path/to/untreated data \
        -I2 /path/to/treated data \
        -O /path/to/output folder/
```
