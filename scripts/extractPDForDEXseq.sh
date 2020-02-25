#!/bin/bash

for d in `cat ./FigAndResults/DEanalysis/allPathDEXseq.txt`; do 
    echo $d
    for pdf in "$d"/BRCA1.pdf
    do
        xt=`echo ${d} | sed -e "s/\//_/g"`
        cp ${pdf} /home/lemrap/serveurs/crihan/local_data/RNASeq/QC_RNAseq/QC-RNAseqAnalysis/FigAndResults/DEanalysis/BRCA1_${xt}.pdf
    done
    for pdf in "$d"/BRCA2.pdf
    do
        xt=`echo ${d} | sed -e "s/\//_/g"`
        cp ${pdf} /home/lemrap/serveurs/crihan/local_data/RNASeq/QC_RNAseq/QC-RNAseqAnalysis/FigAndResults/DEanalysis/BRCA2_${xt}.pdf
    done
done

## merge pdf with https://deftpdf.com/merge-pdf
