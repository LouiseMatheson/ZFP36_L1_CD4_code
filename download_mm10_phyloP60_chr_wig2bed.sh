#!/bin/bash
# requires bedops module to be loaded

wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/phyloP60way/mm10.60way.phyloP60way/${1}.phyloP60way.wigFix.gz

gunzip ${1}.phyloP60way.wigFix.gz
wig2bed < ${1}.phyloP60way.wigFix > ${1}.phyloP60way.bed
gzip ${1}.phyloP60way.bed

rm ${1}.phyloP60way.wigFix
