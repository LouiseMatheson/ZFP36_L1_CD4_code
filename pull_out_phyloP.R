library(rtracklayer)
library(tidyverse)

import("chr15.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

consData %>%
  filter(between(start, 96686500, 96700000)) %>% #Slc38a2
  bind_rows({
    consData %>%
      filter(between(start, 61985000,61991000)) #Myc
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 81870000,81920000)) #Aco2
  }) -> selected_data


import("chr1.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 52180000,52240000)) #Gls
  }) -> selected_data

import("chr2.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 11470000,11510000)) #Pfkfb3
  }) -> selected_data

import("chr4.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 119105000,119140000)) #Slc2a1
  }) -> selected_data

import("chr5.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 77080000,77120000)) #Hopx
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 135775000,135795000)) #Mdh2
  }) -> selected_data

import("chr6.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 122727000,122745000)) #Slc2a3
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 82720000,82780000)) #Hk2
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 83305000,83320000)) #Mthfd2
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 72432000,72440000)) #Mat2a
  }) -> selected_data

import("chr8.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 128680000,128740000)) #Itgb1
  }) -> selected_data

import("chr9.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 54580000,54610000)) #Idh3a
  }) -> selected_data

import("chr11.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 100470000,100530000)) #Acly
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 35765000,35795000)) #Pank3
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 97090000,97120000)) #Tbx21
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 6280000,6360000)) #Ogdh
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 120805000,120825000)) #Fasn
  }) -> selected_data

import("chr12.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 73901000,73948000)) #Hif1a
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 76935000,76965000)) #Max
  }) -> selected_data

import("chr13.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 6580000,6660000)) #Pfkp
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 96645000,96675000)) #Hmgcr
  }) -> selected_data

import("chr14.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 34310000,34350000)) #Glud1
  }) %>%
  bind_rows({
    consData %>%
      filter(between(start, 56250000,56270000)) #Gzmb
  }) -> selected_data

import("chr19.phyloP60way.bed.gz") %>%
  as_tibble() -> consData

selected_data %>%
  bind_rows({
    consData %>%
      filter(between(start, 43495000,43525000)) #Got1
  }) -> selected_data


write_tsv(selected_data, "PhyloP_selected_regions_CD4manuscript.txt")
