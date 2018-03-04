# Exploring the Process of Methylation Induced Mutations in Cancer

## Introduction

The mutations observed in a given cancer genome are the result of different mutational processes, whether endogenous or exogenous, acting upon the cell. These processes generate the mutation space that the cell samples in order to select for those mutations that give it an advantage. These procceses leave distinct patterns or signatures in the genome, such as these:

The goal of this project is to characterize one of the most prevalent mutational proccess observed, which is that of ****C:G to T:A transition in a CpG context.**** The process responsible for this signature is thought to be ****methylation induced deamination of cytosine to thymine.****


![](assets/markdown-img-paste-20180304161030985.png)

## Whole Genome Bisulfite Sequencing(WGBS)

Whole genome bisulfite sequencing is considered the gold standard methylation assay. It provides single base resolution of nearly every cytosine in the genome. It is a sequencing based technology that provides as output the a ratio of methylated reads to total reads at a given CpG. The following diagram illustrates how it works:


 The purpose of this tool is to generate the methylation/mutation profiles of ****protein coding genes****. The profiles integrate(at most) 3 types of data:

- WGBS data from x normal tissues
- Somatic Mutation data from x cancer types (obtained from ICGC)
- WGBS data from 3 cancer types(MALY, PBCA, CLLE)

The profiles come in 3 flavors:

- Normal WGBS + Mutation (Pan-Cancer)

- Normal(specific tissue) + Mutations(matched matched cancer type)

- Normal WGBS + Mutation + Cancer WGBS
