# PlantPromAnalysis

## Description

This package is designed to analyze DNA sequence considered as promoter regions of genes from plants. It uses the published databases from PLACE -- Plant cis-acting regulatory DNA elements¹ and AtcisDB - Arabidopsis cis-regulatory element database². 

## Introduction

The Package includes a Rscript which is the program that must be launched from terminal with the appropiate arguments. 

To make the program run, use this syntax:

	Rscript Rscript_PlantPromAnalysis.R [-f,--file] [-o, --out] [-db, --database] [-h, --help]

The script must take at least the '--file' argument which must be a fasta file. It is recommended to keep the headers of each sequence is as short as possible. Make sure there are no headers duplicated, which may cause the program to overwrite the result of the sequences with the same header.

After running the program, a new directory will be created with the results of the analysis. Each file contains a table with the relative position (closer to 0 mean closer to the gene of interest) of each cis element.  

## References

¹ Higo, K., Y. Ugawa, M. Iwamoto and T. Korenaga (1999) Plant cis-acting regulatory DNA elements (PLACE) database. Nucleic Acids Res. 27 (1): 297-300. [PLACE Database](http://www.dna.affrc.go.jp/htdocs/PLACE/)

² https://agris-knowledgebase.org/AtcisDB/ [AtCIS Database](https://agris-knowledgebase.org/AtcisDB/)
