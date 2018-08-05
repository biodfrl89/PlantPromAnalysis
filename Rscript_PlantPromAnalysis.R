#!/usr/bin/env Rscript
library(optparse)

option_list <-  list(
  make_option( opt_str = c("-f", "--file"), type = "character", default = NULL, 
               help = "The fasta file with the sequence(s) to be analized.", metavar = "filename"),
  make_option( opt_str = c("-o", "--out"), type = "character", default = "out.txt", 
               help = "The filename of the output [default = %default]", metavar = "filename"),
  make_option( opt_str = c("-db", "--database"), type = "character", default = "place", 
               help = "A string defining the database that is going to be used to make the analisys [default = %default]. It can be modified to 'atcis'", 
               metavar = "String")
  ) 

opt_parser <-  OptionParser(option_list = option_list)
opt <-  parse_args(opt_parser)

if (is.null(opt$file)) {
  stop("At least one argument must be supplied (input file).\n\n Look the help section with the argument -h or --help \n\n", call.=FALSE)
}

library(seqinr)

test="test"

#Se carga la base de datos PLACE y se suprime la primer columna que tiene un dato desconocido sin sentido
db <- read.table(file = "DATABASES/place_raw.txt", header = FALSE)
db$V1 <- NULL
#Se nombran las columnas
colnames(db) <- c("Name", "Motif", "Length", "ID")
#Se lee la secuencia a analizar y se cambia a mayusculas para que no haya problema
promoter <- read.fasta(file = opt$file, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
promoter <- read.fasta(file = "atlea45.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
promoter <- toupper(promoter)

#ASIGNACION DE BASES DEGENERADAS
{
  deg_bases <- list(W = c("A", "T"),
                    S = c("C", "G"),
                    M = c("A", "C"),
                    K = c("G", "T"),
                    R = c("A", "G"),
                    Y = c("C", "T"),
                    B = c("C", "G", "T"),
                    D = c("A", "G", "T"),
                    H = c("A", "C", "T"),
                    V = c("A", "C", "G"),
                    N = c("A", "C", "G", "T")
  )
}

#PREPARACIÓN DEL MOTIVO A BUSCAR
#Se declaran las variables vacias a utilizar
pos_motifs <- list(0)
#Primer ciclo: se va a ciclar a traves de todos los elementos de la columna motif de la db
j <- 1
x <- 1
for (j in 1:length(db$Motif)) {
	#El elemento en cuestion se almacena en una variable temporal como un vector
	motif <- as.vector(db$Motif[j])
	#Segundo ciclo: se hace un ciclaje en donde se detecta bases degeneradas y se sustituyen por las bases
	#ATCG formando al final una expresión regular que puede ser evaluada
	for (x in 1:11) {
		motif <- gsub(names(deg_bases[x]), paste(c("[", deg_bases[[x]], "]"), collapse = ""), motif)
	}
	x <- 1
	#Con la expresion regular, se obtiene la posicion de dicha expresion regular en la secuencia a analizar
	result <- list(words.pos(motif, promoter))
	#El resultado se guarda como un elemento de una lista
	pos_motifs[j] <- result
	#Se repite el ciclo y se continuan guardando los resultado encontrados
}

#A cada elemento de la lista se le asigna el nombre del elemento cis que le corresponde
names(pos_motifs) <- db$Name
#Hay elementos cis que no se encuentran y se graban como character(0). Con esta funcion se eliminan
cis_founded <- Filter(length, pos_motifs)
result_tb <- plyr::ldply(cis_founded, rbind)
result_tb[is.na(result_tb)] <- ""
print(result_tb)


if (test == "test") stop("Test done.")



db_filtrada <- db[db$Name %in% names(cis_elements_founded),]


df_cis <-reshape2::melt(cis_elements_founded)
df_cis <- df_cis[c(2,1)]
df_cis_ordered <- df_cis[order(df_cis$value),]
df_cis_ordered$y <- rep(0, length(df_cis_ordered$value))
write.table(unique(df_cis$L1), file = "res_unicos_atlea61.txt", quote = FALSE, 
						row.names = FALSE, col.names = FALSE)

#Para convertir el archivo de resultados en expresiones regulares
#sed -e 's/^/\^/' res_unicos_atlea61.txt >inicio.txt 

#Para preparar la base de datos con la información y que quede lista para su analisis
#cat place.seq | grep '^ID\|^DE\|^KW\|^SQ\|^//' -B1 | grep '^XX\|^--\|^SQ' -v | sed 's/^  /SQ/g'  
	#>char_DE_descriptors.txt
#Se hace lo mismo para los elementos deseados

#Para univer los archivos generados y dejar la base de datos prepara donde se realizara el grep.
#paste -d '\n' char_ID_geneid.txt char_DE_descriptors.txt char_KW_keywords.txt char_SQ_sequence.txt | 
	#sed '0~4 s/$/\n/g' >place_seq_final.txt

#Se hace un grep de las secuencias unicas encontradas sobre la base de datos.
#grep -A3 -wf inicio.txt place_seq_final.txt 





