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

#PLACE database is loaded, first column is eliminated
db <- read.table(file = "DATABASES/place_raw.txt", header = FALSE)
db$V1 <- NULL
#Column names are edited 
colnames(db) <- c("Name", "Motif", "Length", "ID")

#Degenerated nucleotides are designated
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

#Promoter sequence is read according to the filename argument. Typecase is changed to upper
promoter <- read.fasta(file = opt$file, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
#promoter <- read.fasta(file = "atlea45.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
i <- 1
for (i in 1: length(promoter)) promoter[[i]] <- toupper(promoter[[i]])

if ( file.exists( paste0( "./Results_for_",opt$file))){
  unlink(paste0("./Results_for_",opt$file), recursive = TRUE)
}
  
dir.create(path = paste0("./Results_for_",opt$file))
#PREPARACIÓN DEL MOTIVO A BUSCAR
#Se declaran las variables vacias a utilizar
pos_motifs <- list(0)
#Primer ciclo: se va a ciclar a traves de todos los elementos de la columna motif de la db
a <- 1
j <- 1
x <- 1
for (a in 1:length(promoter)) {
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
  	result <- list(words.pos(motif, promoter[[a]]))
  	#El resultado se guarda como un elemento de una lista
  	pos_motifs[j] <- result
  	#Se repite el ciclo y se continuan guardando los resultado encontrados
  }
  #A cada elemento de la lista se le asigna el nombre del elemento cis que le corresponde
  names(pos_motifs) <- db$Name
  #Hay elementos cis que no se encuentran y se graban como character(0). Con esta funcion se eliminan
  cis_founded <- Filter(length, pos_motifs)
  prom_len <- nchar(promoter[[a]])
  result_tb <- lapply(cis_founded, function(n) n-prom_len)
  result_tb <- plyr::ldply(result_tb, rbind)
  result_tb[is.na(result_tb)] <- ""

  setwd(paste0("./Results_for_",opt$file))
  write.table(result_tb, file = paste0(names(promoter[a]),"_results_tb.txt"), quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  write.table(result_tb$.id, file = paste0(names(promoter[a]),"_unique.txt"), quote = FALSE, 
  						row.names = FALSE, col.names = FALSE)
  print(paste0("Printing the table for the promoter of ", names(promoter[a])))
  setwd("..")

print("Printing of DF and unique is done")

db_exp <- readLines(con = "DATABASES/place_database.txt")
res_uni <- scan(paste0("./Results_for_",opt$file,"/",names(promoter[a]),"_unique.txt"), what = "character")
res_uni <- paste0("^", res_uni)

db_exp_filter <- ""
result <- vector(mode = "character")
for (i in 1:length(res_uni)) {
  result <- db_exp[grep(res_uni[i], db_exp) + c(0:3)]
  result <- append(result, "\n")
  db_exp_filter <- append(x = db_exp_filter, values = result)
}

writeLines(db_exp_filter, con = paste0("./Results_for_",opt$file,"/Prom_",names(promoter[a]),"cis_detailed.txt"))
print("Printing of details is done")

}

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