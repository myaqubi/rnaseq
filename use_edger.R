
library(edgeR)
library(biomaRt)
library(readr)
library(tidyverse)
library(magrittr)

## set working directory
### change this based on where your files are
setwd("C:/Users/Moein/Desktop/Jelani")

## prepare data input
### change the ".txt" file name to indicate which raw file you want to work with (e.g. here it is "All ee vs. sh-2.txt")
# counts <- read_delim("All ee vs. sh-2.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  # as.data.frame()
## or if importing from GUI
counts = as.data.frame(Count_table) # change this

View(counts) ## check that data is loaded correctly and 1st column has the IDs
row.names(counts) = counts[[1]] ## makes rownames from 1st column
counts <- counts[, -1] ## removes 1st column because now we have that information as rownames

View(counts) ## check that rownames are correct

## load function ## no change
### if this is not in the working directory specified above, add in full path
### this contains a custom function to use edgeR
source("run_edger_function.R")

## run function
### look at the sourced file above if you want to see what the function does and what changes you can make
dedata <- run_edger(counts_matrix = counts, 
                    genome_name = "hsapiens_gene_ensembl", # change this
                    groups_info = c(rep("control", 3), rep("Treat", 3)),  # change this
                    sample_threshold = 1,
                    file_name = "Elika.txt") # change this (needs file extension)


