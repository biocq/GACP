options(stringsAsFactors = FALSE)
options(warn=-1)

cat("GACP V1.0: identifying genetic alterations that affect cell phenotype from genetic screen data\n\n", file = stdout())
cat("Aim: identifying genetic alterations e.g. mutations that affect cell phenotype, for example, cell proliferation from genetic screen data.\n", file = stdout())
cat("Developed and tested in R 4.1.2.\n", file = stdout())
cat("Copyright: Wenzhou institute, University of Chinese Academy of Sciences. This program is licensed under the MIT License.\n\n", file = stdout())

cat("Install missing packages.\n", file = stdout())
installed_pkgs <- installed.packages()

if(!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	
# The following initializes usage of Bioc devel
if(length(setdiff("fgsea", installed_pkgs)) > 0) {
	BiocManager::install("fgsea")
}
if(length(setdiff("plyr", installed_pkgs)) > 0) {
	install.packages("plyr")
}

if(length(setdiff("dplyr", installed_pkgs)) > 0) {
	install.packages("dplyr")
}

if(length(setdiff("data.table", installed_pkgs)) > 0) {
	install.packages("data.table")
}

if(length(setdiff("boot", installed_pkgs)) > 0) {
	install.packages("boot")
}

if(length(setdiff("snow", installed_pkgs)) > 0) {
	install.packages("snow")
}

if(length(setdiff("parallel", installed_pkgs)) > 0) {
	install.packages("parallel")
}

if(length(setdiff("stringr", installed_pkgs)) > 0) {
	install.packages("stringr")
}

if(length(setdiff("optparse", installed_pkgs)) > 0) {
	install.packages("optparse")
}

if(length(setdiff("collections", installed_pkgs)) > 0) {
	install.packages("collections")
}

if(length(setdiff("RcppAlgos", installed_pkgs)) > 0) {
	install.packages("RcppAlgos")
}

if(length(setdiff("ggpubr", installed_pkgs)) > 0) {
	install.packages("ggpubr")
}

if (length(setdiff("devtools", installed_pkgs)) > 0) {
	install.packages('devtools')
}

if (length(setdiff("cpdetectoR", installed_pkgs)) > 0) {
  devtools::install_github('ontogenerator/cpdetectoR')
}

cat("Load packages.\n", file = stdout())

suppressMessages(library(fgsea))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(snow))
suppressMessages(library(parallel))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(collections))
suppressMessages(library(RcppAlgos))
suppressMessages(library(ggpubr))
suppressMessages(library(cpdetectoR))

option_list = list(
	make_option(c("--folder"), action = "store", default = "ex", type = 'character', help = "Project folder to be used for storing the output for this program."),
	make_option(c("--prefix"), action = "store", default = "ex_", type = 'character', help = "Result file name prefix."),
	make_option(c("--download_depmap_data"), action = "store", default = "FALSE", type = 'logical', help = "If TRUE, download missing DepMap data. [Optional]"),
	make_option(c("--download_PPI_data"), action = "store", default = "FALSE", type = 'logical', help = "If TRUE, download and process missing protein-protein interaction network data for 21Q3 and 19Q3 DepMap data. [Optional]"),
	make_option(c("--ga"), action = "store", default = paste0("./data/mutation_processed.tsv"), type='character', help = "Tabular file of genetic alteration data (i.e. mutations)."),
	make_option(c("--exp"), action = "store", default = "./data/expr_processed.tsv", type = 'character', help= "Tabular file of gene expression data."),
	make_option(c("--dep"), action = "store", default = "./data/dependency_processed.tsv", type = 'character', help = "Tabular file of gene dependency data. [Optional]"),
	make_option(c("--ppi"), action = "store", default = "./data/STRING_PPI_21Q3_processed.txt", type = 'character', help = "Tabular file of two columns of genes, i.e. two-column tab delimited undirected PPI interactions. [Optional]"),
	make_option(c("--start_PPI"), action = "store", default = 1001, type = 'integer', help = "Starting line of --ppi file. Set to -1 if no subsetting is needed [Optional]"),
	make_option(c("--end_PPI"), action = "store", default = 2000, type = 'integer', help = "Ending line of --ppi file. Set to -1 if no subsetting is needed [Optional]"),
	make_option(c("--core"), action = "store", default = 1, type = 'integer', help = "CPU cores used to run this script concurrently. [Optional]"),
	make_option(c("--seed"), action = "store", default = 1, type = 'integer', help = "Seed used in set.seed(). Set to -1 to use rownumbers of PPI file as seeds. [Optional]"),
	make_option(c("--npermut"), action = "store", default = 1000, type = 'integer', help = "Number of permutations for fgsea. [Optional]"),
	make_option(c("--nproc"), action = "store", default = 0, type = 'integer', help = "Number of processes to be used for fgsea. [Optional]"),
	make_option(c("--min_size"), action = "store", default = 5, type = 'integer', help = "Min cell line number used in fgsea. [Optional]"),
	make_option(c("--max_size"), action = "store", default = 200, type = 'integer', help = "Max cell line number used in fgsea. [Optional]"),
	make_option(c("--include_composite_mutation_analysis"), action = "store", default = "TRUE", type = 'logical', help = "If TRUE, all check the composite mutation regulation (mutation number >= 1). [Optional]"),
	make_option(c("--max_combination"), action = "store", default = 4, type = 'integer', help = "Max combination number of genetic alteration combinations. [Optional]"),
	make_option(c("--include_cis"), action = "store", default = TRUE, type = 'logical', help = " Cis- genetic alteration results are also reported. [Optional]"),

	make_option(c("--use_all_cell_lines"), action = "store", default = TRUE, type = 'logical', help = "Whether use all cell-lines. TRUE: all; FALSE: cell-lines with mutations identified in Step 1. [Optional]"),
	make_option(c("--permutation"), action = "store", default = TRUE, type = 'logical', help = "Whether perform permutation to obtain p-values in step 2. TRUE is required for the following Steps. [Optional]"),
	make_option(c("--nperm"), action = "store", default = 100, type = 'integer', help = "Number of permutations. [Optional]"),
	make_option(c("--precise_mediation_pvalue"), action = "store", default = FALSE, type = 'logical', help = "Whether more precise p-values are needed (TRUE: up to 100 x nperm permutations; FALSE: nperm permutations). [Optional]"),


	make_option(c("--remove_cis_confounding"), action = "store", default = FALSE, type = 'logical', help = "Whether to remove trans records with potential cis confounding regulation in CSEA analysis. [Optional]"),
	make_option(c("--cis_confounding_pvalue_cutoff"), action = "store", default = 0.1, type = 'double', help = "Remove trans records that have adjusted p-values above the indicated cutoff in CSEA analysis. [Optional]"),

	make_option(c("--remove_trans_results_from_cis_sig"), action = "store", default = FALSE, type = 'logical', help = "Whether skip Gene-A trans-CSEA calculation if cis-CSEA of Gene A is significant (q-value < 0.1). [Optional]"),
	make_option(c("--include_leading_mutation_only"), action = "store", default = TRUE, type = 'logical', help = "Whether to include the most significant CSEA p-value result for each PPI interaction record. [Optional]"),
	make_option(c("--full_result"), action = "store", default = TRUE, type = 'logical', help = "Full results will be reported. [Optional]"),

	make_option(c("--cna"), action = "store", default = "./data/CNV_processed.tsv", type = 'character', help= "Tabular file of CNA data. [Optional]"),
	make_option(c("--cna_adjust"), action = "store", default = FALSE, type = 'logical', help = "Whether remove CNA confounding effect. Specify --cna also when using this argument. [Optional]"),
	make_option(c("--cna_adjust_Pvalue_cutoff"), action = "store", default = 0.01, type = 'double', help = "Remove CNA confounding effect according to the P-value cutoff. Specify --cna and --cna_adjust also when using this argument. [Optional]"),
	make_option(c("--cna_permutation_number"), action = "store", default = 50, type = 'integer', help = "Calculate empirical p-values based on the predefined permutation size. Specify --cna and --cna_adjust also when using this argument. [Optional]"),
	make_option(c("--cna_discretization"), action= "store", default = FALSE, type = 'logical', help = "Whether discretize the CNA values. Specify --cna and --cna_adjust also when using this argument. [Optional]"),
	make_option(c("--cna_amplification_cutoff"), action = "store", default = 1, type = 'double', help = "CN > this cutoff is considered as amplication. Specify --cna, --cna_adjust and --cna_discretization as well when using this argument. [Optional]"),
	make_option(c("--cna_deletion_cutoff"), action = "store", default = 0.5, type = 'double', help = "CN < this cutoff is considered as deletion. Specify --cna, --cna_adjust and --cna_discretization as well when using this argument. [Optional]"),
	make_option(c("--FDR_cutoff_cis"), action = "store", default = 0.1, type = 'double', help = "The score cutoff for defining cis genetic alterations. [Optional]"),
	make_option(c("--FDR_cutoff_trans"), action = "store", default = 0.1, type = 'double', help = "The score cutoff for defining trans genetic alterations. [Optional]"),
	make_option(c("--full_report"), action = "store", default = TRUE, type = 'logical', help = "Output all mutations and associated information."),
	make_option(c("--BED_report"), action = "store", default = TRUE, type = 'logical', help = "Output BED format file name that can be used in CRAVAT server.")

)


arguments = parse_args(OptionParser(option_list = option_list))

if(arguments$start_PPI > 0 & arguments$end_PPI > 0){
	arguments$prefix <- paste0(arguments$prefix, arguments$start_PPI, "_", arguments$end_PPI, "_")
}


# Print Parameters
cat("PARAMETERS:\n", file = stdout())
for(arg.idx in 1 : length(arguments)) {
	cat(paste0(names(arguments)[arg.idx], "\t", arguments[[arg.idx]], "\n"), file = stdout())
}
cat("------------------\n", file = stdout())

output_dir <- file.path(getwd(), arguments$folder)
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}
cat("Using ", output_dir," as the working directory\n", file = stdout())

######### Step 1: CSEA analysis #########
cat("Step 1: reading data.\n", file = stdout())


if(arguments$download_depmap_data){

	if(!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	if(length(setdiff("biomaRt", installed_pkgs)) > 0) {
		BiocManager::install("biomaRt")
	}
	if(length(setdiff("GenomicFeatures", installed_pkgs)) > 0) {
		BiocManager::install("GenomicFeatures")
	}
	if(length(setdiff("depmap", installed_pkgs)) > 0) {
		BiocManager::install("depmap", version = '3.14')
	}
	if(length(setdiff("RaggedExperiment", installed_pkgs)) > 0) {
		BiocManager::install("RaggedExperiment")
	}

	if(length(setdiff("data.table", installed_pkgs)) > 0) {
		install.packages(pkgs = setdiff("data.table", installed_pkgs))
	}
	if(length(setdiff("igraph", installed_pkgs)) > 0) {
		install.packages(pkgs = setdiff("igraph", installed_pkgs))
	}

	# Load library.
	cat("Load packages.\n", file = stdout())
	suppressMessages(library(depmap))
	suppressMessages(library(ExperimentHub))
	suppressMessages(library(data.table))
	suppressMessages(library(igraph))
	suppressMessages(library(RaggedExperiment))
	suppressMessages(library(stringr))


	eh <- ExperimentHub()


	# Generate 19Q3 CNA, RNAi expression and mutation data


	if(!file.exists("./data/CNV_19Q3_processed.tsv")){
		cat("Download and process DepMap 19Q3 CCLE CNA data.\n", file = stdout())
	  CNA_19Q3_dat <- eh[["EH3082"]] # CNA_19Q3
	  write.table(CNA_19Q3_dat[, c(5, 3, 6)], file = "./data/CNV_19Q3_processed.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}

	if(!file.exists("./data/RNAi_processed.tsv")){
		cat("Download and process DepMap 19Q3 RNAi screen data.\n", file = stdout())
	  RNAi_dat <- eh[["EH3080"]] # RNAi_19Q3
	  RNAi_dat <- RNAi_dat[!is.na(RNAi_dat$dependency), c("gene_name", "dependency", "cell_line")]
	  colnames(RNAi_dat) <- c("gene", "dependency", "cell_line")
	  write.table(RNAi_dat,  file = "./data/RNAi_processed.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}

	if(!file.exists("./data/expr_19Q3_processed.tsv")){
		cat("Download and process DepMap 19Q3 expression data.\n", file = stdout())
		expr_19Q3_dat <- eh[["EH3084"]] # expr_19Q3

		expr_19Q3_dat <- expr_19Q3_dat[, c(5, 3, 6)]
		colnames(expr_19Q3_dat) <- c("gene", "rna_expression", "cell_line")

		expr_19Q3_dat <- expr_19Q3_dat %>% filter(!grepl('MERGED', cell_line))
		expr_19Q3_dat$context <- sapply(expr_19Q3_dat$cell_line, function(z) paste(strsplit(z, split = "_")[[1]][-1], collapse = "_"))
		if(!file.exists("./data/expr_19Q3_processed.tsv")){
			write.table(expr_19Q3_dat, file = "./data/expr_19Q3_processed.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		}
	     
	}

	if(!file.exists("./data/mutation_matching_RNAi_processed.tsv")){
		cat("Download and process DepMap 19Q3 mutation data.\n", file = stdout())
	  mutation_dat <- eh[["EH3085"]] # mutationCalls_19Q3
	  cat("Download 19Q3 meta data.\n", file = stdout())
	  metadata <- eh[["EH3086"]] # metadata of 19Q3
	  cell_line_mapping <- metadata[, c("depmap_id", "cell_line", "primary_or_metastasis", "primary_disease")]
	  mutation_dat <- merge(mutation_dat, cell_line_mapping, by = 'depmap_id', all.x = TRUE)
	  mutation_dat <- mutation_dat %>% filter(!grepl('MERGED', cell_line))
	 
	  mutation_dat$var_class <- mutation_dat$var_class %>% str_replace_all("Missense_Mutation", "Missense")
	  mutation_dat$var_class <- mutation_dat$var_class %>% str_replace_all("Nonsense_Mutation", "Nonsense")

	  mutation_dat2 <- mutation_dat[, c("gene_name", "var_class", "cell_line", "genome_change", "protein_change", "annotation_transcript")]
	  mutation_dat2$var_class[which(mutation_dat2$var_class %in% names(table(mutation_dat2$var_class)[table(mutation_dat2$var_class) < 3000]))] <- ""
	  mutation_dat2 <- mutation_dat2[mutation_dat2$var_class != "" & !is.na(mutation_dat2$var_class),]
	  count <- as.data.frame(table(mutation_dat2$genome_change))
		kept_mut <- as.character(count$Var1[which(count$Freq > 2)])
		mutation_dat2 <- mutation_dat2[mutation_dat2$genome_change %in% kept_mut,]
	  write.table(mutation_dat2, file = "./data/mutation_matching_RNAi_processed.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}

	# Generate 21Q3 CNA, RNAi screen, gene expression and mutation data

	if(file.exists("./data/metadata.rds")){
	  metadata <- readRDS("./data/metadata.rds");
	  
	}else{
		cat("Download meta data.\n", file = stdout())
	  metadata <- eh[["EH6757"]]
	  saveRDS(metadata, file = './data/metadata.rds')
	}

	if(!file.exists("./data/CNA_processed.tsv")){
		cat("Download and process DepMap 21Q3 CCLE CNA data.\n", file = stdout())
		CNA_dat <- eh[["EH6754"]]
		CNA_dat <- CNA_dat[, c(5, 3, 6)]
		CNA_dat <- CNA_dat %>% filter(!grepl('MERGED', cell_line))
		write.table(CNA_dat, file = paste0("./data/CNA_processed.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}


	if(!file.exists("./data/dependency_processed.tsv")){
		cat("Download and process DepMap 21Q3 CRISPR screen data.\n", file = stdout())
		dependency_dat <- eh[["EH6753"]]
		dependency_dat <- dependency_dat[, c(5, 3, 6)]
		colnames(dependency_dat) <- c("gene", "dependency", "cell_line")
		dependency_dat <- dependency_dat %>% filter(!grepl('MERGED', cell_line))
		dependency_dat$context <- sapply(dependency_dat$cell_line, function(z) paste(strsplit(z, split = "_")[[1]][-1], collapse = "_"))
		write.table(dependency_dat, file = paste0("./data/dependency_processed.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	}

	if(!file.exists("./data/expr_processed.tsv")){
		cat("Download and process DepMap 21Q3 expression data.\n", file = stdout())
		exp_tpm_dat <- eh[["EH6755"]]

		exp_tpm_dat <- exp_tpm_dat[, c(5, 3, 6)]
		colnames(exp_tpm_dat) <- c("gene", "rna_expression", "cell_line")

		exp_tpm_dat <- exp_tpm_dat %>% filter(!grepl('MERGED', cell_line))
		exp_tpm_dat$context <- sapply(exp_tpm_dat$cell_line, function(z) paste(strsplit(z, split = "_")[[1]][-1], collapse = "_"))
		if(!file.exists(paste("./data/expr_processed.tsv", sep = ""))){
		    write.table(exp_tpm_dat, file = paste0("./data/expr_processed.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		}
	}



	if(!file.exists("./data/mutation_processed.tsv")){
		cat("Download and process DepMap 21Q3 mutation data.\n", file = stdout())
	  mutation_dat <- eh[["EH6756"]] # mutationCalls_21Q3
	 
	  cell_line_mapping <- metadata[, c("depmap_id", "cell_line", "primary_or_metastasis", "primary_disease")]
	  mutation_dat <- merge(mutation_dat, cell_line_mapping, by = 'depmap_id', all.x = TRUE)
	  mutation_dat <- mutation_dat %>% filter(!grepl('MERGED', cell_line))
	 
	  mutation_dat$var_class<-mutation_dat$var_class %>% str_replace_all("Missense_Mutation", "Missense")
	  mutation_dat$var_class<-mutation_dat$var_class %>% str_replace_all("Nonsense_Mutation", "Nonsense")

	  mutation_dat2 <- mutation_dat[, c("gene_name", "var_class", "cell_line", "genome_change", "protein_change", "annotation_trans")]
	  mutation_dat2$var_class[which(mutation_dat2$var_class %in% names(table(mutation_dat2$var_class)[table(mutation_dat2$var_class) < 3000]))] <- ""
	  mutation_dat2 <- mutation_dat2[mutation_dat2$var_class != "" & !is.na(mutation_dat2$var_class),]
	  count <- as.data.frame(table(mutation_dat2$genome_change))
	  kept_mut <- as.character(count$Var1[which(count$Freq > 2)])
	  mutation_dat2 <- mutation_dat2[mutation_dat2$genome_change %in% kept_mut,]
	  write.table(mutation_dat2, file = paste0("./data/mutation_processed.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
}

if(arguments$download_PPI_data){
	# Generate 21Q3 PPI file
	if(file.exists("./data/STRING_PPI_21Q3_processed.txt")){
		stop(paste0("File ", arguments$ppi, " already exists, no need to download PPI data, set --download_PPI_data=FALSE."))
	}
	if(!dir.exists("./data")){
		dir.create("./data")
	}
	if(length(setdiff("STRINGdb", installed_pkgs)) > 0) {
		BiocManager::install("STRINGdb")
		string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold = 0, input_directory = "./data")
		string_db
	}
	
	if(file.exists("./data/mutation_processed.tsv")){
		mutations_dt <- fread("./data/mutation_processed.tsv", header = TRUE, sep = "\t")
		mut_genes <- unique(mutations_dt$gene_name)
		rm(mutations_dt)
	}else{
		stop(paste0("File path ", "./data/mutation_processed.tsv", " is not correctly specified."))
	}

	if(file.exists("./data/expr_processed.tsv")){
		exp_tpm_dt <- fread("./data/expr_processed.tsv", header = TRUE, sep = "\t")
		expr_genes <- unique(exp_tpm_dt$gene)
		rm(exp_tpm_dt)
	}else{
		stop(paste0("File path ", "./data/expr_processed.tsv", " is not correctly specified."))
	}
	
	query <- data.frame("gene" = expr_genes)
	query_mapped <- string_db$map(query, "gene", removeUnmappedRows = TRUE)
	stringdb <- fread("./data/9606.protein.links.v11.5.txt", header = TRUE, sep = '\t')

	query_interaction <- merge(stringdb, query_mapped, by.x = "protein1", by.y = "STRING_id")

	query_interaction <- merge(query_interaction, query_mapped, by.x = "protein2", by.y = "STRING_id")
	query_interaction <- query_interaction[, c("gene.x", "gene.y")]

	colnames(query_interaction) <- c("V1", "V2")
	query_interaction <- query_interaction[order(query_interaction$V1),]
	query_interaction2 <- query_interaction[which(query_interaction$V1 %in% mut_genes | query_interaction$V2 %in% mut_genes),] # either column should have mutated genes

	query_interaction2 <- query_interaction2[intersect(which(query_interaction2$V1 %in% expr_genes) , which(query_interaction2$V2 %in% expr_genes)),]


	write.table(query_interaction2, file = paste0("./data/STRING_PPI_with_mut_genes.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	rm(query_interaction)
	rm(query_interaction2)
	rm(query_mapped)
	
	# Generate 19Q3 PPI file
	if(file.exists("./data/mutation_matching_RNAi_processed.tsv")){
		mutations_dt <- fread("./data/mutation_matching_RNAi_processed.tsv", header = TRUE, sep = "\t")
		mut_genes <- unique(mutations_dt$gene_name)
		rm(mutations_dt)
	}else{
		stop(paste0("File path ", "./data/mutation_matching_RNAi_processed.tsv", " is not correctly specified."))
	}

	if(file.exists("./data/expr_19Q3_processed.tsv")){
		exp_tpm_dt <- fread("./data/expr_19Q3_processed.tsv", header = TRUE, sep = "\t")
		expr_genes <- unique(exp_tpm_dt$gene)
		rm(exp_tpm_dt)
	}else{
		stop(paste0("File path ", "./data/expr_19Q3_processed.tsv", " is not correctly specified."))
	}
	
	query <- data.frame("gene" = expr_genes)
	query_mapped <- string_db$map(query, "gene", removeUnmappedRows = TRUE)
	stringdb <- fread("./data/9606.protein.links.v11.5.txt", header = TRUE, sep = '\t')

	query_interaction <- merge(stringdb, query_mapped, by.x = "protein1", by.y = "STRING_id")

	query_interaction <- merge(query_interaction, query_mapped, by.x = "protein2", by.y = "STRING_id")
	query_interaction <- query_interaction[, c("gene.x", "gene.y")]

	colnames(query_interaction) <- c("V1", "V2")
	query_interaction <- query_interaction[order(query_interaction$V1),]
	query_interaction2 <- query_interaction[which(query_interaction$V1 %in% mut_genes | query_interaction$V2 %in% mut_genes),] # either column should have mutated genes

	query_interaction2 <- query_interaction2[intersect(which(query_interaction2$V1 %in% expr_genes) , which(query_interaction2$V2 %in% expr_genes)),]


	write.table(query_interaction2, file = paste0("./data/STRING_PPI_with_mut_genes_RNAi.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	rm(query_interaction)
	rm(query_interaction2)
	rm(query_mapped)
	
	
	######
	PPI <- fread(paste0("./data/STRING_PPI_with_mut_genes.txt"), header = TRUE, sep = '\t')
	gene_dict <- dict()
	PPI_nonredundant_tmp <- NULL

	index <- 1
	pairs <- paste(PPI[index, "V1"], PPI[index, "V2"], collapse = ",")
	pairs_rev <- paste(PPI[index, "V2"], PPI[index, "V1"], collapse = ",")
	write.table(PPI[index,], file = "./data/STRING_PPI_with_mut_genes_nonredundant.txt", sep = "\t", col.names = TRUE, row.names = FALSE,quote=FALSE)
	  
	gene_dict$set(pairs, "");
	gene_dict$set(pairs_rev, "");

	for(index in 2:nrow(PPI)){
	  pairs <- paste(PPI[index, "V1"], PPI[index, "V2"], collapse = ",")
	  pairs_rev <- paste(PPI[index, "V2"], PPI[index, "V1"], collapse = ",")
	  if(gene_dict$has(pairs)|gene_dict$has(pairs_rev)){
	  }else{
	    PPI_nonredundant_tmp <- rbind(PPI_nonredundant_tmp, PPI[index,])
	    gene_dict$set(pairs, "");
	    gene_dict$set(pairs_rev, "");
	  }
	  if((index %% 4000) == 0){
	    write.table(PPI_nonredundant_tmp, file = "./data/STRING_PPI_with_mut_genes_nonredundant.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
	    
	    PPI_nonredundant_tmp <- NULL
	    cat(paste0("#", index, "/", nrow(PPI), "\n"), file = stdout())
	  }
	}
	write.table(PPI_nonredundant_tmp, file = "./data/STRING_PPI_with_mut_genes_nonredundant.txt", sep = "\t",
	            col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
	            
	#####
	
	PPI <- fread("./data/STRING_PPI_with_mut_genes_RNAi.txt", header = TRUE, sep = '\t')
	gene_dict <- dict()
	PPI_nonredundant_tmp <- NULL

	index <- 1
	pairs <- paste(PPI[index, "V1"], PPI[index, "V2"], collapse = ",")
	pairs_rev <- paste(PPI[index, "V2"], PPI[index, "V1"], collapse = ",")
	write.table(PPI[index,], file=paste0(Prefix, "./data/STRING_PPI_with_mut_genes_RNAi_nonredundant.txt"),sep="\t", col.names = TRUE, row.names = FALSE,quote=FALSE)
	  
	gene_dict$set(pairs, "");
	gene_dict$set(pairs_rev, "");

	for(index in 2:nrow(PPI)){
	  pairs <- paste(PPI[index, "V1"], PPI[index, "V2"], collapse = ",")
	  pairs_rev <- paste(PPI[index, "V2"], PPI[index, "V1"], collapse = ",")
	  if(gene_dict$has(pairs)|gene_dict$has(pairs_rev)){
	  }else{
	    PPI_nonredundant_tmp <- rbind(PPI_nonredundant_tmp, PPI[index,])
	    gene_dict$set(pairs, "");
	    gene_dict$set(pairs_rev, "");
	  }
	  if((index %% 4000) == 0){
	    write.table(PPI_nonredundant_tmp, file = "./data/STRING_PPI_with_mut_genes_RNAi_nonredundant.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
	    
	    PPI_nonredundant_tmp <- NULL
	    cat(paste0("#", index, "/", nrow(PPI), "\n"), file = stdout())
	  }
	}
	write.table(PPI_nonredundant_tmp, file = "./data/STRING_PPI_with_mut_genes_RNAi_nonredundant.txt", sep = "\t",
	            col.names = FALSE,row.names = FALSE,quote=FALSE, append = TRUE)
	
	
	
	####
	
	if(!dir.exists("./corr")){
	    dir.create("./corr")
	}
	
	corr_gene_expr <- function(index = NA, Gene_A = "", Gene_B = ""){

    tmp_data_frame <- data.frame("index" = vector("numeric", length = 0), "Gene_A" = vector("character", length = 0), "Gene_B" = vector("character", length = 0), "corr" = vector("numeric",length = 0), "P_val" = vector("numeric",length = 0))
    
    if(is.null(key(exp_tpm_dt))){
      setkey(exp_tpm_dt, gene)
    }
    exp_dat_gene_1 <- exp_tpm_dt[.(Gene_A),]
    exp_dat_gene_2 <- exp_tpm_dt[.(Gene_B),]
    exp_dat_gene_12 <-merge(exp_dat_gene_1, exp_dat_gene_2, by="cell_line")
    exp_dat_gene_12 <- exp_dat_gene_12[!is.na(exp_dat_gene_12$cell_line),]
    if(nrow(exp_dat_gene_2) < 1){
      return(tmp_data_frame)
    }
    
    result <- cor.test(exp_dat_gene_12$rna_expression.x, exp_dat_gene_12$rna_expression.y, method = "spearman",exact=FALSE)
    
    tmp_data_frame <- data.frame("index" = index, "Gene_A" = Gene_A, "Gene_B" = Gene_B, "corr" = result$estimate, "P_val" = result$p.value)
    return(tmp_data_frame)
	}
		
	PPI_nonredundant_filename <- c("./data/STRING_PPI_with_mut_genes_nonredundant.txt", "./data/STRING_PPI_with_mut_genes_RNAi_nonredundant.txt")
	PPI_nonredundant_prefix <- c("21Q3_", "19Q3_")
	for(j in 1:2){
		PPI_nonredundant <- fread(PPI_nonredundant_filename[j], header = TRUE, sep ='\t')
		PPI_nonredundant$index <- 1:nrow(PPI_nonredundant)
		PPI_nonredundant <- PPI_nonredundant[,c("index","V1","V2")]
		if(arguments$start_PPI > 0 & arguments$end_PPI > 0 & arguments$end_PPI > arguments$start_PPI){
		    PPI_nonredundant <- PPI_nonredundant[arguments$start_PPI : arguments$end_PPI,]
		}
		if(arguments$start_PPI > 0 & arguments$end_PPI > 0){
		  arguments$prefix <- paste0(arguments$prefix, arguments$start_PPI, "_", arguments$end_PPI, "_")
		}

		exp_tpm_dt <- fread(arguments$exp, header =TRUE, sep ='\t')
		exp_tpm_dt <- exp_tpm_dt[, c("gene", "rna_expression", "cell_line")]
		setkey(exp_tpm_dt, gene)

		corr_gene_expr_result <- NULL


		if(arguments$core < 2){

				for(i in 1:nrow(PPI_nonredundant)){
				  corr_gene_expr_result <- rbind(corr_gene_expr_result, corr_gene_expr(PPI_nonredundant[i, 1], PPI_nonredundant[i,2], PPI_nonredundant[i,3]))
				}
				write.table(corr_gene_expr_result, file=paste0("./corr", "/STRING_PPI_", PPI_nonredundant_prefix[j], "with_mut_genes_nonredundant_with_corr.txt"),sep="\t",col.names = F, row.names = FALSE, quote=FALSE)

		}else{

				options(cl.cores = arguments$core)
				this.cluster <- makeCluster(getOption("cl.cores", 2))
				clusterCall(cl = this.cluster, fun = function(){
				  library(dplyr);
				  library(data.table)
				})

				if(formalArgs(clusterExport)[2] %in% "list"){
				    clusterExport(cl = this.cluster, list = c("PPI_nonredundant", "corr_gene_expr_result", "corr_gene_expr", "exp_tpm_dt"))
				}else{
				    clusterExport(cl = this.cluster, varlist = c("PPI_nonredundant", "corr_gene_expr_result", "corr_gene_expr", "exp_tpm_dt"), envir = environment())
				}

				corr_gene_expr_result <- 
				  parLapply(cl = this.cluster,
				  1:nrow(PPI_nonredundant), # nrow(PPI_nonredundant)
				  function(idx) {
				    rtn <- corr_gene_expr(PPI_nonredundant[idx, 1], PPI_nonredundant[idx, 2], PPI_nonredundant[idx, 3])
				    return(rtn)
				  }
				)

				stopCluster(this.cluster)

				for( j in 1:length(corr_gene_expr_result)){
				  final_corr_gene_expr_result <- rbind(final_corr_gene_expr_result, corr_gene_expr_result[[j]])
				}

				final_corr_gene_expr_result$P_adj <- p.adjust(final_corr_gene_expr_result$P_val, method = "BH")

				write.table(final_corr_gene_expr_result, file = paste0("./corr", "/", PPI_nonredundant_prefix, PPI_nonredundant_prefix[j], "with_mut_genes_nonredundant_with_corr.txt"),sep="\t",col.names = F, row.names = FALSE, quote=FALSE)
				# sbatch -N 1 -n 10 --mem=10000MB --job-name="corr" -t 4-15:00:00 -p wholenodeQ --wrap="Rscript MutCP_exp_correlation.R --core=10"

		}
		
		previouswd <- getwd()
		setwd(paste0(previouswd, "/corr/"))
		if(length(list.files())){
			final_results <- rbindlist(lapply(list.files(), fread))
			colnames(final_results) <- c("index", "V1", "V2", "corr", "P_val")
			final_results$index <- as.numeric(final_results$index)
			final_results$P_val <- as.numeric(final_results$P_val)
			if(is.unsorted(final_results$index)){
	  		final_results <- final_results[order(final_results$index),]
	  	}
	  	
	  	final_results$P_adj <- p.adjust(final_results$P_val, method = "BH")
	  	
	  	# write.table(final_results, file = paste0("../data/", "STRING_PPI_", PPI_nonredundant_prefix, "with_mut_genes_nonredundant_with_corr.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	  	write.table(final_results[final_results$P_adj < 0.1, c("V1", "V2")], file = paste0("../data/", "STRING_PPI_", PPI_nonredundant_prefix, "processed.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		}
		setwd(previouswd)
	}
	cat("PPI files are generated completedly.\n", file = stdout())
	
}
######### Step 1.1: Load omics data #########

Input.net <- NULL
mutations_dt <- NULL
exp_tpm_dt <- NULL
PPI_genes <- NULL
if(file.exists(arguments$ppi)){
	Input.net <- read.table(arguments$ppi, header = TRUE, sep = "\t")
	if(arguments$start_PPI > 0 & arguments$end_PPI > 0 & arguments$end_PPI > arguments$start_PPI){
		Input.net <- Input.net[arguments$start_PPI : arguments$end_PPI,]
	}else{
		if(arguments$end_PPI < arguments$start_PPI){
			stop("End_PPI should be no less than start_PPI and both should be larger than zero.")
		}
	}
	PPI_genes <- unique(c(Input.net$V1, Input.net$V2))
 
}else{
	stop(paste0("File path ", arguments$ppi, " is not correctly specified."))
}

cat("Read", paste0(arguments$ppi, " with ", nrow(Input.net), " rows after filtering.\n"), file = stdout())

if(file.exists(arguments$ga)){
	mutations_dt <- fread(arguments$ga, header = TRUE, sep = "\t")
	if(!("cell_line" %in% colnames(mutations_dt)) & !("-" %in% colnames(mutations_dt))){
		stop(paste0("Data ", arguments$ga, " does not contain columns of cell_line or disease."))
	}
	if(length(colnames(mutations_dt)) < 3){
		stop(paste0("Data ", arguments$ga, " does not contain necessary columns."))
	}
	setkey(mutations_dt, gene_name)
	if(nrow(Input.net) > 0){
		mutations_dt <- mutations_dt[.(PPI_genes)]
	}
	setkey(mutations_dt, gene_name)
}else{
	stop(paste0("File path ", arguments$ga, " is not correctly specified."))
}

cat("Read", paste0(arguments$ga, " with ", nrow(mutations_dt), " rows after filtering\n"), file = stdout())


if(file.exists(arguments$exp)){
	exp_tpm_dt <- fread(arguments$exp, header = TRUE, sep = "\t")
	if(!("cell_line" %in% colnames(exp_tpm_dt)) & !("-" %in% colnames(exp_tpm_dt))){
		stop(paste0("Data ", arguments$exp, " does not contain columns of cell_line or disease."))
	}
	if(grepl("dependency", arguments$exp) | grepl("RNAi", arguments$exp)){
		exp_tpm_dt <- exp_tpm_dt[, c("gene", "dependency", "cell_line")]
		colnames(exp_tpm_dt) <- c("gene", "rna_expression", "cell_line")
	}
	exp_tpm_dt <- exp_tpm_dt[, c("gene", "rna_expression", "cell_line")]
	setkey(exp_tpm_dt, gene)
	if(nrow(Input.net) > 0){
		exp_tpm_dt <- exp_tpm_dt[.(PPI_genes)]
	}
	setkey(exp_tpm_dt, gene)
}else{
	stop(paste0("File path ", arguments$exp, " is not correctly specified."))
}
cat("Read", paste0(arguments$exp, " with ", nrow(exp_tpm_dt), " rows after filtering.\n"), file = stdout())

######### Step 1.2: CSEA analysis #########


getCombinationAtSpecifiedRepetition <- function(z, n) {
	result <- NULL
	tmp <- apply(RcppAlgos::comboGeneral(z, n), 1, paste, collapse = ",")
	result <- c(result, tmp)
	unique(result)
}

CSEA <- function(index = NA, Gene_A = "", Gene_B = ""){
	
	# Gene_A: gene with mutations
	# Gene_B: gene with potential target
	tmp_data_frame <- data.frame("index" = vector("numeric", length = 0), "Gene_A" = vector("character", length = 0), "Gene_B" = vector("character", length = 0), "Mutation_type" = vector("character", length = 0), "Mutation_genome_change" = vector("character", length = 0), "Mutation_protein_change" = vector("character", length = 0), "ES" = vector("numeric", length = 0), "NES" = vector("numeric", length = 0), "P_val" = vector("numeric", length = 0), "cell_line" = vector("character", length = 0))

	if(is.null(key(exp_tpm_dt))){
		setkey(exp_tpm_dt, gene)
	}
	if(is.null(key(mutations_dt))){
		setkey(mutations_dt, gene_name)
	}
	mutation_df <- mutations_dt[.(Gene_A), c("var_class", "cell_line", "genome_change", "protein_change")]
	mutation_unique_Gene_A <- unique(mutation_df$genome_change)

	exp_dat_gene_2 <- exp_tpm_dt[.(Gene_B),]
	exp_dat_gene_2 <- exp_dat_gene_2[!duplicated(exp_dat_gene_2),]
	if(nrow(exp_dat_gene_2) < 1){
		return(tmp_data_frame)
	}

	exp_dat_gene_2 <- exp_dat_gene_2[!is.na(exp_dat_gene_2$cell_line),]
	RNK2 <- exp_dat_gene_2$rna_expression
	names(RNK2) <- exp_dat_gene_2$cell_line
	RNK2 <- sort(RNK2, decreasing = FALSE) #Ranking from low to high; ES > 0, enriched in higher expression


	# for expressed cell lines
	mutations_genes <- mutation_df[which(mutation_df$cell_line %in% names(RNK2)),]
	ctx_mut_genes <- list()

	for(i in 1 : length(mutation_unique_Gene_A)){
  
	  tmp <- mutations_genes[mutations_genes$genome_change %in% mutation_unique_Gene_A[i],]
	  if(length(tmp$cell_line) > 4){
		 	ctx_mut_genes[[mutation_unique_Gene_A[i]]] <- tmp$cell_line
	  }
  
	}
	if(length(ctx_mut_genes) < 1 & length(mutation_unique_Gene_A) > 1 & arguments$include_composite_mutation_analysis){ # 2 combination
	 
		for(i in 2:length(mutation_unique_Gene_A)){
			if(i > arguments$max_combination){
				break;
			}
			comb <- getCombinationAtSpecifiedRepetition(mutation_unique_Gene_A, i)

			for(i in 1 : length(comb)){
			  if(grepl(",", comb[i])){
					split <- strsplit(comb[i], split = ",")[[1]];
					combined <- c()
					for(j in 1 : length(split)){
						combined <- c(combined, split[j]);
					}
					tmp <- mutations_genes[mutations_genes$genome_change %in% combined,]
					if(length(tmp$cell_line) > 4){
						ctx_mut_genes[[comb[i]]] <- tmp$cell_line
					}
			    
			  }
		  
			}
			
			if(length(ctx_mut_genes) > 0){
				break;
			}
		}
	 
	}


	if(length(RNK2) > 0){
		set.seed(arguments$seed)
		res_std <- fgseaMultilevel(pathways = ctx_mut_genes, stats = RNK2, minSize = arguments$min_size, maxSize = arguments$max_size, eps = 1e-8, sampleSize = 101, nPermSimple = arguments$npermut, nproc = arguments$nproc, scoreType = "std")


		if(nrow(res_std) > 0){
		  for(j in 1 : nrow(res_std)){
				if(length(ctx_mut_genes[j][[1]]) > 4){
				  if(grepl(",", res_std$pathway[j])){
			    	split <- strsplit(res_std$pathway[j], split = ",")[[1]];
			    	info <- mutation_df[which(mutation_df$genome_change %in% split),]
				  }else{
					 	info <- mutation_df[which(mutation_df$genome_change == res_std$pathway[j]),]
				  }
				  
				  mut_type <- unique(info$var_class)
				  mut_protein <- unique(info$protein_change)
				  if(length(mut_type) > 1){
						mut_type <- paste(mut_type, collapse = ",")
				  }
				  if(length(mut_protein) > 1){
						mut_protein <- paste(mut_protein, collapse = ",")
				  }
				  
			    tmp_data_frame <- rbind(tmp_data_frame, data.frame("index" = index, "Gene_A" = Gene_A, "Gene_B" = Gene_B, "Mutation_type" = mut_type, "Mutation_genome_change" = res_std$pathway[j], "Mutation_protein_change" = mut_protein, "ES" = round(res_std$ES[j], 3), "NES" = round(res_std$NES[j], 3), "P_val" = res_std$pval[j], "cell_line" = paste(ctx_mut_genes[j][[1]], collapse = ',')))

			  }
			}
		}
	}
	if(arguments$include_leading_mutation_only){
	  tmp_data_frame <- tmp_data_frame[which(tmp_data_frame$P_val == min(tmp_data_frame$P_val, na.rm = T)),]
	}

	return(tmp_data_frame)
}



params <- NULL
gene_dict <- dict()

row_ppi <- nrow(Input.net)


params_tmp <- NULL
omit_genes_for_trans_detection <- NA

CSEA_results <- NULL

if(arguments$include_cis){
	cat("Construct a parameter table for cis-genetic alteration CSEA analysis.\n", file = stdout())
	for(index in 1:length(PPI_genes)){
	  
  	if(length(mutations_dt$gene_name %in% PPI_genes[index]) > 0){ # try gene A as mutation gene
	 	 	params <- rbind(params,data.frame("index" = 0, "Gene_A" = PPI_genes[index], "Gene_B" = PPI_genes[index]))
		}
   
	}
	
	

	cat("Performing cis- genetic alteration CSEA analysis. It may take a long time.\n", file = stdout())

	if(arguments$core < 2){

		for(index in 1 : nrow(params)){
			result <- CSEA(params[index, "index"], params[index, "Gene_A"], params[index, "Gene_B"])
			if(nrow(result) > 0){
		 		CSEA_results <- rbind(CSEA_results, result)
		 	}
		 	if((index %% 100) == 0){
		 		cat(paste0("#", index, "/", nrow(params), "\n"), file = stdout())
		 	}
		}

	}else{
		if(arguments$core > detectCores()){
			stop("Error: specified CPU cores are larger than available!")
		}
		
		options(cl.cores = arguments$core)
		this.cluster <- makeCluster(getOption("cl.cores", 2))
		clusterCall(cl = this.cluster, fun = function(){
		  library(data.table);
		  library(fgsea)
		})
		if(formalArgs(clusterExport)[2] %in% "list"){
			clusterExport(cl = this.cluster, list = c("CSEA", "params", "exp_tpm_dt", "mutations_dt", "arguments"))
		}else{
			clusterExport(cl = this.cluster, varlist = c("CSEA", "params", "exp_tpm_dt", "mutations_dt", "arguments"))
		}

		CSEA_results1 <- 
		parLapply(cl = this.cluster,
			1 : nrow(params),
			function(idx) {
			  args1 <- as.list(params[idx,])
			  formals(CSEA) <- args1
			  rtn <- replicate(1, CSEA(), simplify = FALSE)
			  return(rtn)
			}
		)
		stopCluster(this.cluster)

		for(j in 1:length(CSEA_results1)){
			CSEA_results <- rbind(CSEA_results, CSEA_results1[[j]][[1]])
		}

	}
	
	if(nrow(CSEA_results) > 0){
		CSEA_results <- CSEA_results[!duplicated(CSEA_results[, c("Gene_A", "Gene_B", "Mutation_genome_change")]),]
	}
	
	
	if(arguments$remove_trans_results_from_cis_sig){
		CSEA_results$P_val_adjusted <- p.adjust(CSEA_results$P_val, method = "BH")
		omit_genes_for_trans_detection <- as.character(CSEA_results$Gene_A[which(CSEA_results$P_val_adjusted < 0.1 & CSEA_results$Gene_A == CSEA_results$Gene_A)])
		CSEA_results$P_val_adjusted <- NULL
		
	}
	
	cat("cis- genetic alteration CSEA analysis complete.\n", file = stdout())
}


params <- NULL
params_tmp <- NULL
# trans


cat("Construct a parameter table for trans- genetic alteration CSEA analysis.\n", file = stdout())
if(arguments$start_PPI > 0 & arguments$end_PPI > 0){
	
	for(index in 1 : row_ppi){
	  
		if((index %% 10000) == 0) cat(paste0("#", index, "/", row_ppi, "\n"), file = stdout())
	  
	  if(length(mutations_dt$gene_name %in% Input.net[index, 1]) > 0){ # try gene A as mutation gene
		 	if(!is.na(omit_genes_for_trans_detection)){
			 	if(!Input.net[index, 1] %in% omit_genes_for_trans_detection){
				 	params_tmp <- rbind(params_tmp, data.frame("index" = (index + arguments$start_PPI - 1), "Gene_A" = Input.net[index, 1], "Gene_B" = Input.net[index, 2]))
				    
			 	}
		 	}else{
		 		params_tmp <- rbind(params_tmp, data.frame("index" = (index + arguments$start_PPI - 1), "Gene_A" = Input.net[index, 1], "Gene_B" = Input.net[index, 2]))
		 	}

		}
		  
		if(length(mutations_dt$gene_name %in% Input.net[index, 2]) > 0){ # try gene B as mutation gene
	    if(!is.na(omit_genes_for_trans_detection)){
		    if(!Input.net[index, 2] %in% omit_genes_for_trans_detection){
			  	params_tmp <- rbind(params_tmp, data.frame("index" = (index + arguments$start_PPI - 1), "Gene_A" = Input.net[index, 2], "Gene_B" = Input.net[index, 1]))	
				}
			}else{
				params_tmp <- rbind(params_tmp, data.frame("index" = (index + arguments$start_PPI - 1), "Gene_A" = Input.net[index, 2], "Gene_B" = Input.net[index, 1]))
			}
		}
		  
		if((index %% 1000) == 0){
			params <- rbind(params, params_tmp)
			params_tmp <- NULL
		}
	}
	params <- rbind(params, params_tmp)
	
}else{
	params1 <- data.frame("index" = rownames(Input.net), "Gene_A" = Input.net[,1], "Gene_B" = Input.net[,2])
	params2 <- data.frame("index" = rownames(Input.net), "Gene_A" = Input.net[,2], "Gene_B" = Input.net[,1])
	params <- rbind(params1, params2)
	params <- params[which(params$Gene_A %in% mutations_dt$gene_name),]
	if(!is.na(omit_genes_for_trans_detection)){
		params <- params[which(!(params$Gene_A %in% omit_genes_for_trans_detection)),]
	}
	rm(params1)
	rm(params2)
}

cat("Performing trans- genetic alteration CSEA analysis. It may take a long time.\n", file = stdout())
if(arguments$core < 2){

	for(index in 1 : nrow(params)){
		result <- CSEA(params[index, "index"], params[index, "Gene_A"], params[index, "Gene_B"])
		if(nrow(result) > 0){
	    CSEA_results <- rbind(CSEA_results, result)
		}
		if((index %% 1000) == 0){
			cat(paste0("#", index, "/", nrow(params), "\n"), file = stdout())
		}
	}

}else{
	if(arguments$core > detectCores()){
    stop("Error: specified CPU cores are larger than available!")
	}

	options(cl.cores = arguments$core)
	this.cluster <- makeCluster(getOption("cl.cores", 2))
	clusterCall(cl = this.cluster, fun = function(){
	 	library(data.table);
	 	library(fgsea)
	})
	if(formalArgs(clusterExport)[2] %in% "list"){
		clusterExport(cl = this.cluster, list = c("CSEA", "params", "exp_tpm_dt", "mutations_dt", "arguments"))
	}else{
		clusterExport(cl = this.cluster, varlist = c("CSEA", "params", "exp_tpm_dt", "mutations_dt", "arguments"))
	}

	CSEA_results1 <- 
	parLapply(cl = this.cluster,
		1 : nrow(params),
		function(idx) {
			args1 <- as.list(params[idx,])
			formals(CSEA) <- args1
			rtn <- replicate(1, CSEA(), simplify = FALSE)
			return(rtn)
		}
	)
	stopCluster(this.cluster)

	for(j in 1 : length(CSEA_results1)){
		CSEA_results <- rbind(CSEA_results, CSEA_results1[[j]][[1]])
	}

}

if(nrow(CSEA_results) > 0){
	CSEA_results <- CSEA_results[!duplicated(CSEA_results[, c("Gene_A", "Gene_B", "Mutation_genome_change")]),]

	CSEA_results$P_val_adjusted <- p.adjust(CSEA_results$P_val, method = "BH")
	CSEA_results$ID <- (1 : nrow(CSEA_results))

	CSEA_results <- CSEA_results[, c("ID", "index", "Gene_A", "Gene_B", "Mutation_type", "Mutation_genome_change", "Mutation_protein_change", "ES", "NES", "P_val", "cell_line", "P_val_adjusted")]
		
}
	
cat("trans- genetic alteration CSEA analysis complete.\n", file = stdout())


arguments$Step1_result_filename <- paste0(arguments$prefix, "analysis_results_Step1.tsv")

if(is.null(nrow(CSEA_results))){
	write.table(CSEA_results, file = paste0("./", arguments$folder, "/", arguments$Step1_result_filename), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}else if(nrow(CSEA_results) > 0){
	CSEA_results <- CSEA_results[!duplicated(CSEA_results[, c("Gene_A", "Gene_B", "Mutation_genome_change")]),]

	CSEA_results$P_val_adjusted <- p.adjust(CSEA_results$P_val, method = "BH")
	CSEA_results$ID <- (1 : nrow(CSEA_results))

	CSEA_results <- CSEA_results[, c("ID", "index", "Gene_A", "Gene_B", "Mutation_type", "Mutation_genome_change", "Mutation_protein_change", "ES", "NES", "P_val", "cell_line", "P_val_adjusted")]

	if(arguments$full_result){
		write.table(CSEA_results, file = paste0("./", arguments$folder, "/", arguments$Step1_result_filename), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}else{
		write.table(CSEA_results[CSEA_results$P_val_adjusted < 0.1,], file = paste0("./", arguments$folder, "/", arguments$Step1_result_filename), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
}

rm(params)
rm(params_tmp)

cat("Step 1 complete.\n\n", file = stdout())

######### Step 2: mediation analysis #########

cat("Step 2: reading data\n", file = stdout())

# Step 1: Load omics data
arguments$Step2_result_filename <- paste0(arguments$prefix, "analysis_results_Step2.tsv")


if(file.exists(paste0("./", arguments$folder, "/", arguments$Step1_result_filename))){
	CSEA_results <- fread(paste0("./", arguments$folder, "/", arguments$Step1_result_filename), header = TRUE, sep = "\t")
	
 	if(is.unsorted(CSEA_results$index)){
		CSEA_results <- CSEA_results[order(CSEA_results$index),]
		CSEA_results <- CSEA_results[!duplicated(CSEA_results[, c("Gene_A", "Gene_B", "Mutation_genome_change")]),]
 		CSEA_results$ID <- seq(1, nrow(CSEA_results), by = 1)
 	}
}else{
	stop(paste0("Step 2 relies on results from Step 1, please perform Step 1 CSEA analysis first! The provided file is ", paste0("./", arguments$folder, "/", arguments$Step1_result_filename), "\n"))
	
}

dependency_dt <- NULL
exp_tpm_dt <- NULL

if(file.exists(arguments$dep)){
  dependency_dt <- fread(arguments$dep, header = TRUE, sep = "\t")
  dependency_dt <- dependency_dt[, c("gene", "dependency", "cell_line")]
  setkey(dependency_dt, gene)
  
  if(nrow(CSEA_results) > 0){
    dependency_dt <- dependency_dt[.(unique(c(CSEA_results$Gene_A, CSEA_results$Gene_B)))]
  }
  setkey(dependency_dt, gene)
}else{
	stop(paste0("File path ", arguments$dep, " is not correctly specified."))
}

cat("Read", paste0(arguments$dep, " with ", nrow(dependency_dt), " rows after filtering.\n"), file = stdout())


if(file.exists(arguments$exp)){
  exp_tpm_dt <- fread(arguments$exp, header = TRUE, sep = "\t")
  if(!("cell_line" %in% colnames(exp_tpm_dt)) & !("-" %in% colnames(exp_tpm_dt))){
    stop(paste0("Data ", arguments$exp, " does not contain columns of cell_line or disease."))
  }
  setkey(exp_tpm_dt, gene)
  
  if(nrow(CSEA_results) > 0){
    exp_tpm_dt <- exp_tpm_dt[.(unique(c(CSEA_results$Gene_A, CSEA_results$Gene_B)))]
  }
  setkey(exp_tpm_dt, gene)
}else{
	stop(paste0("File path ", arguments$exp, " is not correctly specified."))
}
cat("Read", paste0(arguments$exp, " with ", nrow(exp_tpm_dt), " rows after filtering.\n"), file = stdout())


# 2.2 Mediation

lm_fit_fast <- function(dat, X, out) {
	tryCatch({
		beta = solve(crossprod(X), t(crossprod(dat, X)))
		return(beta[out, 1])
	},
	error = function(e){
		return(NA)
	})
}


# Freedman Lane Function of smaller p-value

freedman_lane_sim_refine_gamma1_p_value <- function(adj.design, adj.residual, adj.reduced.fitted, adj.reduced.residual, nperms, adj.stat, gamma1.orig) {
	n.greater.gamma1 <- 0
	# number of simulations with a stat greater than the observed.
	nc <- ncol(adj.reduced.residual)

	adj.residual <- t(adj.residual)


	adj.reduced.residual <- t(adj.reduced.residual)
	adj.reduced.fitted <- t(adj.reduced.fitted)

	for(i in 1:nperms){

		ystar.null <- NULL
		gamma1.null <- NULL
		
    set.seed(i)
		perm <- sample(1:nc)
  
		for(j in 1:nc){
	    ystar.null[j] <- adj.reduced.fitted[j] + adj.reduced.residual[perm[j]]
		}
  
		gamma1.star.null <- lm_fit_fast(as.matrix(ystar.null), adj.design, adj.stat)
		
		if(length(gamma1.star.null) > 0 & !is.na(gamma1.star.null)){
			n.greater.gamma1 <- n.greater.gamma1 + (abs(gamma1.star.null) >= abs(gamma1.orig))
		}else{
			n.greater.gamma1 <- n.greater.gamma1 + 1
		}
		
	}

	return(n.greater.gamma1)

}

freedman_lane_sim_refine_delta_p_value <- function(delta.design,
	   delta.residual, nperms, delta.stat, delta.orig) {
	n.greater.delta <- 0
	# number of simulations with a stat greater than the observed.
	nc <- ncol(delta.residual)

	delta.residual <- t(delta.residual)

	for(i in 1 : nperms){

		deltastar.null <- NULL
		
	 	set.seed(i)
		perm <- sample(1 : nc)
  
		for(j in 1:nc){ 
  		deltastar.null[j] <- delta.residual[perm[j]]
		}
  
		delta.star.null <- lm_fit_fast(as.matrix(deltastar.null), delta.design, delta.stat)

	 	if(length(delta.star.null) > 0 & !is.na(delta.star.null)){
			n.greater.delta <- n.greater.delta + (abs(delta.star.null) >= abs(delta.orig))
	 	}else{
		 	n.greater.delta <- n.greater.delta + 1
		}
	}

	return(n.greater.delta)
}

mediation_analysis <- 
function(dat = NULL, nperm = 100, Permut1 = FALSE, precise_pvalue = FALSE, cis = FALSE){

	y <- dat$y
	x1 <- dat$x1
	x2 <- dat$x2
	m <- dat$m

	# model fitting
	crude.lr <- lm(y ~ x2 + x1, data = dat, model = FALSE)
	adj.lr  <- update(crude.lr, formula = . ~ m + .)
	ie.lr  <- lm(m ~ x1 + x2, data = dat, model = FALSE)
	x.lr <- lm(x2 ~ x1, data = dat, model = FALSE)

	# Data needs to be in matrix form for the following functions
	dat <- dat[, c(4, 1, 2, 3)]
	my.matrix <- as.matrix(dat)
	y <- as.matrix(my.matrix[, 1]) # y
	adj.variables_null <- as.matrix(my.matrix[, 2:3]) # x1, x2
	adj.variables <- as.matrix(my.matrix[, 2:4]) # x1, x2, m
	m <- as.matrix(my.matrix[, 4]) # m
	ie.variables_null <- as.matrix(my.matrix[, 3]) # x2
	ie.variables  <- as.matrix(my.matrix[, 2:3]) # x1, x2
	x1 <- as.matrix(my.matrix[, 2]) # x1
	x2 <- as.matrix(my.matrix[, 3]) # x2
	gamma2.variables <- as.matrix(my.matrix[, 3:4]) # x2, m
	gamma3.variables <- as.matrix(my.matrix[, c(2,4)]) # x1, m

	colnames(x1) <- "x1" # X
	colnames(y) <- "y" # Y
	colnames(m) <- "m" # M
	colnames(ie.variables_null) <- "x2" # C

	# Set design matrices
	# gamma 1
	adj.design <- cbind(1, adj.variables) # 1 x1, x2, m
	adj.reduced.design <- cbind(1, adj.variables_null) # 1 x1, x2

	# alpha 1
	ie.design <- cbind(1, ie.variables) # 1 x1, x2
	ie.reduced.design <- cbind(1, ie.variables_null) # 1 x2

	# alpha 2
	alpha2.design <- cbind(1, ie.variables) # 1 x1, x2
	alpha2.reduced.design <- cbind(1, x1) # 1 x1

	# gamma 2
	gamma2.design <- cbind(1, adj.variables) # 1 x1, x2, m
	gamma2.reduced.design <- cbind(1, gamma2.variables) # 1 x2, m

	# gamma 3
	gamma3.design <- cbind(1, adj.variables) # 1 x1, x2, m
	gamma3.reduced.design <- cbind(1, gamma3.variables) # 1 x1, m
	
	# delta
	delta.design <- cbind(1, x1) # 1 x1

	# Obtain fitted outcomes and residuals
	# gamma 1
	adj.fit <- lm.fit(adj.design, y)
	adj.fitted <- t(adj.fit$fitted.values)
	adj.residual <- t(adj.fit$residuals)

	adj.fit0 <- lm.fit(adj.reduced.design, y)
	adj.reduced.fitted <- t(adj.fit0$fitted.values)
	adj.reduced.residual <- t(adj.fit0$residuals)

	# alpha 1
	ie.fit  <- lm.fit(ie.design, m)
	ie.fitted <- t(ie.fit$fitted.values)
	ie.residual <- t(ie.fit$residuals)

	ie.fit0 <- lm.fit(ie.reduced.design, m)
	ie.reduced.fitted <- t(ie.fit0$fitted.values)
	ie.reduced.residual <- t(ie.fit0$residuals)

	# alpha 2
	alpha2.fit  <- lm.fit(alpha2.design, m)
	alpha2.fitted <- t(alpha2.fit$fitted.values)
	alpha2.residual <- t(alpha2.fit$residuals)

	alpha2.fit0 <- lm.fit(alpha2.reduced.design, m)
	alpha2.reduced.fitted <- t(alpha2.fit0$fitted.values)
	alpha2.reduced.residual <- t(alpha2.fit0$residuals)

	#gamma 2
	gamma2.fit <- lm.fit(gamma2.design, y)
	gamma2.fitted <- t(gamma2.fit$fitted.values)
	gamma2.residual <- t(gamma2.fit$residuals)

	gamma2.reduced.fit <- lm.fit(gamma2.reduced.design, y)
	gamma2.reduced.fitted <- t(gamma2.reduced.fit$fitted.values)
	gamma2.reduced.residual <- t(gamma2.reduced.fit$residuals)

	#gamma 3
	gamma3.fit <- lm.fit(gamma3.design, y)
	gamma3.fitted <- t(gamma3.fit$fitted.values)
	gamma3.residual <- t(gamma3.fit$residuals)

	gamma3.reduced.fit <- lm.fit(gamma3.reduced.design, y)
	gamma3.reduced.fitted <- t(gamma3.reduced.fit$fitted.values)
	gamma3.reduced.residual <- t(gamma3.reduced.fit$residuals)

	#delta
	delta.fit <- lm.fit(delta.design, x2)
	delta.fitted <- t(delta.fit$fitted.values)
	delta.residual <- t(delta.fit$residuals)

	# Calculate test statistic
	adj.stat <- setdiff(colnames(adj.variables), colnames(adj.variables_null))
	gamma1.orig <- adj.fit$coef[adj.stat] #gamma 1

	ie.stat	<- setdiff(colnames(ie.variables), colnames(ie.variables_null))
	# ie.stat	<- "x1"
	alpha1.orig 	<- ie.fit$coef[ie.stat] # alpha 1

	alpha2.stat	<- setdiff(colnames(ie.variables), colnames(x1))
	alpha2.orig 	<- alpha2.fit$coef[alpha2.stat] # alpha 2

	gamma2.stat <- setdiff(colnames(adj.variables), colnames(gamma2.variables)) # x2
	gamma2.orig 	<- ie.fit$coef[gamma2.stat] # gamma 2

	gamma3.stat <- setdiff(colnames(adj.variables), colnames(gamma3.variables)) # x1
	gamma3.orig 	<- ie.fit$coef[gamma3.stat] # gamma 3

	delta.stat <- colnames(x1) # x1
	delta.orig 	<- delta.fit$coef[delta.stat] # delta

	med.orig.lr <- alpha1.orig * gamma1.orig

	# Permutations 
	n.greater.gamma1 <- 0
	n.greater.delta	<- 0
	p.gamma1.null <- -1
	p.delta.null <- -1

	if(Permut1){
		if(cis){
			n.greater.delta <- freedman_lane_sim_refine_delta_p_value(delta.design, delta.residual, nperm, delta.stat, delta.orig)

			p.delta.null <- n.greater.delta / nperm

			if(precise_pvalue & p.delta.null < 0.1){
			 
				refined_p_delta <- NA
				refined_N_delta <- freedman_lane_sim_refine_delta_p_value(delta.design, delta.residual, nperm * 100, delta.stat, delta.orig)
				refined_p_delta <- refined_N_delta / (nperm * 100)
				
				if(refined_p_delta < 0.001){
					refined_N_delta <- freedman_lane_sim_refine_delta_p_value(delta.design, delta.residual, nperm * 1000, delta.stat, delta.orig)
					refined_p_delta <- refined_N_delta / (nperm * 1000)
				}
				if(refined_p_delta > 0){
					p.delta.null <- refined_p_delta
				}
				
			}
		}else{
			n.greater.gamma1 <- freedman_lane_sim_refine_gamma1_p_value(adj.design, adj.residual, adj.reduced.fitted, adj.reduced.residual, nperm, adj.stat, gamma1.orig)

			p.gamma1.null <- n.greater.gamma1 / nperm

			if(precise_pvalue & p.gamma1.null < 0.1){
			 
				refined_p_gamma1 <- NA
				refined_N_gamma1 <- freedman_lane_sim_refine_gamma1_p_value(adj.design, adj.residual, adj.reduced.fitted, adj.reduced.residual, nperm * 100, adj.stat, gamma1.orig)
				refined_p_gamma1 <- refined_N_gamma1 / (nperm * 100)
				if(refined_p_gamma1 < 0.001){
					refined_N_gamma1 <- freedman_lane_sim_refine_gamma1_p_value(adj.design, adj.residual, adj.reduced.fitted, adj.reduced.residual, nperm * 1000, adj.stat, gamma1.orig)
					refined_p_gamma1 <- refined_N_gamma1 / (nperm * 1000)
				}
				if(refined_p_gamma1 > 0){
					p.gamma1.null <- refined_p_gamma1
				}
		  
			}
		}
	 
	}

	tmp<-data.frame("gamma_1" = round(coef(adj.lr)[2], digits = 6), "gamma_2" = round(coef(adj.lr)[3], digits = 6), "gamma_3" = round(coef(adj.lr)[4], digits = 6), "alpha_1" = round(coef(ie.lr)[2], digits = 6), "alpha_2" = round(coef(ie.lr)[3], digits = 6), "delta" = round(coef(x.lr)[2], digits = 6), "P_gamma_1" = p.gamma1.null, "P_delta" = p.delta.null)

	return(tmp)
}

mediation_analysis_exp_dep <-
function(ID = 1, Gene_A = "", Gene_B = "", Cell_lines = ""){
	
	nperm <- arguments$nperm
	Permut1 <- arguments$permutation
	precise_pvalue <- arguments$precise_mediation_pvalue
	use_all_cell_lines <- arguments$use_all_cell_lines

	res <- list("result" = data.frame("ID" = ID, "Gene_A" = Gene_A, "Gene_B" = Gene_B, "gamma_1" = NA, "gamma_2" = NA, "gamma_3" = NA, "alpha_1" = NA, "alpha_2" = NA, "delta" = NA, "P_gamma_1" = NA, "P_delta" = NA))

	cls_split <- strsplit(Cell_lines, split = ",")[[1]];

	expr_dat_gene_1 <- exp_tpm_dt[.(Gene_A)]
	expr_dat_gene_2 <- exp_tpm_dt[.(Gene_B)]
	if(nrow(expr_dat_gene_1) < 1 | nrow(expr_dat_gene_2) < 1){
		return(res)
	}

	if(use_all_cell_lines){

	}else{
		expr_dat_gene_1 <- expr_dat_gene_1[which(expr_dat_gene_1$cell_line %in% cls_split),]
		expr_dat_gene_2 <- expr_dat_gene_2[which(expr_dat_gene_2$cell_line %in% cls_split),]
	}

	if(nrow(expr_dat_gene_1) < 1 | nrow(expr_dat_gene_2) < 1){
		return(res)
	}

	expr_dat_gene_1 <- expr_dat_gene_1[, c("rna_expression", "cell_line")]
	colnames(expr_dat_gene_1) <- c("rna_expression_1", "cell_line")
	expr_dat_gene_2 <- expr_dat_gene_2[, c("rna_expression", "cell_line")]
	colnames(expr_dat_gene_2) <- c("rna_expression_2", "cell_line")


	dep_dat_gene_1 <- dependency_dt[.(Gene_A)]
	if(use_all_cell_lines){

	}else{
		dep_dat_gene_1 <- dep_dat_gene_1[which(dep_dat_gene_1$cell_line %in% cls_split),]
	}
	if(nrow(dep_dat_gene_1) < 1){
		return(res)
	}

	if(ncol(dep_dat_gene_1) > 3){
		colnames(dep_dat_gene_1) <- c("gene_1", "dependency_1", "cell_line", "context")
	}else{
		colnames(dep_dat_gene_1) <- c("gene_1", "dependency_1", "cell_line")
	}


	setkey(dep_dat_gene_1, 'cell_line')
	setkey(expr_dat_gene_1, 'cell_line')

	dep_dat_gene_2 <- dependency_dt[.(Gene_B)]
	if(use_all_cell_lines){

	}else{
		dep_dat_gene_2 <- dep_dat_gene_2[which(dep_dat_gene_2$cell_line %in% cls_split),]
	}
	if(nrow(dep_dat_gene_2) < 1){
		return(res)
	}
	if(ncol(dep_dat_gene_2) > 3){
		colnames(dep_dat_gene_2) <- c("gene_2", "dependency_2", "cell_line", "context")
	}else{
		colnames(dep_dat_gene_2) <- c("gene_2", "dependency_2", "cell_line")
	}


	setkey(dep_dat_gene_2, 'cell_line')
	setkey(expr_dat_gene_2, 'cell_line')

	DTlist <- list(expr_dat_gene_1, dep_dat_gene_1, expr_dat_gene_2, dep_dat_gene_2)
	merged_dat_dep_exp_gene_1_2 <-Reduce(function(X,Y) X[Y], DTlist)
	merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2[!is.na(merged_dat_dep_exp_gene_1_2$cell_line),]
	merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2 %>% filter(!is.na(rna_expression_1) & !is.na(rna_expression_2) & !is.na(dependency_1) & !is.na(dependency_2))
	if(nrow(merged_dat_dep_exp_gene_1_2) < 3){
		return(res)
	}

	# data merged_dat_dep_exp_gene_1_2
	# col 1: rna_expression # gene 1's expr
	# col 2: dependency # gene 1's dep
	# col 3: i.rna_expression # gene 2's expr
	# col 4: i.dependency # gene 2's dep

	merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2[, c("rna_expression_1", "dependency_1", "rna_expression_2", "dependency_2")]
	colnames(merged_dat_dep_exp_gene_1_2) <- c("x1", "x2", "m", "y")
	cis <- FALSE
	if(Gene_A == Gene_B){
		merged_dat_dep_exp_gene_1_2$m <- rep(0, nrow(merged_dat_dep_exp_gene_1_2))
		merged_dat_dep_exp_gene_1_2$y <- rep(0, nrow(merged_dat_dep_exp_gene_1_2))
		cis <- TRUE
	}
	if(nrow(merged_dat_dep_exp_gene_1_2) > arguments$max_size){
		set.seed(arguments$seed)
		merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2[sample(1:arguments$max_size),]

	}
		

	if(nrow(merged_dat_dep_exp_gene_1_2) >= arguments$min_size){ # lr requires at least five records
		result_gene1_to_2 <- mediation_analysis(merged_dat_dep_exp_gene_1_2, nperm, Permut1, precise_pvalue, cis)
		if(result_gene1_to_2$P_gamma_1 < 0){
		 	result_gene1_to_2$P_gamma_1 <- NA
		}
		if(result_gene1_to_2$P_delta < 0){
		 	result_gene1_to_2$P_delta <- NA
		}
		res <- list("result" = data.frame("ID" = ID, "Gene_A" = Gene_A, "Gene_B" = Gene_B, "gamma_1" = result_gene1_to_2[1], "gamma_2" = result_gene1_to_2[2], "gamma_3" = result_gene1_to_2[3], "alpha_1" = result_gene1_to_2[4], "alpha_2" = result_gene1_to_2[5], "delta" = result_gene1_to_2[6], "P_gamma_1" = result_gene1_to_2[7], "P_delta" = result_gene1_to_2[8]))

		return(res)
	}else{
	 	return(res)
	}
}

final_med_results <- data.frame("ID" = vector("numeric", length = 0), "Gene_A" = vector("character", length = 0), "Gene_B" = vector("character", length = 0), "gamma_1" = vector("numeric", length = 0), "gamma_2" = vector("numeric", length = 0), "gamma_3" = vector("numeric", length = 0), "alpha_1" = vector("numeric", length = 0), "alpha_2" = vector("numeric", length = 0), "delta" = vector("numeric", length = 0), "P_gamma_1" = vector("numeric", length = 0), "P_delta" = vector("numeric", length = 0))

cat("Construct a parameter table for mediation analysis.\n", file = stdout())

records <- data.frame("ID" = CSEA_results[, "ID"], CSEA_results[, "Gene_A"], CSEA_results[, "Gene_B"], "Cell_lines" = CSEA_results[, "cell_line"])
colnames(records) <- c("ID", "Gene_A", "Gene_B", "Cell_lines")


cat("Perform mediation analysis. It may take a long time.\n", file = stdout())

if(arguments$core < 2){
  tmp_results <- NULL
  
	for(idx in 1 : nrow(records)){
		med <- mediation_analysis_exp_dep(records[idx, 1], records[idx, 2], records[idx, 3], records[idx, 4])
		tmp_results <- rbind(tmp_results, med[["result"]])
		if((idx %% 100) == 0){
		 	# at most 100 records are stored in tmp_results to keep the data frame not too big
		 	final_med_results <- rbind(final_med_results, tmp_results)
	  
		 	tmp_results <- NULL
		}
		if((idx %% (floor(nrow(records)/100))) == 0) cat(paste0("#", idx, "/", nrow(records), "\n"), file = stdout())
	}
	final_med_results <- rbind(final_med_results, tmp_results)
}else{
  if(arguments$core > detectCores()){
    stop("Error: specified CPU cores are larger than available!")
  }
 
  options(cl.cores = arguments$core)
  this.cluster <- makeCluster(getOption("cl.cores", 2))
  clusterCall(cl = this.cluster, fun = function(){
    library(dplyr);
    library(data.table)
  })
  
	if(formalArgs(clusterExport)[2] %in% "list"){
    clusterExport(cl = this.cluster, list = c("arguments", "records", "mediation_analysis_exp_dep", "mediation_analysis", "freedman_lane_sim", "lm_fit_fast", "dependency_dt", "freedman_lane_sim_refine_gamma1_p_value", "freedman_lane_sim_refine_delta_p_value", "exp_tpm_dt"))
  }else{
    clusterExport(cl = this.cluster, varlist = c("arguments", "records", "mediation_analysis_exp_dep", "mediation_analysis", "freedman_lane_sim", "lm_fit_fast", "dependency_dt", "freedman_lane_sim_refine_gamma1_p_value", "freedman_lane_sim_refine_delta_p_value", "exp_tpm_dt"), envir = environment())
  }
  
  mediation_analysis_results <- 
	  parLapply(cl = this.cluster,
	  1:nrow(records),
	  function(idx) {
	    rtn <- mediation_analysis_exp_dep(records[idx, 1], records[idx, 2], records[idx, 3], records[idx, 4])
	    return(rtn)
	  }
  )
  
  stopCluster(this.cluster)
  
  for( j in 1:length(mediation_analysis_results)){
		final_med_results <- rbind(final_med_results, mediation_analysis_results[[j]][["result"]])
  }
}

write.table(final_med_results, file = paste0("./", arguments$folder, "/", arguments$Step2_result_filename), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Step 2 completed.\n\n", file = stdout())

######### Step 3: merge data #########


cat("Step 3: Merging data.\n", file = stdout())
######### Step 3.1: Load data #########
if(!exists("final_med_results")){
	final_med_results <- fread(paste0("./", arguments$folder_Step2, "/", arguments$Step2_result_filename), header = TRUE, sep = "\t")
	strcol <- c("ID", "gamma_1", "gamma_2", "gamma_3", "alpha_1", "alpha_2", "delta", "P_gamma_1", "P_delta")
	for(col in strcol){
		set(final_med_results, j = col, value = as.numeric(final_med_results[[col]]))
	}
}

if(!exists("CSEA_results")){
	CSEA_results <- fread(paste0("./", arguments$folder_Step1, "/", arguments$Step1_result_filename), header = TRUE, sep = "\t")
}


merged_CSEA_Mediation <- merge(CSEA_results, final_med_results, by = 'ID', all.x = T)
colnames(merged_CSEA_Mediation)[3] <- "Gene_A"
colnames(merged_CSEA_Mediation)[4] <- "Gene_B"
merged_CSEA_Mediation$Gene_A.y <- NULL
merged_CSEA_Mediation$Gene_B.y <- NULL
  
if(length(which(duplicated(CSEA_results[, c("Gene_A", "Gene_B", "Mutation_genome_change")]))) > 0){
  
	merged_CSEA_Mediation <- merged_CSEA_Mediation[which(!duplicated(CSEA_results[, c("Gene_A", "Gene_B", "Mutation_genome_change")])),]
	cat(paste0(length(which(duplicated(CSEA_results[, c("Gene_A", "Gene_B", "Mutation_genome_change")]))), " records are duplicated and will be deleted.\n"), file = stdout())
}

if(is.null(merged_CSEA_Mediation)){
	stop(paste0("No record is available for filtering."))
}

rm(CSEA_results)
rm(final_med_results)

original_nrow <- nrow(merged_CSEA_Mediation)
merged_CSEA_Mediation <- merged_CSEA_Mediation[!(is.na(merged_CSEA_Mediation$P_delta) & is.na(merged_CSEA_Mediation$P_gamma_1)),]
cat(paste0(nrow(merged_CSEA_Mediation), " records with non-NA cis-/trans-mediation P-values from ", original_nrow , " records were kept.\n"), file = stdout())


merged_CSEA_Mediation_cis <- merged_CSEA_Mediation[which(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B), c("ID", "Gene_A", "Gene_B", "Mutation_genome_change", "NES", "P_val_adjusted", "delta")]

merged_CSEA_Mediation_trans <- merged_CSEA_Mediation[which(merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B), c("ID", "Gene_A", "Gene_B", "Mutation_genome_change", "NES", "P_val_adjusted", "gamma_1")]

######### Step 3.2: Filtering data #########
cat("Step 3: filtering trans-mutations affected by cis confounding effects.\n", file = stdout())
removed_trans_line_index <- NULL
	
cis_values <- function(ID = NULL) {
	id = ID
	tmp <- merged_CSEA_Mediation_trans[which(merged_CSEA_Mediation_trans$ID == id), c("ID", "Gene_A", "Gene_B", "NES", "P_val_adjusted")]

	index_select <- which(merged_CSEA_Mediation_cis$Gene_B == tmp$Gene_B[1])

	if(length(index_select) < 1){
		res <- data.frame("ID" = ID, "Gene_A" = tmp$Gene_A, "Gene_B" = tmp$Gene_B, "NES_cis" = 0, "CSEA_pvaladj_cis" = 1, "NES" = tmp$NES, "P_val_adjusted" = tmp$P_val_adjusted)
		return(res)
	}
	NES <- merged_CSEA_Mediation_cis$NES[index_select]
	pvaladj_tmp <- merged_CSEA_Mediation_cis$P_val_adjusted[index_select]
	if(length(NES) < 1){
		NES <- NA
	}else{
		if(length(NES) > 1){
			NES <- NES[which.min(merged_CSEA_Mediation_cis$P_val)]
			pvaladj_tmp <- pvaladj_tmp[which.min(merged_CSEA_Mediation_cis$P_val)]
			if(length(NES) > 1){
				NES <- NES[1]
				pvaladj_tmp <- pvaladj_tmp[1]
			}
		}

	}

	res <- data.frame("ID" = ID, "Gene_A" = tmp$Gene_A, "Gene_B" = tmp$Gene_B, "NES_cis" = NES, "CSEA_pvaladj_cis" = pvaladj_tmp, "NES" = tmp$NES, "P_val_adjusted" = tmp$P_val_adjusted)
	return(res)
}

if(arguments$remove_cis_confounding){
	final_cis_results <- NULL
	if(arguments$core < 2){
		params <- merged_CSEA_Mediation_trans$ID[1 : nrow(merged_CSEA_Mediation_trans)]
		cis_results <- lapply(params, cis_values)
			 
	}else{
		options(cl.cores = arguments$core)
		this.cluster <- makeCluster(getOption("cl.cores", 2))
		clusterCall(cl = this.cluster, fun = function(){
			library(data.table)
		})
		if(formalArgs(clusterExport)[2] %in% "list"){
			clusterExport(cl = this.cluster, list = c("cis_values", "merged_CSEA_Mediation_trans", "merged_CSEA_Mediation_cis"))
		}else{
			clusterExport(cl = this.cluster, varlist = c("cis_values", "merged_CSEA_Mediation_trans", "merged_CSEA_Mediation_cis"))
		}

		cis_results <- 
		parSapply(cl = this.cluster,
		1 : nrow(merged_CSEA_Mediation_trans), # nrow(merged_CSEA_Mediation_trans)
		function(idx) {
		 	args1 <- as.list(merged_CSEA_Mediation_trans[idx, "ID"])
		 	formals(cis_values) <- args1
		 	rtn <- replicate(1, cis_values(), simplify = FALSE)
		 	return(rtn)
		})
		stopCluster(this.cluster)

	}

	final_cis_results <- data.frame(matrix(unlist(cis_results), nrow=length(cis_results), byrow = TRUE))
	colnames(final_cis_results) <- c("ID", "Gene_A", "Gene_B", "NES_cis", "CSEA_pvaladj_cis", "NES", "P_val_adjusted")
	final_cis_results$NES_cis <- as.numeric(final_cis_results$NES_cis)
	final_cis_results$NES <- as.numeric(final_cis_results$NES)
	final_cis_results$CSEA_pvaladj_cis <- as.numeric(final_cis_results$CSEA_pvaladj_cis)

	removed_trans_line_index <- final_cis_results$ID[which(final_cis_results$CSEA_pvaladj_cis <= arguments$cis_confounding_pvalue_cutoff & abs(final_cis_results$NES_cis) > abs(final_cis_results$NES))]
	if(length(removed_trans_line_index) > 0){
		cat(paste0((nrow(merged_CSEA_Mediation) - length(removed_trans_line_index)), " of ", nrow(merged_CSEA_Mediation), " records were kept without potential cis regulation confounding.\n"), file = stdout())
	}
	merged_CSEA_Mediation <- merged_CSEA_Mediation[!(merged_CSEA_Mediation$ID %in% removed_trans_line_index),]
  
}else{
	cat("Step 3: \"Filtering trans regulation results that can be explained by cis regulation\" is disabled.\n", file = stdout())
}


write.table(merged_CSEA_Mediation, file = paste0("./", arguments$folder, "/", arguments$prefix, "combined_analysis_results_Step3.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Step 3 completed.\n\n", file = stdout())


######### Step 4: filter results based on CNA data #########
if(is.null(merged_CSEA_Mediation)){
	if(file.exists(paste0("./", arguments$folder, "/", arguments$prefix, "combined_analysis_results_Step3.tsv"))){
		merged_CSEA_Mediation <- fread(paste0("./", arguments$folder, "/", arguments$prefix, "combined_analysis_results_Step3.tsv"), header = TRUE, sep = "\t")
		
	}else{
		stop(paste0("Step 3 results are not available."))
	}
}


cat("Step 4: reading data.\n", file = stdout())

CNA_dt <- NULL


######### Step 4.1: Filtering data #########

CNA_confounder_exp_glm_analysis <- function(ID = NA, Gene_A = "", Gene_B = "", Cell_lines = ""){
	# Gene_A: gene with mutations
	# Gene_B: gene with CNA and exp
	# Cell_lines: Cell line with specific mutations in record with ID

	Permut1 <- arguments$cna_permutation_number

	cls_split <- strsplit(Cell_lines, split = ",")[[1]];

	exp_dat_gene_2 <- exp_tpm_dt[.(Gene_B), c("gene", "rna_expression", "cell_line")]

	CNA_dat_gene_2 <- CNA_dt[.(Gene_B), c("log_copy_number", "cell_line")]
	colnames(CNA_dat_gene_2) <- c("CNA", "cell_line")

	merged_dat_CNA_exp_gene_2 <- merge(exp_dat_gene_2, CNA_dat_gene_2, by = "cell_line")
	merged_dat_CNA_exp_gene_2 <- merged_dat_CNA_exp_gene_2[!is.na(merged_dat_CNA_exp_gene_2$cell_line),]
	rm(exp_dat_gene_2)
	rm(CNA_dat_gene_2)

	mutations_dat <- mutations_dt[.(Gene_A), c("gene_name", "cell_line")]
	merged_dat_CNA_exp_mut <- merge(merged_dat_CNA_exp_gene_2, mutations_dat, by = "cell_line", all = TRUE)
	merged_dat_CNA_exp_mut$Mutation <- 0

	merged_dat_CNA_exp_mut$Mutation[merged_dat_CNA_exp_mut$cell_line %in% cls_split] <- 1
	mutations_dat <- merged_dat_CNA_exp_mut
	rm(merged_dat_CNA_exp_gene_2)

	if(nrow(mutations_dat) < 1){
		tmp_data_frame <- data.frame("ID" = ID, "P_value_CNA_confounding" = NA)
		return(tmp_data_frame)
	}

	if(arguments$cna_discretization){
		CN <- rep(NA, nrow(mutations_dat))
		CN[mutations_dat$CNA >= arguments$cna_amplification_cutoff] <- 1
		CN[mutations_dat$CNA <= arguments$cna_deletion_cutoff] <- -1
		CN[mutations_dat$CNA > arguments$cna_deletion_cutoff & mutations_dat$CNA < arguments$cna_amplification_cutoff] <- 0
		mutations_dat$CNA <- CN # log2CN to binary CN
	}

	if(nrow(mutations_dat) < 2){
		tmp_data_frame <- data.frame("ID" = ID, "P_value_CNA_confounding" = NA)
		return(tmp_data_frame)
	}
	coef_mut <- NULL
	coef_CN <- NULL
   
	tryCatch({
		model <- glm(rna_expression ~ CNA + Mutation, data = mutations_dat)
		coef_mut <- coef(summary(model))[, 1]["Mutation"]
		coef_CN <- coef(summary(model))[, 1]["CNA"]
  },
  error = function(e){
		tmp_data_frame <- data.frame("ID" = ID, "P_value_CNA_confounding" = NA)
		return(tmp_data_frame)
  })

  smaller_num <- 0
   
  P_permutation <- NA
  if(anyNA(coef_CN)){
	 	 tmp_data_frame <- data.frame("ID" = ID, "P_value_CNA_confounding" = NA)
     return(tmp_data_frame)
	}
	 
	if(is.null(coef_CN)){
		tmp_data_frame <- data.frame("ID" = ID, "P_value_CNA_confounding" = NA)
		return(tmp_data_frame)
	}

	if(coef_CN >= 0) {
		coef_mut_noCNA <- NULL
		tryCatch({
		  model_mut_noCNA <- glm(rna_expression ~ Mutation, data = mutations_dat)
		 	coef_mut_noCNA <- coef(summary(model_mut_noCNA))[, 1]["Mutation"]
		},
		error = function(e){
			tmp_data_frame <- data.frame("ID" = ID, "P_value_CNA_confounding" = NA)
			return(tmp_data_frame)
		})
		 
		decrease_rate <- abs((coef_mut_noCNA - coef_mut)/coef_mut) # https://nbisweden.github.io/excelerate-scRNAseq/session-normalization/confounding-factors.pdf
		# higher decrease_rate indicates more contribution of CNA to expression
		for(index in 1 : Permut1){
			set.seed(index)
			perm <- sample(1 : nrow(mutations_dat))
			merged_dat1 <- mutations_dat[, c("rna_expression", "CNA")]
			merged_dat1$Mutation <- mutations_dat$Mutation[perm]

			tryCatch({
				model <- glm(rna_expression ~ CNA + Mutation, data = merged_dat1)
				coef_mut <- coef(summary(model))[, 1]["Mutation"]
				model_mut_noCNA <- glm(rna_expression ~ Mutation, data = merged_dat1)
				coef_mut_noCNA <- coef(summary(model_mut_noCNA))[, 1]["Mutation"]
				decrease_rate_p <- (coef_mut_noCNA - coef_mut)/coef_mut
				if(abs(decrease_rate_p) <= decrease_rate){
				  smaller_num <- smaller_num + 1
				}
			},
			error = function(e){
			})

		}

		P_permutation <- smaller_num/Permut1

		if(P_permutation < 0.2){ # Further refine p-values if close to cutoff cna_adjust_Pvalue_cutoff
			for(index in (1 + Permut1) : (Permut1 * 2)){
				set.seed(index)
				perm <- sample(1 : nrow(mutations_dat))
				merged_dat1 <- mutations_dat[, c("rna_expression", "CNA")]
				merged_dat1$Mutation <- mutations_dat$Mutation[perm]
				tryCatch({
			    model <- glm(rna_expression ~ CNA + Mutation, data = merged_dat1)
			    coef_mut <- coef(summary(model))[, 1]["Mutation"]
			   
			    model_mut_noCNA <- glm(rna_expression ~ Mutation, data = merged_dat1)
			    coef_mut_noCNA <- coef(summary(model_mut_noCNA))[, 1]["Mutation"]
			    decrease_rate_p <- (coef_mut_noCNA - coef_mut)/coef_mut
			    if(abs(decrease_rate_p) <= decrease_rate){
			      smaller_num <- smaller_num + 1
			    }
				},
				error = function(e){
				})

			}
				
			P_permutation <- smaller_num/(Permut1 * 2)

		}

		if(P_permutation < 0.05){
			for(index in (1 + Permut1 * 2) : (Permut1 * 10)){
				set.seed(index)
				perm <- sample(1 : nrow(mutations_dat))
				merged_dat1 <- mutations_dat[, c("rna_expression", "CNA")]
				merged_dat1$Mutation <- mutations_dat$Mutation[perm]
				tryCatch({
				  model <- glm(rna_expression ~ CNA + Mutation, data = merged_dat1)
			    coef_mut <- coef(summary(model))[, 1]["Mutation"]
			   
			    model_mut_noCNA <- glm(rna_expression ~ Mutation, data = merged_dat1)
			    coef_mut_noCNA <- coef(summary(model_mut_noCNA))[, 1]["Mutation"]
			    decrease_rate_p <- (coef_mut_noCNA - coef_mut)/coef_mut
			    if(abs(decrease_rate_p) <= decrease_rate){
			      smaller_num <- smaller_num + 1
			    }
				},
				error = function(e){
				 	
				})

			 }
			 P_permutation <- smaller_num/(Permut1 * 10)
		}
		if(P_permutation < 0.005){
			for(index in (1 + Permut1 * 10) : (Permut1 * 20)){
				set.seed(index)
				perm <- sample(1 : nrow(mutations_dat))
				merged_dat1 <- mutations_dat[, c("rna_expression", "CNA")]
				merged_dat1$Mutation <- mutations_dat$Mutation[perm]
				tryCatch({
					model <- glm(rna_expression ~ CNA + Mutation, data = merged_dat1)
					coef_mut <- coef(summary(model))[, 1]["Mutation"]

					model_mut_noCNA <- glm(rna_expression ~ Mutation, data = merged_dat1)
					coef_mut_noCNA <- coef(summary(model_mut_noCNA))[, 1]["Mutation"]
					decrease_rate_p <- (coef_mut_noCNA - coef_mut)/coef_mut
					if(abs(decrease_rate_p) <= decrease_rate){
					  smaller_num <- smaller_num + 1
					}
				},
				error = function(e){
					 	
				})

			}
			P_permutation <- smaller_num/(Permut1 * 20)
		}
	}

	# P-value indicates significance of CNA contribution to expression
	tmp_data_frame <- data.frame("ID" = ID, "P_value_CNA_confounding" = P_permutation)
	  
	return(tmp_data_frame)
}


CNA_confounder_analysis_data_frame <- data.frame("ID" = rep(NA, nrow(merged_CSEA_Mediation)), "P_value_CNA_confounding" = rep(NA, nrow(merged_CSEA_Mediation)))
removed_index_CNA_confounder <- NULL
records <- NULL

if(arguments$cna_adjust){
   gene_A_filtered <- unique(merged_CSEA_Mediation$Gene_A)
	 gene_B_filtered <- unique(merged_CSEA_Mediation$Gene_B)
	 cat("Step 4.1: loading CNA data.\n", file = stdout())
	 if(file.exists(arguments$cna)){
	   CNA_dt <- fread(arguments$cna, header = TRUE, sep = "\t")
	   if(!("cell_line" %in% colnames(CNA_dt)) & !("-" %in% colnames(CNA_dt))){
		   stop(paste0("Data ", arguments$cna, " does not contain columns of cell_line or disease."))
	   }
	   setkey(CNA_dt, gene_name)
	   if(!is.null(gene_B_filtered)){
	     CNA_dt <- CNA_dt[.(gene_B_filtered)]
		 }
		 setkey(CNA_dt, gene_name)

	 }else{
	 	 stop(paste0("File path ", arguments$cna, " is not correctly specified."))
	 }
	 cat("Read ", paste0(arguments$cna, " with ", nrow(CNA_dt), " rows after filtering.\n"), file = stdout())

	 
	 records <- merged_CSEA_Mediation[, c("ID", "Gene_A", "Gene_B", "cell_line")]

	 cat("Step 4.2: removing CNA confounding.\n", file = stdout())
	 if(arguments$core < 2){
	   for(index in 1:nrow(records)){
		   CNA_confounder_analysis_data_frame[index,] <- CNA_confounder_exp_glm_analysis("ID" = records[index, "ID"], "Gene_A" = records[index, "Gene_A"], "Gene_B" = records[index, "Gene_B"], "Cell_lines" = records[index, "cell_line"][[1]])
	   }
	   removed_index_CNA_confounder <- CNA_confounder_analysis_data_frame$ID[which(CNA_confounder_analysis_data_frame$P_value_CNA_confounding < arguments$cna_adjust_Pvalue_cutoff)]
	 }else{
	   options(cl.cores = arguments$core)
	   this.cluster <- makeCluster(getOption("cl.cores", 2))
	   clusterCall(cl = this.cluster, fun = function(){
	     library(data.table); library(dplyr)
	   })
	   if(formalArgs(clusterExport)[2] %in% "list"){
	     clusterExport(cl = this.cluster, list = c("CNA_confounder_exp_glm_analysis", "exp_tpm_dt", "mutations_dt", "CNA_dt", "records", "arguments"))
	   }else{
	     clusterExport(cl = this.cluster, varlist = c("CNA_confounder_exp_glm_analysis", "exp_tpm_dt", "mutations_dt", "CNA_dt", "records", "arguments"))
	   }
	   
	   CNA_confounder_pval_results <- 
	   parLapply(cl = this.cluster,
	     1:nrow(records),
	     function(idx) {
	       args1 <- as.list(c(records[idx, c("ID", "Gene_A", "Gene_B", "cell_line")]))
	       formals(CNA_confounder_exp_glm_analysis) <- args1
	       rtn <- replicate(1, CNA_confounder_exp_glm_analysis(), simplify = FALSE)
	       return(rtn)
	     }
	   )
	   stopCluster(this.cluster)
	   for(j in 1:length(CNA_confounder_pval_results)){
	     CNA_confounder_analysis_data_frame[j,] <- CNA_confounder_pval_results[[j]][[1]]
	   }
	   removed_index_CNA_confounder <- CNA_confounder_analysis_data_frame$ID[which(CNA_confounder_analysis_data_frame$P_value_CNA_confounding < arguments$cna_adjust_Pvalue_cutoff)]
	   rm(CNA_confounder_pval_results)
	 }
	 
	 merged_CSEA_Mediation <- merge(merged_CSEA_Mediation, CNA_confounder_analysis_data_frame, by = "ID")
	 
	 cat(paste0((nrow(CNA_confounder_analysis_data_frame) - length(removed_index_CNA_confounder)), " of ", nrow(CNA_confounder_analysis_data_frame), " records were kept after removing potential CNA confounding.\n"), file = stdout())
	 rm(CNA_dt)
	 rm(mutations_dt)
	 rm(exp_tpm_dt)
}else{
	cat("Step 4: \"Filtering trans regulation results that can be explained by cis CNA regulation\" is disabled.\n", file = stdout())
}

#merged_CSEA_Mediation <- merged_CSEA_Mediation[!(merged_CSEA_Mediation$ID %in% removed_index_CNA_confounder),]

write.table(merged_CSEA_Mediation, file = paste0("./", arguments$folder, "/", arguments$prefix, "combined_analysis_results_Step4.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Step 4 completed.\n\n", file = stdout())


######### Step 5: Combining results to GACP scores #########

cat("Step 5: combining the results and generating final results.\n", file = stdout())

######### Step 5.1: calculate weight for combining p-values #########
liptak <- function(p, w) {
	 p[which(p == 1)] <- 0.99999
	 
	 if(missing(w)) {
	 	 w <- rep(1, length(p))/length(p)
	 }else{
	 	 if(length(w) != length(p)){
     	 stop("Length of p and w must equal!")
     }
	 }
	 
	 Zi <- -qnorm(p) 
	 Z <- sum(w*Zi) / sqrt(sum(w^2))
	 p.val <- pnorm(-Z)
	 return(p.val)
}

if(!exists("merged_CSEA_Mediation")){
	merged_CSEA_Mediation <- fread(paste0(
	"./", arguments$folder, "/", arguments$prefix, "combined_analysis_results_Step4.tsv"), header = TRUE, sep = "\t")
	if("P_value_CNA_confounding" %in% colnames(merged_CSEA_Mediation)){
		merged_CSEA_Mediation <- merged_CSEA_Mediation[which(merged_CSEA_Mediation$P_value_CNA_confounding >= 0.01),]
		merged_CSEA_Mediation$P_value_CNA_confounding <- NULL
	}
}

cis_line <- which(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B)
trans_line <- which(merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B)
pvalues_cis_df <- merged_CSEA_Mediation[cis_line, 
                    c("P_val", "P_delta", "ID", "Mutation_genome_change", "Gene_A", "Gene_B", "Mutation_protein_change", "Mutation_type")]
pvalues_trans_df <- merged_CSEA_Mediation[trans_line, 
                     c("P_val", "P_gamma_1", "ID", "Mutation_genome_change", "Gene_A", "Gene_B", "Mutation_protein_change", "Mutation_type")]
colnames(pvalues_cis_df) <- c("p_CSEA", "p_cis", "ID", "Mutation_genome_change", "Gene_A", "Gene_B", "Mutation_protein_change", "Mutation_type")
colnames(pvalues_trans_df) <- c("p_CSEA", "p_trans", "ID", "Mutation_genome_change", "Gene_A", "Gene_B", "Mutation_protein_change", "Mutation_type")

pvalues_cis_df$p_CSEA[which(is.na(pvalues_cis_df$p_CSEA))] <- 1
pvalues_cis_df$p_CSEA[which(pvalues_cis_df$p_CSEA < 1e-8)] <- 1e-8
pvalues_cis_df$p_cis[which(is.na(pvalues_cis_df$p_cis))] <- 1
pvalues_cis_df$p_cis[which(pvalues_cis_df$p_cis < 0.00001)] <- 0.00001

pvalues_trans_df$p_CSEA[which(is.na(pvalues_trans_df$p_CSEA))] <- 1
pvalues_trans_df$p_CSEA[which(pvalues_trans_df$p_CSEA < 1e-8)] <- 1e-8
pvalues_trans_df$p_trans[which(is.na(pvalues_trans_df$p_trans))] <- 1
pvalues_trans_df$p_trans[which(pvalues_trans_df$p_trans < 0.00001)] <- 0.00001

combined <- NULL
unique_trans_genes <- unique(pvalues_trans_df$Gene_A)
for(index in 1:length(unique_trans_genes)){
	pvalues_trans_df_new <- pvalues_trans_df[which(pvalues_trans_df$Gene_A==unique_trans_genes[index]),]
	unique_mut_genome_change <- unique(pvalues_trans_df_new$Mutation_genome_change)
	for(idx in 1:length(unique_mut_genome_change)){
		subset <- pvalues_trans_df_new[pvalues_trans_df_new$Mutation_genome_change %in% unique_mut_genome_change[idx]]
		combined <- rbind(combined, data.frame("Gene_A" = unique_trans_genes[index], "Mutation_genome_change" = unique_mut_genome_change[idx], "p_CSEA" = median(subset$p_CSEA), "p_trans" = median(subset$p_trans), "Mutation_protein_change" = subset$Mutation_protein_change[1], "ID" = paste0(subset$ID, collapse = ','), "Gene_B" = paste0(subset$Gene_B, collapse = ','), "Mutation_type" = subset$Mutation_type[1]))
	}
}
rm(pvalues_trans_df_new)

combined_cis <- NULL
if(length(which(duplicated(pvalues_cis_df$Mutation_genome_change))) > 0){
	pvalues_cis_df <- unique(pvalues_cis_df$Gene_A)
	for(index in 1 : length(pvalues_cis_df)){
	  if((index %% (floor(length(pvalues_cis_df)/10))) == 0) cat(paste0("#", index, "/", length(pvalues_cis_df), "\n"), file = stdout())
	  pvalues_cis_df_new <- pvalues_cis_df[which(pvalues_cis_df$Gene_A == pvalues_cis_df[index]),]
	  unique_mut_genome_change <- unique(pvalues_cis_df_new$Mutation_genome_change)
	  for(idx in 1:length(unique_mut_genome_change)){
		 	subset <- pvalues_cis_df_new[pvalues_cis_df_new$Mutation_genome_change %in% unique_mut_genome_change[idx]]
		 	combined_cis <- rbind(combined_cis, data.frame("Gene_A" = pvalues_cis_df[index], "Mutation_genome_change" = unique_mut_genome_change[idx], "p_CSEA" = median(subset$p_CSEA), "p_cis" = median(subset$p_cis), "Mutation_protein_change" = subset$Mutation_protein_change[1], "ID" = paste0(subset$ID, collapse = ','), "Gene_B" = paste0(subset$Gene_B, collapse = ','), "Mutation_type" = subset$Mutation_type[1]))
	   
	  }
	}
 	rm(pvalues_cis_df_new)

}else{
 	combined_cis <- pvalues_cis_df[, c("Gene_A", "Mutation_genome_change", "p_CSEA", "p_cis", "Mutation_protein_change", "ID", "Gene_B", "Mutation_type")]
}


cat(paste0("Calculating weights of mediation (scores from step 2) scores.\n"), file = stdout())
weight_mediations <- c(1/10, 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 1, seq(2, 10, by = 1))


MLFC_metrics <- function(weight_mediation) {
  
  MLFC_dat <- NULL
  # cat(paste0("Calculating ", " ", weight_mediation, "\n"), file = stdout())
	 
  combined_cis$combined_P_vals <- apply(combined_cis[, c("p_CSEA", "p_cis")], 1, liptak, c(1, weight_mediation))
  
  combined_cis$combined_P_vals_adjusted <- p.adjust(combined_cis$combined_P_vals, method = "BH")
  
  combined$combined_P_vals <- apply(combined[, c("p_CSEA", "p_trans")], 1, liptak, c(1, weight_mediation))
  
  combined$combined_P_vals_adjusted <- p.adjust(combined$combined_P_vals, method = "BH")

  combined_cis <- combined_cis[order(combined_cis$combined_P_vals),]
  
  running_sum <- 0
  for(i in 1:nrow(combined_cis)){
		pval_i <- combined_cis$combined_P_vals_adjusted[i]
		quantile_i <- i / nrow(combined_cis)
		tmp <- abs(log2(pval_i / quantile_i))
		if(length(tmp) > 0){
			if(!is.nan(tmp)){
				running_sum <- running_sum + tmp
			}
		}
  }
  if(nrow(combined_cis) < 1){
		MLFC_val_cis <- 0
  }else{
		MLFC_val_cis <- running_sum / nrow(combined_cis)
  }
  MLFC_val_cis <- running_sum / nrow(combined_cis)
  
  # pvalues_trans_df_tmp <- pvalues_trans_df_tmp[pvalues_trans_df_tmp$combined_P_vals_adjusted >= cutoff_FDR,]
  combined <- combined[order(combined$combined_P_vals),]
  
  running_sum <- 0
  if(nrow(combined) > 0){
	  for(i in 1:nrow(combined)){
	    pval_i <- combined$combined_P_vals_adjusted[i]
	    quantile_i <- i / nrow(combined)
	    tmp <- abs(log2(pval_i / quantile_i))
	    if(length(tmp) > 0){
		     if(!is.nan(tmp)){
		       running_sum <- running_sum + tmp
		     }
		  }
	  }
  }
  if(nrow(combined) < 1){
		MLFC_val_trans <- 0
  }else{
		MLFC_val_trans <- running_sum / nrow(combined)
  }
  
  MLFC_dat <- data.frame("weight_mediation" = weight_mediation, "MLFC_cis" = MLFC_val_cis, "MLFC_trans" = MLFC_val_trans)
  return(MLFC_dat)
}

final_MLFC_results <- NULL

for(i in 1:length(weight_mediations)){

	 MLFC_results <- MLFC_metrics(weight_mediation = weight_mediations[i])
	 final_MLFC_results <- rbind(final_MLFC_results, MLFC_results)
}


# Find turning points
# Determine the number and the position of extrema (change points, either peaks or pits) in a regular time series.
# see https://github.com/ontogenerator/cpdetectoR

cis_values <- final_MLFC_results$MLFC_cis
names(cis_values) <- 1:length(cis_values)
sorted_cis_values <- sort(cis_values)

trans_values <- final_MLFC_results$MLFC_trans
names(trans_values) <- 1 : length(trans_values)
sorted_trans_values <- sort(trans_values)
set.seed(arguments$seed)

cat(paste0("Weights used for Liptak's p-value combination are chosen based on MLFC metric turning points.\n"), file = stdout())
cat(paste0("After calculation, the weights are shown below:\n"), file = stdout())
	
if(length(sorted_cis_values) > 0){
	tp_cis <- 0
	if(as.numeric(names(sorted_cis_values)[1]) > 1){
		tp_cis <- as.numeric(names(sorted_cis_values)[1])
	}else{
		time <- cp_wrapper(sorted_cis_values[1:10], FALSE, "KS", 2)$Time
		time <- time[time > 0]
		if(length(time) > 1){
			time <- min(time)
		}
		time <- floor(time)
		tp_cis <- as.numeric(names(sorted_cis_values)[time])
	}

	tp_trans <- 0
	if(as.numeric(names(sorted_trans_values)[1]) > 1){
		tp_trans <- as.numeric(names(sorted_trans_values)[1])
	}else{
		time <- cp_wrapper(sorted_trans_values[1:10], FALSE, "KS", 2)$Time
		time <- time[time > 0]
		if(length(time) > 1){
			time <- min(time)
		}
		time <- floor(time)
		tp_trans <- as.numeric(names(sorted_trans_values)[time])

	}

	cat(paste0("weights [CSEA, mediation] for cis- mutations are [1, ", final_MLFC_results$weight_mediation[tp_cis], "].\n"), file = stdout())
	
	pvalues_cis_df$weight_mediation <- final_MLFC_results$weight_mediation[tp_cis]
	
}else{
	cat(paste0("Weights of cis-mutation calculation failed, assign 1 as weight_mediation value.\n\n"), file = stdout())
	
	pvalues_cis_df$weight_mediation <- 1

}

if(length(sorted_trans_values) > 0){
	tp_trans <- 0
	if(as.numeric(names(sorted_trans_values)[1]) > 1){
	 	tp_trans <- as.numeric(names(sorted_trans_values)[1])
	}else{
		time <- cp_wrapper(sorted_trans_values[1:10], FALSE, "KS", 2)$Time
		time <- time[time > 0]
		if(length(time) > 1){
		  time <- min(time)
		}
		time <- floor(time)
		tp_trans <- as.numeric(names(sorted_trans_values)[time])
	}

	cat(paste0("weights [CSEA, mediation] for trans- mutations are [1, ", final_MLFC_results$weight_mediation[tp_trans], "].\n"), file = stdout())

	pvalues_trans_df$weight_mediation <- final_MLFC_results$weight_mediation[tp_trans]
}else{
	cat(paste0("Weights of trans-mutation calculation failed, assign 1 as weight_mediation value.\n\n"), file = stdout())
	pvalues_trans_df$weight_mediation <- 1
}

MLFC_dat_long <- melt(setDT(final_MLFC_results[, c("weight_mediation", "MLFC_cis", "MLFC_trans")]), 
 id.vars = c("weight_mediation"), measure.vars = c("MLFC_cis", "MLFC_trans"), value.name = "MLFC")
MLFC_dat_long$group <- apply(MLFC_dat_long[, c("variable"), drop = FALSE], MARGIN = 1, 
FUN = function(i) paste(i, collapse = "_"))

p <- ggline(MLFC_dat_long, "weight_mediation", "MLFC", color = "variable", group = "variable", 
shape = "variable", size = 0.2, palette = "npg", title = "Mediation score weight estimation for p-value combination") + 
scale_x_discrete(breaks = unique(as.numeric( MLFC_dat_long$weight_mediation)), labels = c("1/10", "1/9", "1/8", "1/7", "1/6", "1/5", "1/4", "1/3", "1/2", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + theme_bw() + 
theme(axis.ticks.y = element_line(size = 0.5), axis.ticks.x = element_line(size=0.5), axis.text.x = element_text(angle = 30, hjust = 0.8, colour = "black"), axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), legend.position = c(0.55, 0.75), legend.background = element_blank(), legend.title = element_blank(), legend.key = element_blank(), legend.key.size = unit(0.35, 'cm'))

ggsave(paste0("./", arguments$folder, "/", arguments$prefix, "weight_mediation_distribution.pdf"), p, width = 5, height = 5)

liptak_2 <- function(p) {
	w <- c(1, p[3])
	p <- p[1:2]
	if(missing(w)) {
	  w <- rep(1, length(p)) / length(p)
	}else{
	  if(length(w) != length(p))
		  stop("Length of p and w must equal!")
	}
	Zi <- -qnorm(p) 
	Z <- sum(w*Zi) / sqrt(sum(w^2))
	p.val <- pnorm(-Z)
	return(p.val)
}

ID_sig_cis <- c()
ID_sig_trans <- c()


######### Step 5.2: combine p-values #########

if(nrow(pvalues_cis_df) > 0){
	pvalues_cis_df$combined_P_vals <- apply(pvalues_cis_df[, c("p_CSEA", "p_cis", "weight_mediation")], 1 , liptak_2)
	pvalues_cis_df$combined_P_vals_adjusted <- p.adjust(pvalues_cis_df$combined_P_vals, method = "BH")

	ID_sig_cis <- pvalues_cis_df$ID[pvalues_cis_df$combined_P_vals_adjusted < arguments$FDR_cutoff_cis] # 
	cat(paste0("Number of final cis-composite mutations is: ", length(ID_sig_cis), "\n"), file = stdout())
}else{
	pvalues_cis_df$combined_P_vals <- 1
	pvalues_cis_df$combined_P_vals_adjusted <- p.adjust(pvalues_cis_df$combined_P_vals, method = "BH")
	cat(paste0("Number of final cis-composite mutations is: ", length(ID_sig_cis), "\n"), file = stdout())
}

if(nrow(pvalues_trans_df) > 0){
  pvalues_trans_df$combined_P_vals <- apply(pvalues_trans_df[, c("p_CSEA", "p_trans", "weight_mediation")], 1, liptak_2)
	pvalues_trans_df$combined_P_vals_adjusted <- p.adjust(pvalues_trans_df$combined_P_vals, method = "BH")

  ID_sig_trans <- pvalues_trans_df$ID[pvalues_trans_df$combined_P_vals_adjusted < arguments$FDR_cutoff_trans]
  cat(paste0("Number of final trans-composite mutations is: ", length(ID_sig_trans), "\n"), file = stdout())
}else{
	pvalues_trans_df$combined_P_vals <- 1
	pvalues_trans_df$combined_P_vals_adjusted <- p.adjust(pvalues_trans_df$combined_P_vals, method = "BH")
	cat(paste0("Number of final trans-composite mutations is: ", length(ID_sig_trans), "\n"), file = stdout())
}

######### Step 5.3: obtain the phenotypes based on scores #########

merged_CSEA_Mediation$Effect_on_cell_phenotype <- rep("", nrow(merged_CSEA_Mediation))

merged_CSEA_Mediation$Effect_on_cell_phenotype[(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$ES * merged_CSEA_Mediation$delta > 0) | (merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$ES * merged_CSEA_Mediation$gamma_1 > 0)] <- "Decrease_dependency"

merged_CSEA_Mediation$Effect_on_cell_phenotype[(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$ES * merged_CSEA_Mediation$delta < 0) | (merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$ES * merged_CSEA_Mediation$gamma_1 < 0)] <- "Increase_dependency"

merged_CSEA_Mediation$Role_in_cancer <- rep("", nrow(merged_CSEA_Mediation))

merged_CSEA_Mediation$Role_in_cancer[(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$ES * merged_CSEA_Mediation$delta > 0) | (merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$ES * merged_CSEA_Mediation$gamma_1 > 0)] <- "Tumor_suppressive"

merged_CSEA_Mediation$Role_in_cancer[(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$ES * merged_CSEA_Mediation$delta < 0) | (merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$ES * merged_CSEA_Mediation$gamma_1 < 0)] <- "Oncogenic"


merged_CSEA_Mediation$GACP_score <- NA
merged_CSEA_Mediation$GACP_score[merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B] <- merged_CSEA_Mediation$ES[merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B] * merged_CSEA_Mediation$delta[merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B]

merged_CSEA_Mediation$GACP_score[merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B] <- merged_CSEA_Mediation$ES[merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B] * merged_CSEA_Mediation$gamma_1[merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B]

merged_CSEA_Mediation_cis_sig <- merged_CSEA_Mediation[merged_CSEA_Mediation$ID %in% ID_sig_cis, c("ID", "index", "Gene_A", "Gene_B", "Mutation_type", "Mutation_genome_change", "Mutation_protein_change", "ES", "NES", "P_val", "P_val_adjusted", "Effect_on_cell_phenotype", "Role_in_cancer", "GACP_score")]
merged_CSEA_Mediation_trans_sig <- merged_CSEA_Mediation[merged_CSEA_Mediation$ID %in% ID_sig_trans, c("ID", "index", "Gene_A", "Gene_B", "Mutation_type", "Mutation_genome_change", "Mutation_protein_change", "ES", "NES", "P_val", "P_val_adjusted", "Effect_on_cell_phenotype", "Role_in_cancer", "GACP_score")]


######### Step 5.4: Generate composite cis- and trans- genetic alterations (GAs) and associated scores #########
cat(paste0("Generate composite cis- and trans- genetic alterations (GAs) and associated scores.\n"), file = stdout())

write.table(merged_CSEA_Mediation_cis_sig, file = paste0("./", arguments$folder, "/", arguments$prefix, "identified_cis_composite_genetic_alteration_and_cell_phenotype_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(merged_CSEA_Mediation_trans_sig, file = paste0("./", arguments$folder, "/", arguments$prefix, "identified_trans_composite_genetic_alteration_and_cell_phenotype_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

if(!exists("mutations_dt")){
	
  mutations_dt <- fread(arguments$ga, header = TRUE, sep = "\t")
  setkey(mutations_dt, gene_name)
  if(!("cell_line" %in% colnames(mutations_dt)) & !("-" %in% colnames(mutations_dt))){
    stop(paste0("Data ", arguments$ga, " does not contain columns of cell_line or disease."))
  }
	 
}

# merged_CSEA_Mediation_cis_sig is split into weighted individual mutations.


######### Step 5.5: Generate statistically significant cis GA reports #########

if(nrow(pvalues_cis_df) > 0){
	cat(paste0("Generate statistically significant cis GA reports.\n"), file = stdout())
	merged2_CSEA_Mediation_cis_sig <- merge(merged_CSEA_Mediation_cis_sig, pvalues_cis_df[, c("p_CSEA", "p_cis", "ID", "combined_P_vals", "combined_P_vals_adjusted")], by= "ID", all.x = TRUE)
	expanded_CSEA_Mediation_cis_sig <- NULL
	unq_gene_A <- unique(merged2_CSEA_Mediation_cis_sig$Gene_A)
	for(idx in 1:length(unq_gene_A)){
		tmp2 <- merged2_CSEA_Mediation_cis_sig[which(merged2_CSEA_Mediation_cis_sig$Gene_A %in% unq_gene_A[idx]),]

		tmp <- tmp2[which(tmp2$Effect_on_cell_phenotype %in% "Increase_dependency"),]
		if(nrow(tmp) > 0){
			# tmp <- merged2_CSEA_Mediation_cis_sig[which(merged2_CSEA_Mediation_cis_sig$Gene_A %in% "EFNA5"),] 
			Gene_A <- unq_gene_A[idx]
			mutpos_split <- strsplit(tmp$Mutation_genome_change, split = ",")[[1]];
			mutations_dat <- mutations_dt[.(Gene_A), ]
			mutations_dat <- mutations_dat[which(mutations_dat$genome_change %in% mutpos_split),] # all mutations not limited to expressed genes
			freq <- table(mutations_dat$genome_change)
			mutations_vector <- names(freq)
			mutations_vector_percentage <- freq/sum(freq)
			freq2 <- table(mutations_dat$protein_change)
			mutations_vector2 <- names(freq2)
			mutations_vector_percentage2 <- freq2/sum(freq2)

			for(mut_genomics in mutations_vector){
			  pval <- unique(tmp$combined_P_vals)
			  pval_adj <- unique(tmp$combined_P_vals_adjusted)
			  if(length(pval) > 1){
			    pval <- median(tmp$combined_P_vals)
			    pval_adj <- median(tmp$combined_P_vals_adjusted)
			  }
			  expanded_CSEA_Mediation_cis_sig_tmp <- data.frame("ID" = paste(tmp$ID, collapse = ','), "index" = paste(unique(tmp$index), collapse = ','), "Gene_A" = Gene_A, "Gene_B" = paste(unique(tmp$Gene_B), collapse = ','), "Mutation_type" = paste(unique(mutations_dat$var_class[which(mutations_dat$genome_change %in% mut_genomics)]), collapse = ','), "Mutation_genome_change" = mut_genomics, "Mutation_genome_change_count" = as.character(freq[mut_genomics]), "Mutation_genome_change_percentage" = round(mutations_vector_percentage[mut_genomics], digits = 3), "Mutation_protein_change" = paste(unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Mutation_protein_change_count" = as.character(sum(freq2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])])), "Mutation_protein_change_percentage" = sum(round(mutations_vector_percentage2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])], digits = 3)), "Cell_lines" = paste(unique(mutations_dat$cell_line[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Effect_on_cell_phenotype" = paste(unique(tmp$Effect_on_cell_phenotype), collapse = ','), "Role_in_cancer"= paste(unique(tmp$Role_in_cancer), collapse = ','), "combined_P_vals"= pval, "combined_P_vals_adjusted" = pval_adj, "GACP_score" = median(tmp$GACP_score))
			  rownames(expanded_CSEA_Mediation_cis_sig_tmp) <- ""
			  expanded_CSEA_Mediation_cis_sig <- rbind(expanded_CSEA_Mediation_cis_sig, expanded_CSEA_Mediation_cis_sig_tmp)
			}
		}

		tmp <- tmp2[which(tmp2$Effect_on_cell_phenotype %in% "Decrease_dependency"),]
		if(nrow(tmp) > 0){
			# tmp <- merged2_CSEA_Mediation_cis_sig[which(merged2_CSEA_Mediation_cis_sig$Gene_A %in% "EFNA5"),] 
			Gene_A <- unq_gene_A[idx]
			mutpos_split <- strsplit(tmp$Mutation_genome_change, split = ",")[[1]];
			mutations_dat <- mutations_dt[.(Gene_A), ]
			mutations_dat <- mutations_dat[which(mutations_dat$genome_change %in% mutpos_split),] # all mutations not limited to expressed genes
			freq <- table(mutations_dat$genome_change)
			mutations_vector <- names(freq)
			mutations_vector_percentage <- freq/sum(freq)
			freq2 <- table(mutations_dat$protein_change)
			mutations_vector2 <- names(freq2)
			mutations_vector_percentage2 <- freq2/sum(freq2)

			for(mut_genomics in mutations_vector){
				pval <- unique(tmp$combined_P_vals)
				pval_adj <- unique(tmp$combined_P_vals_adjusted)
				if(length(pval) > 1){
				  pval <- median(tmp$combined_P_vals)
				  pval_adj <- median(tmp$combined_P_vals_adjusted)
				}
				expanded_CSEA_Mediation_cis_sig_tmp <- data.frame("ID" = paste(tmp$ID, collapse = ','), "index" = paste(unique(tmp$index), collapse = ','), "Gene_A" = Gene_A, "Gene_B" = paste(unique(tmp$Gene_B), collapse = ','), "Mutation_type" = paste(unique(mutations_dat$var_class[which(mutations_dat$genome_change %in% mut_genomics)]), collapse = ','), "Mutation_genome_change" = mut_genomics, "Mutation_genome_change_count" = as.character(freq[mut_genomics]), "Mutation_genome_change_percentage" = round(mutations_vector_percentage[mut_genomics], digits = 3), "Mutation_protein_change" = paste(unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Mutation_protein_change_count" = as.character(sum(freq2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])])), "Mutation_protein_change_percentage" = sum(round(mutations_vector_percentage2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])], digits = 3)), "Cell_lines" = paste(unique(mutations_dat$cell_line[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Effect_on_cell_phenotype" = paste(unique(tmp$Effect_on_cell_phenotype), collapse = ','), "Role_in_cancer"= paste(unique(tmp$Role_in_cancer), collapse = ','), "combined_P_vals"= pval, "combined_P_vals_adjusted" = pval_adj, "GACP_score" = median(tmp$GACP_score))
				rownames(expanded_CSEA_Mediation_cis_sig_tmp) <- ""
				expanded_CSEA_Mediation_cis_sig <- rbind(expanded_CSEA_Mediation_cis_sig, expanded_CSEA_Mediation_cis_sig_tmp)
		 	}
		}
	 
	}
	if(length(expanded_CSEA_Mediation_cis_sig) > 0){
		expanded_CSEA_Mediation_cis_sig <- expanded_CSEA_Mediation_cis_sig %>% group_by(Gene_A, Mutation_genome_change) %>% slice_min(combined_P_vals_adjusted) %>% distinct()
		write.table(expanded_CSEA_Mediation_cis_sig, file = paste0("./", arguments$folder,"/", arguments$prefix, "identified_cis_genetic_alteration_and_cell_phenotype_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
	rm(expanded_CSEA_Mediation_cis_sig_tmp)
}


# for trans sig

######### Step 5.5: Generate statistically significant trans GA reports #########
if(nrow(pvalues_trans_df) > 0){
	cat(paste0("Generate statistically significant trans GA reports.\n"), file = stdout())

	merged2_CSEA_Mediation_trans_sig <- merge(merged_CSEA_Mediation_trans_sig, pvalues_trans_df[, c("p_CSEA", "p_trans", "ID", "combined_P_vals", "combined_P_vals_adjusted")], by= "ID", all.x = TRUE)
	expanded_CSEA_Mediation_trans_sig <- NULL
	unq_gene_A <- unique(merged2_CSEA_Mediation_trans_sig$Gene_A)


	for(idx in 1:length(unq_gene_A)){

		tmp2 <- merged2_CSEA_Mediation_trans_sig[which(merged2_CSEA_Mediation_trans_sig$Gene_A %in% unq_gene_A[idx]),]

		tmp <- tmp2[which(tmp2$Effect_on_cell_phenotype %in% "Increase_dependency"),]
		if(nrow(tmp) > 0){
			# tmp <- merged2_CSEA_Mediation_trans_sig[which(merged2_CSEA_Mediation_trans_sig$Gene_A %in% "EFNA5"),] 
			Gene_A <- unq_gene_A[idx]
			mutpos_split <- strsplit(tmp$Mutation_genome_change, split = ",")[[1]];
			mutations_dat <- mutations_dt[.(Gene_A), ]
			mutations_dat <- mutations_dat[which(mutations_dat$genome_change %in% mutpos_split),] # all mutations not limited to expressed genes
			freq <- table(mutations_dat$genome_change)
			mutations_vector <- names(freq)
			mutations_vector_percentage <- freq/sum(freq)
			freq2 <- table(mutations_dat$protein_change)
			mutations_vector2 <- names(freq2)
			mutations_vector_percentage2 <- freq2/sum(freq2)

			for(mut_genomics in mutations_vector){
			  pval <- unique(tmp$combined_P_vals)
			  pval_adj <- unique(tmp$combined_P_vals_adjusted)
			  if(length(pval) > 1){
			 		pval <- median(tmp$combined_P_vals)
			 		pval_adj <- median(tmp$combined_P_vals_adjusted)
			  }
			  expanded_CSEA_Mediation_trans_sig_tmp <- data.frame("ID" = paste(tmp$ID, collapse = ','), "index" = paste(unique(tmp$index), collapse = ','), "Gene_A" = Gene_A, "Gene_B" = paste(unique(tmp$Gene_B), collapse = ','), "Mutation_type" = paste(unique(mutations_dat$var_class[which(mutations_dat$genome_change %in% mut_genomics)]), collapse = ','), "Mutation_genome_change" = mut_genomics, "Mutation_genome_change_count" = as.character(freq[mut_genomics]), "Mutation_genome_change_percentage" = round(mutations_vector_percentage[mut_genomics], digits = 3), "Mutation_protein_change" = paste(unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Mutation_protein_change_count" = as.character(sum(freq2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])])), "Mutation_protein_change_percentage" = sum(round(mutations_vector_percentage2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])], digits = 3)), "Cell_lines" = paste(unique(mutations_dat$cell_line[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Effect_on_cell_phenotype" = paste(unique(tmp$Effect_on_cell_phenotype), collapse = ','), "Role_in_cancer"= paste(unique(tmp$Role_in_cancer), collapse = ','), "combined_P_vals"= pval, "combined_P_vals_adjusted" = pval_adj, "GACP_score" = median(tmp$GACP_score))
			  rownames(expanded_CSEA_Mediation_trans_sig_tmp) <- ""
			  expanded_CSEA_Mediation_trans_sig <- rbind(expanded_CSEA_Mediation_trans_sig, expanded_CSEA_Mediation_trans_sig_tmp)
			}
		}

		tmp <- tmp2[which(tmp2$Effect_on_cell_phenotype %in% "Decrease_dependency"),]
		if(nrow(tmp) > 0){
			# tmp <- merged2_CSEA_Mediation_trans_sig[which(merged2_CSEA_Mediation_trans_sig$Gene_A %in% "EFNA5"),] 
			Gene_A <- unq_gene_A[idx]
			mutpos_split <- strsplit(tmp$Mutation_genome_change, split = ",")[[1]];
			mutations_dat <- mutations_dt[.(Gene_A), ]
			mutations_dat <- mutations_dat[which(mutations_dat$genome_change %in% mutpos_split),] # all mutations not limited to expressed genes
			freq <- table(mutations_dat$genome_change)
			mutations_vector <- names(freq)
			mutations_vector_percentage <- freq/sum(freq)
			freq2 <- table(mutations_dat$protein_change)
			mutations_vector2 <- names(freq2)
			mutations_vector_percentage2 <- freq2/sum(freq2)

			for(mut_genomics in mutations_vector){
				pval <- unique(tmp$combined_P_vals)
				pval_adj <- unique(tmp$combined_P_vals_adjusted)
				if(length(pval) > 1){
					pval <- median(tmp$combined_P_vals)
					pval_adj <- median(tmp$combined_P_vals_adjusted)
				}
				expanded_CSEA_Mediation_trans_sig_tmp <- data.frame("ID" = paste(tmp$ID, collapse = ','), "index" = paste(unique(tmp$index), collapse = ','), "Gene_A" = Gene_A, "Gene_B" = paste(unique(tmp$Gene_B), collapse = ','), "Mutation_type" = paste(unique(mutations_dat$var_class[which(mutations_dat$genome_change %in% mut_genomics)]), collapse = ','), "Mutation_genome_change" = mut_genomics, "Mutation_genome_change_count" = as.character(freq[mut_genomics]), "Mutation_genome_change_percentage" = round(mutations_vector_percentage[mut_genomics], digits = 3), "Mutation_protein_change" = paste(unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Mutation_protein_change_count" = as.character(sum(freq2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])])), "Mutation_protein_change_percentage" = sum(round(mutations_vector_percentage2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])], digits = 3)), "Cell_lines" = paste(unique(mutations_dat$cell_line[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Effect_on_cell_phenotype" = paste(unique(tmp$Effect_on_cell_phenotype), collapse = ','), "Role_in_cancer"= paste(unique(tmp$Role_in_cancer), collapse = ','), "combined_P_vals"= pval, "combined_P_vals_adjusted" = pval_adj, "GACP_score" = median(tmp$GACP_score))
				rownames(expanded_CSEA_Mediation_trans_sig_tmp) <- ""
				expanded_CSEA_Mediation_trans_sig <- rbind(expanded_CSEA_Mediation_trans_sig, expanded_CSEA_Mediation_trans_sig_tmp)
			}
		}
	 
	}
	if(length(expanded_CSEA_Mediation_trans_sig) > 0){
		expanded_CSEA_Mediation_trans_sig <- expanded_CSEA_Mediation_trans_sig %>% group_by(Gene_A, Mutation_genome_change) %>% slice_min(combined_P_vals_adjusted) %>% distinct()

		write.table(expanded_CSEA_Mediation_trans_sig, file = paste0("./", arguments$folder,"/", arguments$prefix, "identified_trans_genetic_alteration_and_cell_phenotype_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
	rm(expanded_CSEA_Mediation_trans_sig_tmp)
}


if(arguments$full_report){
	######### Step 5.6: Generate cis GA full reports #########
	cat(paste0("Generate cis GA full reports.\n"), file = stdout())

	merged_CSEA_Mediation_cis <- merged_CSEA_Mediation[which(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B),]
	merged2_CSEA_Mediation_cis <- merge(merged_CSEA_Mediation_cis, pvalues_cis_df[, c("p_CSEA", "p_cis", "ID", "combined_P_vals", "combined_P_vals_adjusted")], by = "ID", all.x = TRUE)
	expanded_CSEA_Mediation_cis <- NULL
	unq_gene_A <- unique(merged2_CSEA_Mediation_cis$Gene_A)
	for(idx in 1:length(unq_gene_A)){
		tmp2 <- merged2_CSEA_Mediation_cis[which(merged2_CSEA_Mediation_cis$Gene_A %in% unq_gene_A[idx]),]

		tmp <- tmp2[which(tmp2$Effect_on_cell_phenotype %in% "Increase_dependency"),]
		if(nrow(tmp) > 0){
			# tmp <- merged2_CSEA_Mediation_cis[which(merged2_CSEA_Mediation_cis$Gene_A %in% "EFNA5"),] 
			Gene_A <- unq_gene_A[idx]
			mutpos_split <- strsplit(tmp$Mutation_genome_change, split = ",")[[1]];
			mutations_dat <- mutations_dt[.(Gene_A), ]
			mutations_dat <- mutations_dat[which(mutations_dat$genome_change %in% mutpos_split),] # all mutations not limited to expressed genes
			freq <- table(mutations_dat$genome_change)
			mutations_vector <- names(freq)
			mutations_vector_percentage <- freq/sum(freq)
			freq2 <- table(mutations_dat$protein_change)
			mutations_vector2 <- names(freq2)
			mutations_vector_percentage2 <- freq2/sum(freq2)

			for(mut_genomics in mutations_vector){
			  pval <- unique(tmp$combined_P_vals)
			  pval_adj <- unique(tmp$combined_P_vals_adjusted)
			  if(length(pval) > 1){
			 		pval <- median(tmp$combined_P_vals)
			 		pval_adj <- median(tmp$combined_P_vals_adjusted)
			  }
			  expanded_CSEA_Mediation_cis_tmp <- data.frame("ID" = paste(tmp$ID, collapse = ','), "index" = paste(unique(tmp$index), collapse = ','), "Gene_A" = Gene_A, "Gene_B" = paste(unique(tmp$Gene_B), collapse = ','), "Mutation_type" = paste(unique(mutations_dat$var_class[which(mutations_dat$genome_change %in% mut_genomics)]), collapse = ','), "Mutation_genome_change" = mut_genomics, "Mutation_genome_change_count" = as.character(freq[mut_genomics]), "Mutation_genome_change_percentage" = round(mutations_vector_percentage[mut_genomics], digits = 3), "Mutation_protein_change" = paste(unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Mutation_protein_change_count" = as.character(sum(freq2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])])), "Mutation_protein_change_percentage" = sum(round(mutations_vector_percentage2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])], digits = 3)), "Cell_lines" = paste(unique(mutations_dat$cell_line[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Effect_on_cell_phenotype" = paste(unique(tmp$Effect_on_cell_phenotype), collapse = ','), "Role_in_cancer"= paste(unique(tmp$Role_in_cancer), collapse = ','), "combined_P_vals"= pval, "combined_P_vals_adjusted" = pval_adj, "GACP_score" = median(tmp$GACP_score))
			  rownames(expanded_CSEA_Mediation_cis_tmp) <- ""
			  expanded_CSEA_Mediation_cis <- rbind(expanded_CSEA_Mediation_cis, expanded_CSEA_Mediation_cis_tmp)
			}
		}

		tmp <- tmp2[which(tmp2$Effect_on_cell_phenotype %in% "Decrease_dependency"),]
		if(nrow(tmp) > 0){
			# tmp <- merged2_CSEA_Mediation_cis[which(merged2_CSEA_Mediation_cis$Gene_A %in% "EFNA5"),] 
			Gene_A <- unq_gene_A[idx]
			mutpos_split <- strsplit(tmp$Mutation_genome_change, split = ",")[[1]];
			mutations_dat <- mutations_dt[.(Gene_A), ]
			mutations_dat <- mutations_dat[which(mutations_dat$genome_change %in% mutpos_split),] # all mutations not limited to expressed genes
			freq <- table(mutations_dat$genome_change)
			mutations_vector <- names(freq)
			mutations_vector_percentage <- freq/sum(freq)
			freq2 <- table(mutations_dat$protein_change)
			mutations_vector2 <- names(freq2)
			mutations_vector_percentage2 <- freq2/sum(freq2)

			for(mut_genomics in mutations_vector){
			pval <- unique(tmp$combined_P_vals)
			pval_adj <- unique(tmp$combined_P_vals_adjusted)
			if(length(pval) > 1){
				pval <- median(tmp$combined_P_vals)
				pval_adj <- median(tmp$combined_P_vals_adjusted)
			}
			expanded_CSEA_Mediation_cis_tmp <- data.frame("ID" = paste(tmp$ID, collapse = ','), "index" = paste(unique(tmp$index), collapse = ','), "Gene_A" = Gene_A, "Gene_B" = paste(unique(tmp$Gene_B), collapse = ','), "Mutation_type" = paste(unique(mutations_dat$var_class[which(mutations_dat$genome_change %in% mut_genomics)]), collapse = ','), "Mutation_genome_change" = mut_genomics, "Mutation_genome_change_count" = as.character(freq[mut_genomics]), "Mutation_genome_change_percentage" = round(mutations_vector_percentage[mut_genomics], digits = 3), "Mutation_protein_change" = paste(unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Mutation_protein_change_count" = as.character(sum(freq2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])])), "Mutation_protein_change_percentage" = sum(round(mutations_vector_percentage2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])], digits = 3)), "Cell_lines" = paste(unique(mutations_dat$cell_line[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Effect_on_cell_phenotype" = paste(unique(tmp$Effect_on_cell_phenotype), collapse = ','), "Role_in_cancer"= paste(unique(tmp$Role_in_cancer), collapse = ','), "combined_P_vals"= pval, "combined_P_vals_adjusted" = pval_adj, "GACP_score" = median(tmp$GACP_score))
			rownames(expanded_CSEA_Mediation_cis_tmp) <- ""
			expanded_CSEA_Mediation_cis <- rbind(expanded_CSEA_Mediation_cis, expanded_CSEA_Mediation_cis_tmp)
			}
		}
	 
	}

	rm(expanded_CSEA_Mediation_cis_tmp)
	if(length(expanded_CSEA_Mediation_cis) > 0){
	 	expanded_CSEA_Mediation_cis <- expanded_CSEA_Mediation_cis %>% group_by(Gene_A, Mutation_genome_change) %>% slice_min(combined_P_vals_adjusted) %>% distinct()
	 	write.table(expanded_CSEA_Mediation_cis, file = paste0("./", arguments$folder,"/", arguments$prefix, "cis_genetic_alteration_and_cell_phenotype_full_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
	
}


if(arguments$full_report){
	######### Step 5.7: Generate tran GA full reports #########
	cat(paste0("Generate tran GA full reports.\n"), file = stdout())
	write.table("ID	index	Gene_A	Gene_B	Mutation_type	Mutation_genome_change	Mutation_genome_change_count	Mutation_genome_change_percentage	Mutation_protein_change	Mutation_protein_change_count	Mutation_protein_change_percentage	Cell_lines	Effect_on_cell_phenotype	Role_in_cancer	combined_P_vals	combined_P_vals_adjusted	GACP_score", file = paste0("./", arguments$folder,"/", arguments$prefix, "trans_genetic_alteration_and_cell_phenotype_full_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	merged_CSEA_Mediation_trans <- merged_CSEA_Mediation[which(merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B),]
	merged2_CSEA_Mediation_trans <- merge(merged_CSEA_Mediation_trans, pvalues_trans_df[, c("p_CSEA", "p_trans", "ID", "combined_P_vals", "combined_P_vals_adjusted")], by = "ID", all.x = TRUE)
	expanded_CSEA_Mediation_trans <- NULL
	unq_gene_A <- unique(merged2_CSEA_Mediation_trans$Gene_A)
	for(idx in 1:length(unq_gene_A)){

		tmp2 <- merged2_CSEA_Mediation_trans[which(merged2_CSEA_Mediation_trans$Gene_A %in% unq_gene_A[idx]),]

		tmp <- tmp2[which(tmp2$Effect_on_cell_phenotype %in% "Increase_dependency"),]
		if(nrow(tmp) > 0){
			# tmp <- merged2_CSEA_Mediation_trans[which(merged2_CSEA_Mediation_trans$Gene_A %in% "EFNA5"),] 
			Gene_A <- unq_gene_A[idx]
			mutpos_split <- strsplit(tmp$Mutation_genome_change, split = ",")[[1]];
			mutations_dat <- mutations_dt[.(Gene_A), ]
			mutations_dat <- mutations_dat[which(mutations_dat$genome_change %in% mutpos_split),] # all mutations not limited to expressed genes
			freq <- table(mutations_dat$genome_change)
			mutations_vector <- names(freq)
			mutations_vector_percentage <- freq/sum(freq)
			freq2 <- table(mutations_dat$protein_change)
			mutations_vector2 <- names(freq2)
			mutations_vector_percentage2 <- freq2/sum(freq2)

			for(mut_genomics in mutations_vector){
				pval <- unique(tmp$combined_P_vals)
				pval_adj <- unique(tmp$combined_P_vals_adjusted)
				if(length(pval) > 1){
					pval <- median(tmp$combined_P_vals)
					pval_adj <- median(tmp$combined_P_vals_adjusted)
				}
				expanded_CSEA_Mediation_trans_tmp <- data.frame("ID" = paste(tmp$ID, collapse = ','), "index" = paste(unique(tmp$index), collapse = ','), "Gene_A" = Gene_A, "Gene_B" = paste(unique(tmp$Gene_B), collapse = ','), "Mutation_type" = paste(unique(mutations_dat$var_class[which(mutations_dat$genome_change %in% mut_genomics)]), collapse = ','), "Mutation_genome_change" = mut_genomics, "Mutation_genome_change_count" = as.character(freq[mut_genomics]), "Mutation_genome_change_percentage" = round(mutations_vector_percentage[mut_genomics], digits = 3), "Mutation_protein_change" = paste(unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Mutation_protein_change_count" = as.character(sum(freq2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])])), "Mutation_protein_change_percentage" = sum(round(mutations_vector_percentage2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])], digits = 3)), "Cell_lines" = paste(unique(mutations_dat$cell_line[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Effect_on_cell_phenotype" = paste(unique(tmp$Effect_on_cell_phenotype), collapse = ','), "Role_in_cancer"= paste(unique(tmp$Role_in_cancer), collapse = ','), "combined_P_vals"= pval, "combined_P_vals_adjusted" = pval_adj, "GACP_score" = median(tmp$GACP_score))
				rownames(expanded_CSEA_Mediation_trans_tmp) <- ""
				expanded_CSEA_Mediation_trans <- rbind(expanded_CSEA_Mediation_trans, expanded_CSEA_Mediation_trans_tmp)
			}
		}

		tmp <- tmp2[which(tmp2$Effect_on_cell_phenotype %in% "Decrease_dependency"),]
		if(nrow(tmp) > 0){
			# tmp <- merged2_CSEA_Mediation_trans[which(merged2_CSEA_Mediation_trans$Gene_A %in% "EFNA5"),] 
			Gene_A <- unq_gene_A[idx]
			mutpos_split <- strsplit(tmp$Mutation_genome_change, split = ",")[[1]];
			mutations_dat <- mutations_dt[.(Gene_A), ]
			mutations_dat <- mutations_dat[which(mutations_dat$genome_change %in% mutpos_split),] # all mutations not limited to expressed genes
			freq <- table(mutations_dat$genome_change)
			mutations_vector <- names(freq)
			mutations_vector_percentage <- freq/sum(freq)
			freq2 <- table(mutations_dat$protein_change)
			mutations_vector2 <- names(freq2)
			mutations_vector_percentage2 <- freq2/sum(freq2)

			for(mut_genomics in mutations_vector){
				pval <- unique(tmp$combined_P_vals)
				pval_adj <- unique(tmp$combined_P_vals_adjusted)
				if(length(pval) > 1){
				  pval <- median(tmp$combined_P_vals)
				  pval_adj <- median(tmp$combined_P_vals_adjusted)
				}
				expanded_CSEA_Mediation_trans_tmp <- data.frame("ID" = paste(tmp$ID, collapse = ','), "index" = paste(unique(tmp$index), collapse = ','), "Gene_A" = Gene_A, "Gene_B" = paste(unique(tmp$Gene_B), collapse = ','), "Mutation_type" = paste(unique(mutations_dat$var_class[which(mutations_dat$genome_change %in% mut_genomics)]), collapse = ','), "Mutation_genome_change" = mut_genomics, "Mutation_genome_change_count" = as.character(freq[mut_genomics]), "Mutation_genome_change_percentage" = round(mutations_vector_percentage[mut_genomics], digits = 3), "Mutation_protein_change" = paste(unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Mutation_protein_change_count" = as.character(sum(freq2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])])), "Mutation_protein_change_percentage" = sum(round(mutations_vector_percentage2[unique(mutations_dat$protein_change[mutations_dat$genome_change == mut_genomics])], digits = 3)), "Cell_lines" = paste(unique(mutations_dat$cell_line[mutations_dat$genome_change == mut_genomics]), collapse = ','), "Effect_on_cell_phenotype" = paste(unique(tmp$Effect_on_cell_phenotype), collapse = ','), "Role_in_cancer"= paste(unique(tmp$Role_in_cancer), collapse = ','), "combined_P_vals"= pval, "combined_P_vals_adjusted" = pval_adj, "GACP_score" = median(tmp$GACP_score))
				rownames(expanded_CSEA_Mediation_trans_tmp) <- ""
				expanded_CSEA_Mediation_trans <- rbind(expanded_CSEA_Mediation_trans, expanded_CSEA_Mediation_trans_tmp) 
			}
		}
		if(nrow(expanded_CSEA_Mediation_trans) > 999){
			write.table(expanded_CSEA_Mediation_trans, file = paste0("./", arguments$folder, "/", arguments$prefix, "trans_genetic_alteration_and_cell_phenotype_full_results.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
			expanded_CSEA_Mediation_trans <- NULL
		}
	}

	rm(expanded_CSEA_Mediation_trans_tmp)

	if(length(expanded_CSEA_Mediation_trans) > 0){
		write.table(expanded_CSEA_Mediation_trans, file = paste0("./", arguments$folder, "/", arguments$prefix, "trans_genetic_alteration_and_cell_phenotype_full_results.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
		expanded_CSEA_Mediation_trans <- fread(paste0("./", arguments$folder, "/", arguments$prefix, "trans_genetic_alteration_and_cell_phenotype_full_results.txt"),
	              header = TRUE, sep = "\t")
		expanded_CSEA_Mediation_trans <- expanded_CSEA_Mediation_trans %>% group_by(Gene_A, Mutation_genome_change) %>% slice_min(combined_P_vals_adjusted) %>% distinct()

		write.table(expanded_CSEA_Mediation_trans, file = paste0("./", arguments$folder, "/", arguments$prefix, "trans_genetic_alteration_and_cell_phenotype_full_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)
	}
}


# split composite mutations into single mutations per line and computed the relative percentage of single muations in composite mutations. Output can be submitted to CRAVAT server.
if(arguments$BED_report){
	######### Step 5.8: Generate a cis missense BED file (hg19) as the input for CRAVAT server #########
	cat(paste0("Generate a cis missense BED file (hg19) as the input for CRAVAT server.\n"), file = stdout())
	mutation_cis_list <- expanded_CSEA_Mediation_cis_sig[expanded_CSEA_Mediation_cis_sig$Mutation_type == "Missense", c("ID", "Mutation_genome_change", "Gene_A", "Mutation_genome_change_count", "Gene_B")]
	mutation_cis_list$chr <- ""
	mutation_cis_list$pos <- ""
	mutation_cis_list$chain <- "+"
	mutation_cis_list$ref <- ""
	mutation_cis_list$alt <- ""
	mutation_cis_list$label <- paste(mutation_cis_list$ID, mutation_cis_list$Gene_A, mutation_cis_list$Mutation_genome_change_count, mutation_cis_list$Gene_B, sep = "_")
	for(idx in 1:nrow(mutation_cis_list)){
		tmp <- strsplit(mutation_cis_list$Mutation_genome_change, "[^[:alnum:] ]")[[idx]]
		spl <- strsplit(tmp[3], "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl = TRUE)[[1]]
		mutation_cis_list$chr[idx] <- tmp[2]
		mutation_cis_list$pos[idx] <- spl[1]
		mutation_cis_list$ref[idx] <- spl[2]
		mutation_cis_list$alt[idx] <- tmp[4]
	}
	mutation_cis_list$ID <- 1 : nrow(mutation_cis_list)

	write.table(mutation_cis_list[, c("ID", "chr", "pos", "chain", "ref", "alt", "label")], file = paste0("./", arguments$folder, "/", arguments$prefix, "cis_genetic_alteration_pass_threshold.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}


if(arguments$BED_report){
  ######### Step 5.8: Generate a trans missense BED file (hg19) as the input for CRAVAT server #########
	cat(paste0("Generate a trans missense BED file (hg19) as the input for CRAVAT server.\n"), file = stdout())
	mutation_trans_list <- expanded_CSEA_Mediation_trans_sig[expanded_CSEA_Mediation_trans_sig$Mutation_type == "Missense", c("ID", "Mutation_genome_change", "Gene_A", "Mutation_genome_change_count", "Gene_B")]
	mutation_trans_list$chr <- ""
	mutation_trans_list$pos <- ""
	mutation_trans_list$chain <- "+"
	mutation_trans_list$ref <- ""
	mutation_trans_list$alt <- ""
	mutation_trans_list$label <- paste(mutation_trans_list$ID, mutation_trans_list$Gene_A, mutation_trans_list$Mutation_genome_change_count, mutation_trans_list$Gene_B, sep = "_")
	for(idx in 1:nrow(mutation_trans_list)){
	  tmp <- strsplit(mutation_trans_list$Mutation_genome_change, "[^[:alnum:] ]")[[idx]]
	  spl <- strsplit(tmp[3], "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl = TRUE)[[1]]
	  mutation_trans_list$chr[idx] <- tmp[2]
	  mutation_trans_list$pos[idx] <- spl[1]
	  mutation_trans_list$ref[idx] <- spl[2]
	  mutation_trans_list$alt[idx] <- tmp[4]
	}
	mutation_trans_list$ID <- 1 : nrow(mutation_trans_list)

	write.table(mutation_trans_list[, c("ID", "chr", "pos", "chain", "ref", "alt", "label")], file = paste0("./", arguments$folder, "/", arguments$prefix, "trans_genetic_alteration_pass_threshold.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

cat(paste0("Step 5 completed.\n"), file = stdout())
cat(paste0("Program successfully completed.\n"), file = stdout())