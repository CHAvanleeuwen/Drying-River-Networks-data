library(metabaR)
library(stringr)

## this function keep only the rank containning in the rank_vector and format the taxonomic lineage
filter_on_rank <- function(taxo_path, rank_vector){
  rank_filter <- lapply(taxo_path, str_detect, rank_vector);
  filtered_matrix <- matrix(unlist(rank_filter), 
                            nrow=length(unlist(rank_filter[1])))
  filtered_res <- as.vector(taxo_path[apply(filtered_matrix, 2, any)])
  available_rank <- apply(filtered_matrix, 1, any)
  new_path <- lapply(seq_along(available_rank),
                     function(ind){
                       ranknames <- gsub(".*@", "@", filtered_res)
                       if (rank_vector[ind] %in% ranknames) {
                         filtered_res[which(rank_vector[ind]==ranknames)]}
                       else {
                         indn <- 1
                         for(index in ind:1){
                           if(rank_vector[index] %in% ranknames){
                             indn <- which(rank_vector[index]==ranknames)
                             break
                           }
                         }
                         # paste the last taxon name with the active rank and undefined
                         paste(gsub("@.*", "", filtered_res[indn]), 
                               paste("(", 
                                     paste(gsub("@", "", rank_vector[ind]),
                                           "undefined", 
                                           sep=" "),
                                     ")", 
                                     sep=""), 
                               sep=" ")
                       }
                     })
  new_path
}

## format all ecotag path
format_ecotag_path <- function(path, ranks=rank){
  rank_vector <- paste("@", ranks, sep="")
  tmp_data <- strsplit(path, ":")
  filtered_data <- lapply(tmp_data, filter_on_rank, rank_vector)
  
  formated_data <- lapply(filtered_data, 
                          function(taxo_value){
                            gsub("@.*", "", as.character(taxo_value))})
  out <- data.frame(t(matrix(unlist(formated_data), 
                             nrow=length(unlist(formated_data[1])))), stringsAsFactors = F)
  names(out) <- ranks
  out
}

lib_dir <- "~/travail/document/DRYVER/bact02"

mtbl <- readRDS(file.path(lib_dir, "RDS_file", "bact02_clean_clustered_8.rds"))

fasta_generator(mtbl, output_file = file.path(lib_dir, "fasta_to_silva_bact02_dryver.fasta"))

## add silva annotation

mtbl <- silva_annotator(mtbl, 
                        file.path(lib_dir, "results_silva/ssu/exports/bact02_dryver---ssu---otus.csv"), 
                        file.path(lib_dir, "results_silva/ssu/stats/sequence_cluster_map/data/bact02_dryver---ssu---sequence_cluster_map---fasta_to_silva_bact02_dryver.clstr"), 
                        "~/travail/dev/project/R_project_pipeline/metabaR_external_data/tax_slv_ssu_138.1.txt")

mtbl$motus$scientific_name_Silva <- unlist(lapply(as.character(mtbl$motus$lineage_silva), function(x) { 
  v <- unlist(strsplit(x, ";"))
  taxname <- tail(v, n=1)
  if(taxname == "uncultured"){
    taxname <- paste(tail(v, n=2), collapse = "_")
  }
  return(taxname)
} ))

silva_mtd <- read.csv2("~/travail/dev/project/R_project_pipeline/metabaR_external_data/tax_slv_ssu_138.1.txt", 
                       header = F, sep = "\t", 
                       row.names = 1, stringsAsFactors = FALSE)

mtbl$motus$rank_Silva <- silva_mtd[mtbl$motus$lineage_silva, "V3"]

saveRDS(mtbl, file.path(lib_dir, "RDS_file", "bact02_clean_clustered_8_silva.rds"))

clust_mtbl <- readRDS(file.path(lib_dir, "RDS_file", "bact02_clean_motus_9.rds"))
clust_mtbl$motus <- mtbl$motus[rownames(clust_mtbl$motus), ]
saveRDS(clust_mtbl, file.path(lib_dir, "RDS_file", "bact02_clean_motus_9_silva.rds"))
rm(clust_mtbl)


samples <- subset_metabarlist(mtbl, "pcrs", mtbl$pcrs$type=="sample")
samples <- aggregate_pcrs(samples, FUN = FUN_agg_pcrs_sum)


paths <- format_ecotag_path(samples$motus$path,
                            c("superkingdom", "kingdom", "phylum", "class", "order",
                              "family", "genus", "species"))

colnames(paths) <- paste0(colnames(paths), "_ecotag")
samples$motus[, colnames(paths)] <- paths
samples$motus[,"init_count"] <- samples$motus[, "count"]
samples$motus[,"count"] <- colSums(samples$reads)
samples$motus[!is.na(samples$motus$species_name), "species_ecotag"] <- samples$motus[!is.na(samples$motus$species_name), "species_name"]

write.csv2(samples$motus[, !(colnames(samples$motus) %in% c("definition", "x", 
                                                            "init_count", "Kolmo.contaminant", 
                                                            "Kolmo.contaminant.p.value") | 
                               grepl("obiclean", colnames(samples$motus)) |
                               grepl("control_abundance", colnames(samples$motus)) |
                               grepl("sequencing_control", colnames(samples$motus)))], 
           file = file.path(lib_dir, "csv_files", "dryver-bact02-samples_motus_silva.csv"))

clust_mtbl <- readRDS(file.path(lib_dir, "RDS_file", "bact02_clean_motus_9_silva.rds"))
clust_mtbl <- subset_metabarlist(clust_mtbl, "pcrs", clust_mtbl$pcrs$type=="sample")
clust_samples <- aggregate_pcrs(clust_mtbl)

paths <- format_ecotag_path(clust_samples$motus$path,
                            c("superkingdom", "kingdom", "phylum", "class", "order",
                              "family", "genus", "species"))

colnames(paths) <- paste0(colnames(paths), "_ecotag")
clust_samples$motus[, colnames(paths)] <- paths
clust_samples$motus[,"init_count"] <- clust_samples$motus[, "count"]
clust_samples$motus[,"count"] <- colSums(clust_samples$reads)
clust_samples$motus[!is.na(clust_samples$motus$species_name), "species_ecotag"] <- clust_samples$motus[!is.na(clust_samples$motus$species_name), "species_name"]


write.csv2(clust_samples$motus[, !(colnames(clust_samples$motus) %in% c("definition", "x", 
                                                                        "init_count", "Kolmo.contaminant", 
                                                                        "Kolmo.contaminant.p.value") | 
                                     grepl("obiclean", colnames(clust_samples$motus)) |
                                     grepl("control_abundance", colnames(clust_samples$motus)) |
                                     grepl("sequencing_control", colnames(clust_samples$motus)))], 
           file = file.path(lib_dir, "csv_files", "dryver-bact02-clust_samples_motus_silva.csv"))

clust_rds <- list("reads"=t(clust_samples$reads),
                  "motus"=clust_samples$motus[, !(colnames(clust_samples$motus) %in% c("definition", "x", 
                                                                                       "init_count", "Kolmo.contaminant", 
                                                                                       "Kolmo.contaminant.p.value") | 
                                                    grepl("obiclean", colnames(clust_samples$motus)) |
                                                    grepl("control_abundance", colnames(clust_samples$motus)) |
                                                    grepl("sequencing_control", colnames(clust_samples$motus)))],
                  "samples"=clust_samples$samples)

saveRDS(clust_rds, file = file.path(lib_dir, "RDS_file", "dryver-bact02-clust_motus_samples_silva.rds"))

rm(clust_mtbl)
rm(clust_rds)
rm(clust_samples)
rm(paths)
gc()

final_mtbl <- samples
final_mtbl$motus <- final_mtbl$motus[, !(colnames(samples$motus) %in% c("definition", "x") | 
                                           grepl("obiclean", colnames(samples$motus)) |
                                           grepl("control_abundance", colnames(samples$motus)) |
                                           grepl("sequencing_control", colnames(samples$motus)))]
saveRDS(final_mtbl, file = file.path(lib_dir, "RDS_file", "dryver-bact02-result-samples-metabarlist_silva.rds"))

rm(final_mtbl)
