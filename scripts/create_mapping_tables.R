# Load libraries, or install if not present

suppressPackageStartupMessages({
  if(!require("data.table")) {
    install.packages("data.table")
    require("data.table")
  }
  if(!require("dplyr")) {
    install.packages("dplyr")
    require("dplyr")
  }
  if(!require("tidyr")) {
    install.packages("tidyr")
    require("tidyr")
  }
})

setDTthreads(as.integer(snakemake@threads))

dbname <- snakemake@params[[1]]

#read taxonomy
taxDB <- fread(snakemake@input[[1]],
               header = FALSE,
               sep = "\t",
               colClasses = "character",
               col.names = c("OTU", "tax"))

#read mappings
mappings<-data.frame()
files<-list.files(path = paste0(snakemake@input[[2]],"/mapped/"), pattern = ".idmapped.txt")

for (file in files){
  mapping <- fread(
    paste0(snakemake@input[[2]],"/mapped/",file),
    header = FALSE,
    sep = "\t",
    col.names = c("readID", "read_time", "add_info")
  ) %>% separate(add_info, c("barcode", "SAMflag", "OTU", "Qlen", "alnlen", "MapID", "NMtag", "alnscore", "MinimapID"), " ") %>% mutate(SeqID=sub(".idmapped.txt","",file))
  mappings<-rbind(mappings,mapping)
}

# Subset mappings based on alignment length, then write them out
mappings <- mappings %>% mutate(Qr = as.numeric(Qlen)/as.numeric(alnlen)) 
mappings_s <- mappings %>% subset(., Qr < 1.25 & Qr > 0.75)
fwrite(mappings_s, snakemake@output[[4]], quote = F, sep = "\t", row.names = F, col.names = T)
fwrite(mappings, snakemake@output[[3]], quote = F, sep = "\t", row.names = F, col.names = T)

# Load in previous mappings in last folder, if there are any
dirs <- list.dirs(path = ".", recursive = FALSE) %>% grep("*_processed", ., value = TRUE)
last_dir <- dirs[length(dirs)-1]
last_map_dir <- paste0(last_dir, "/mapped/temp/all_mappings_filt25.txt")
print(paste0("Found last used directory to be ", last_dir, " - looking for previous mappings here (mapped/temp/all_mappings_filt25.txt)."))
if (file.exists(last_map_dir)) {
  last_map_file <- fread(last_map_dir, header = TRUE, sep = "\t")
  
  # Bind mappings files together and make sure that each row is unique (readIDs should be unique)
  all_mappings <- rbind(last_map_file, mappings_s) %>% .[!duplicated(.), ]
  print("Loaded previous mappings and outputting combined OTUtable + mappings file")
} else {
  print("Could not find previous mappings. Do you have any previous mappings?")
  print("Resuming without previous mappings...")
  all_mappings <- mappings_s %>% .[!duplicated(.), ]
}


# Write out mapping file with all mappings
fwrite(all_mappings, snakemake@output[[2]], quote = F, sep = "\t", row.names = F, col.names = T)

# Convert mappings to OTU counts
mappings_to_otu <- all_mappings %>% 
  select(c("SeqID", "OTU"))

# Define function for 
#join taxonomy with mapping
joined <- taxDB[mappings_to_otu, on = "OTU"]

#transform into "OTU table", where rows are OTU's, columns are sample AND taxonomy
BIOMotutable <- dcast(joined, OTU + tax ~ SeqID, fun.aggregate = length) %>% setDT()
BIOMotutable[,taxonomy := gsub("[a-z]__", "", tax)]
BIOMotutable[,tax := NULL]
BIOMotutable[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") := tstrsplit(taxonomy, ";", fixed=FALSE)]
BIOMotutable[,taxonomy := NULL]
#write out
fwrite(BIOMotutable, snakemake@output[[1]], sep = "\t", col.names = TRUE, na = "NA", quote = FALSE)