library(downloader)
library(PharmacoGxPrivate)
library(devtools)
library(data.table)

options(stringsAsFactors=FALSE)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]|[(]|[)]"

# cell_all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
# drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))

# curationCell <- cell_all[which(!is.na(cell_all[ , "CTRPv2.cellid"])),]
# curationCell <- curationCell[ , c("unique.cellid", "CTRPv2.cellid")]

# rownames(curationCell) <- curationCell[ , "unique.cellid"]

# curationDrug <- drug_all[which(!is.na(drug_all[ , "CTRPv2.drugid"])),]
# curationDrug <- curationDrug[ , c("unique.drugid", "CTRPv2.drugid","smiles")]
# rownames(curationDrug) <- curationDrug[ , "unique.drugid"]


tmpdir=tempdir()
path.data=file.path("/pfs/out")
path.drug=file.path(path.data, "drug")
path.cell=file.path(path.data, "celline")
path.sens=file.path(path.data, "sens")
saveres=file.path(path.data,"saveres")
verbose=FALSE

options(stringsAsFactors=FALSE)


if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.sens)) { dir.create(path.sens, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(file.path(path.data, "dwl"))) { dir.create(file.path(path.data, "dwl"), showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }

myfn <- "CTRPv2.0_2015_ctd2_ExpandedDataset.zip"
if(!file.exists(file.path(path.data, "dwl", myfn))){
  dwlresult <- download.file(url="ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip", destfile=file.path(path.data, "dwl", myfn), quiet=TRUE)

  res <- unzip(zipfile=file.path(path.data, "dwl", "CTRPv2.0_2015_ctd2_ExpandedDataset.zip"), exdir=file.path(path.data, "dwl"), overwrite=FALSE)

  file.copy(from=file.path(path.data, "dwl", "v20.meta.per_compound.txt"), to=file.path(path.drug, "v20.meta.per_compound.txt"))
  file.copy(from=file.path(path.data, "dwl", "v20.meta.per_cell_line.txt"), to=file.path(path.cell, "v20.meta.per_cell_line.txt"))
  file.copy(from=file.path(path.data, "dwl", "v20.meta.media_comp.txt"), to=file.path(path.cell, "v20.meta.media_comp.txt"))
  file.copy(from=file.path(path.data, "dwl", "v20.data.per_cpd_post_qc.txt"), to=file.path(path.sens, "v20.data.per_cpd_post_qc.txt"))
  file.copy(from=file.path(path.data, "dwl", "v20.data.curves_post_qc.txt"), to=file.path(path.sens, "v20.data.curves_post_qc.txt"))
  file.copy(from=file.path(path.data, "dwl", "v20.meta.per_experiment.txt"), to=file.path(path.sens, "v20.meta.per_experiment.txt"))

}



ctrp.drugs <- read.delim(file.path(path.drug,"v20.meta.per_compound.txt"), sep = "\t", header = TRUE)
ctrp.cells <- read.delim(file.path(path.cell,"v20.meta.per_cell_line.txt"), sep = "\t", header = TRUE)


####################################################
##### Cells
####################################################


# matchToIDTableCELL <- function(ids,tbl, column) {
#   sapply(ids, function(x) {
#     myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
#     if(length(myx) > 1){
#       stop("Something went wrong in curating cell ids")
#     }
#     return(tbl[myx, "unique.cellid"])
#   })
# }

  ctrp.cells[, "cellid"] <- ctrp.cells[,"ccl_name"]

### Using CTRP tissue types for now because they are more precise than COSMIC it seems (verifying with Cellosaurus)
  ctrp.cells[, "tissueid"] <- ctrp.cells$ccle_primary_site

  # matchUnique <- match(tolower(gsub(badchars, "", ctrp.cells[,"ccl_name"])), tolower(gsub(badchars, "", curationCell$CTRPv2.cellid)))
  # ctrp.cells <- ctrp.cells[!is.na(matchUnique),]
  # cell_match <- matchToIDTableCELL(ctrp.cells[,"ccl_name"], curationCell, "CTRPv2.cellid")
  # cell_match <- as.character(unlist(cell_match))
  # ctrp.cells$cellid <- cell_match



####################################################
##### Drugs
####################################################

  # matchToIDTableDRUG <- function(ids,tbl, column) {
  #   sapply(ids,function(x) {
  #     myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
  #     if(length(myx) > 1){
  #       stop("Something went wrong in curating drug ids")
  #     }
  #     return(tbl[myx, "unique.drugid"])
  #   })
  # }
  
  # matchUnique <- match(tolower(gsub(badchars, "", ctrp.drugs[,"cpd_name"])), tolower(gsub(badchars, "", curationDrug$CTRPv2.drugid)))
  # No_drug_match <- which(is.na(matchUnique))
  # No_drug_match <- ctrp.drugs[,"cpd_name"][No_drug_match]
  # print(No_drug_match)
  # curationDrug_temp <- curationDrug
  # curationDrug_temp$CTRPv2.drugid <- tolower(gsub(badchars, "", curationDrug_temp$CTRPv2.drugid))
  # ctrp.drugs_temp <- tolower(gsub(badchars, "", ctrp.drugs[,"cpd_name"]))
  
  # drug_match <- matchToIDTableDRUG(ctrp.drugs_temp, curationDrug_temp, "CTRPv2.drugid")
  # drug_match <- as.character(unlist(unique(drug_match)))
  # drug_match <- append(drug_match, "Tipifarnib-P1", after=218) #accounts for duplicate drug - will be removed below afterwords.
  # drug_match[218] <- "Tipifarnib-P2"
  ctrp.drugs$drugid <- ctrp.drugs[,"cpd_name"]
  
####################################################
##### Viability
####################################################



ctrp.sensitivityRaw <- read.delim(file.path(path.sens,"v20.data.per_cpd_post_qc.txt"))

ctrp.sensitivityInfo <- read.delim(file.path(path.sens,"v20.meta.per_experiment.txt"))

repExps <- unique(ctrp.sensitivityInfo$experiment_id[duplicated(ctrp.sensitivityInfo$experiment_id)])

## Just collapsing the dates, everything else is identical. 
for(exp in repExps){

  myx <- ctrp.sensitivityInfo$experiment_id == exp
  duplicates <- duplicated(ctrp.sensitivityInfo$experiment_id) & myx
  first <- myx & !duplicates
  # print(ctrp.sensitivityInfo[myx,])
  # browser()
  ctrp.sensitivityInfo[first,] <- apply(ctrp.sensitivityInfo[myx,], 2, function(x) paste(unique(x), collapse="//"))

  ctrp.sensitivityInfo <- ctrp.sensitivityInfo[!duplicates,]
}

ctrp.sensitivityInfo[,"cellid"] <- ctrp.cells$cellid[match(ctrp.sensitivityInfo$master_ccl_id, ctrp.cells$master_ccl_id)]
ctrp.sensitivityRaw[,"cellid"] <- ctrp.sensitivityInfo[match(ctrp.sensitivityRaw$experiment_id, ctrp.sensitivityInfo$experiment_id), "cellid"] 
ctrp.sensitivityRaw[,"drugid"] <- ctrp.drugs$drugid[match(ctrp.sensitivityRaw$master_cpd_id, ctrp.drugs$master_cpd_id)] 
ctrp.sensitivityRaw[,"culture_media"] <- ctrp.sensitivityInfo[match(ctrp.sensitivityRaw$experiment_id, ctrp.sensitivityInfo$experiment_id), "culture_media"] 

#### There exist two instances of two different experiment ids with all other annotations exactly the same!!!! Maybe this was an internal control? Should ask[[?
#### Resolved. This was an accidental duplication. 
## need to keep experiment ids in the column for this case?
experimentIds <- paste(ctrp.sensitivityRaw[,"cellid"], ctrp.sensitivityRaw[,"drugid"],ctrp.sensitivityRaw[,"culture_media"], ctrp.sensitivityRaw[,"experiment_id"],sep="_")
ctrp.sensitivityRaw$experimentIds <- experimentIds


# ctrp.sensitivityInfo[,"drugid"] <- drugExps$drugid[match(ctrp.sensitivityInfo$experiment_id, drugExps$experiment_id)]


# repExps <- unique(ctrp.sensitivityInfo$experiment_id[duplicated(ctrp.sensitivityInfo$experiment_id)])

# expIds <- paste( ctrp.sensitivityInfo$cellid, ctrp.sensitivityInfo$drugid, ctrp.sensitivityInfo$culture_media, sep="_")

sensitivityInfo <- data.frame("experimentIds" = unique(experimentIds))

sensitivityInfo[,c("cellid", "drugid","culture_media", "experiment_id")] <- do.call(rbind, strsplit(sensitivityInfo$experimentIds, "_"))


mediaInfo <- read.delim(file.path(path.cell,"v20.meta.media_comp.txt"))

sensitivityInfo[,"media_composition"] <- mediaInfo$media_composition[match(sensitivityInfo$culture_media, mediaInfo$culture_media)]

rownames(sensitivityInfo) <- sensitivityInfo$experimentIds





sensitivityInfo <- data.frame("experimentIds" = unique(experimentIds))

sensitivityInfo[,c("cellid", "drugid","culture_media", "experiment_id")] <- do.call(rbind, strsplit(sensitivityInfo$experimentIds, "_"))

mediaInfo <- read.delim(file.path(path.cell, "v20.meta.media_comp.txt"))

sensitivityInfo[,"media_composition"] <- mediaInfo$media_composition[match(sensitivityInfo$culture_media, mediaInfo$culture_media)]

rownames(sensitivityInfo) <- sensitivityInfo$experimentIds

# all_conc <- unique(ctrp.sensitivityRaw$cpd_conc_umol)



sensRaw.dt <- as.data.table(ctrp.sensitivityRaw[,c("experimentIds", "cpd_conc_umol","cpd_avg_pv")])

library(parallel)



# options("mc.cores"=nthread)

# values <- mclapply(sensitivityInfo$experimentIds[1:10], function(x, sensRaw.dt) {
    
#     data <- sensRaw.dt[experimentIds ==x]

#   # data <- data[order(data$Concentration),]

#     return(data)

#   }, sensRaw.dt=sensRaw.dt)

setorder(sensRaw.dt, experimentIds, cpd_conc_umol)

values <- split(sensRaw.dt, by="experimentIds")


concList <- (unlist(lapply(values, nrow)))

ncon <- max(concList)


sensitivityRaw = array(data=NA_real_, dim=c(length(values), ncon, 2), dimnames=list(names(values), paste(rep("dose", times=ncon), 1:ncon, sep=""), c("Dose","Viability")))

# names(values) <- sensitivityInfo$experimentIds

for(n in 1:length(values)){
  if (n %% 5000 == 0){

    print(paste("Experiment",n,"of",length(values)))

  }

  values[[n]][,experimentIds:=NULL]

  sensitivityRaw[n, 1:(nrow(values[[n]])), ] <- as.matrix(values[[n]])

}

  sensitivityRaw[,,"Viability"] <- sensitivityRaw[,,"Viability"] *100

  save(sensitivityRaw, file=file.path(saveres, "ctrpv2.sens.raw.RData"))
  
  
raw.sensitivity <- sensitivityRaw
                                        
sensitivityInfo[,"Number of Doses Tested"] <- concList[rownames(sensitivityInfo)]
                                        
  
save(raw.sensitivity,sensitivityInfo, ctrp.drugs, ctrp.cells, concList, file="/pfs/out/drug_post.RData")


raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))

dir.create("/pfs/out/slices/")

for(i in seq_along(raw.sensitivity.x)){

  slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/ctrpv2_raw_sens_", i, ".rds"))

}
