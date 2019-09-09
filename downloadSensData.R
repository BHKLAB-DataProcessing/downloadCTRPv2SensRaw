library(downloader)
library(PharmacoGxPrivate)
library(devtools)


options(stringsAsFactors=FALSE)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]|[(]|[)]"

cell_all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))

curationCell <- cell_all[which(!is.na(cell_all[ , "CTRPv2.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "CTRPv2.cellid")]

rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[ , "CTRPv2.drugid"])),]
curationDrug <- curationDrug[ , c("unique.drugid", "CTRPv2.drugid","smiles")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]


ctrp.drugs <- read.delim("/pfs/ctrpv2raw/v20.meta.per_compound.txt", sep = "\t", header = TRUE)
ctrp.cells <- read.delim("/pfs/ctrpv2raw/v20.meta.per_cell_line.txt", sep = "\t", header = TRUE)


####################################################
##### Cells
####################################################


matchToIDTableCELL <- function(ids,tbl, column) {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating cell ids")
    }
    return(tbl[myx, "unique.cellid"])
  })
}

  ctrp.cells[, "cellid"] <- NA

### Using CTRP tissue types for now because they are more precise than COSMIC it seems (verifying with Cellosaurus)
  ctrp.cells[, "tissueid"] <- ctrp.cells$ccle_primary_site

  matchUnique <- match(tolower(gsub(badchars, "", ctrp.cells[,"ccl_name"])), tolower(gsub(badchars, "", curationCell$CTRPv2.cellid)))
  ctrp.cells <- ctrp.cells[!is.na(matchUnique),]
  cell_match <- matchToIDTableCELL(ctrp.cells[,"ccl_name"], curationCell, "CTRPv2.cellid")
  cell_match <- as.character(unlist(cell_match))
  ctrp.cells$cellid <- cell_match



####################################################
##### Drugs
####################################################

  matchToIDTableDRUG <- function(ids,tbl, column) {
    sapply(ids,function(x) {
      myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
      if(length(myx) > 1){
        stop("Something went wrong in curating drug ids")
      }
      return(tbl[myx, "unique.drugid"])
    })
  }
  
  matchUnique <- match(tolower(gsub(badchars, "", ctrp.drugs[,"cpd_name"])), tolower(gsub(badchars, "", curationDrug$CTRPv2.drugid)))
  No_drug_match <- which(is.na(matchUnique))
  No_drug_match <- ctrp.drugs[,"cpd_name"][No_drug_match]
  print(No_drug_match)
  curationDrug_temp <- curationDrug
  curationDrug_temp$CTRPv2.drugid <- tolower(gsub(badchars, "", curationDrug_temp$CTRPv2.drugid))
  ctrp.drugs_temp <- tolower(gsub(badchars, "", ctrp.drugs[,"cpd_name"]))
  
  drug_match <- matchToIDTableDRUG(ctrp.drugs_temp, curationDrug_temp, "CTRPv2.drugid")
  drug_match <- as.character(unlist(unique(drug_match)))
  drug_match <- append(drug_match, "Tipifarnib-P1", after=218) #accounts for duplicate drug - will be removed below afterwords.
  drug_match[218] <- "Tipifarnib-P2"
  ctrp.drugs$drugid <- drug_match
  
####################################################
##### Viability
####################################################



ctrp.sensitivityRaw <- read.delim("/pfs/ctrpv2raw/v20.data.per_cpd_post_qc.txt")

ctrp.sensitivityInfo <- read.delim("/pfs/ctrpv2raw/v20.meta.per_experiment.txt")

repExps <- unique(ctrp.sensitivityInfo$experiment_id[duplicated(ctrp.sensitivityInfo$experiment_id)])


for(exp in repExps){

  myx <- ctrp.sensitivityInfo$experiment_id == exp
  duplicates <- duplicated(ctrp.sensitivityInfo$experiment_id) & myx
  first <- myx & !duplicates
  
  ctrp.sensitivityInfo[first,] <- apply(ctrp.sensitivityInfo[myx,], 2, function(x) paste(unique(x), collapse="//"))

  ctrp.sensitivityInfo <- ctrp.sensitivityInfo[!duplicates,]
}

ctrp.sensitivityInfo[,"cellid"] <- ctrp.cells$cellid[match(ctrp.sensitivityInfo$master_ccl_id, ctrp.cells$master_ccl_id)]
ctrp.sensitivityRaw[,"cellid"] <- ctrp.sensitivityInfo[match(ctrp.sensitivityRaw$experiment_id, ctrp.sensitivityInfo$experiment_id), "cellid"] ## matches are sometimes duplicated, but the values for cells and drugs are unique either way, so take the first match
ctrp.sensitivityRaw[,"drugid"] <- ctrp.drugs$drugid[match(ctrp.sensitivityRaw$master_cpd_id, ctrp.drugs$master_cpd_id)] 
ctrp.sensitivityRaw[,"culture_media"] <- ctrp.sensitivityInfo[match(ctrp.sensitivityRaw$experiment_id, ctrp.sensitivityInfo$experiment_id), "culture_media"] 

#### There exist two instances of two different experiment ids with all other annotations exactly the same!!!! Maybe this was an internal control? Should ask[[?

experimentIds <- paste(ctrp.sensitivityRaw[,"cellid"], ctrp.sensitivityRaw[,"drugid"],ctrp.sensitivityRaw[,"culture_media"], ctrp.sensitivityRaw[,"experiment_id"],sep="_")
ctrp.sensitivityRaw$experimentIds <- experimentIds


# ctrp.sensitivityInfo[,"drugid"] <- drugExps$drugid[match(ctrp.sensitivityInfo$experiment_id, drugExps$experiment_id)]


# repExps <- unique(ctrp.sensitivityInfo$experiment_id[duplicated(ctrp.sensitivityInfo$experiment_id)])

# expIds <- paste( ctrp.sensitivityInfo$cellid, ctrp.sensitivityInfo$drugid, ctrp.sensitivityInfo$culture_media, sep="_")

sensitivityInfo <- data.frame("experimentIds" = unique(experimentIds))

sensitivityInfo[,c("cellid", "drugid","culture_media", "experiment_id")] <- do.call(rbind, strsplit(sensitivityInfo$experimentIds, "_"))


mediaInfo <- read.delim("/pfs/ctrpv2raw/v20.meta.media_comp.txt")

sensitivityInfo[,"media_composition"] <- mediaInfo$media_composition[match(sensitivityInfo$culture_media, mediaInfo$culture_media)]

rownames(sensitivityInfo) <- sensitivityInfo$experimentIds


load("/pfs/ctrpv2raw/ctrp_raw.RData")
  
  
raw.sensitivity <- ctrp.sensitivityRaw
                                        
#sensitivityInfo[,"Number of Doses Tested"] <- concList[rownames(sensitivityInfo)]
                                        
  
#save(raw.sensitivity,sensitivityInfo, ctrp.drugs, ctrp.cells, file="/pfs/out/drug_post.RData")


raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))

dir.create("/pfs/out/slices/")

for(i in seq_along(raw.sensitivity.x)){

  slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/ctrpv2_raw_sens_", i, ".rds"))

}
