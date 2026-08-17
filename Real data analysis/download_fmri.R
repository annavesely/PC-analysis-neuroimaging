require(httr)
require(jsonlite)
require(dplyr)
require(stringr)
require(RNifti)


# --------------------------------------------
# DOWNLOAD T-MAPS
# --------------------------------------------

collection <- 2447

url <- paste0("https://neurovault.org/api/collections/", collection, "/images/?limit=100")
res <- GET(url)

images <- fromJSON(content(res, "text", encoding="UTF-8"), flatten=TRUE)$results
images <- as.data.frame(images)

names(images)

tmaps <- images[images$contrast_definition == "Left Hand > Right Hand", ]

dir.create("~/neurovault_2447", showWarnings=FALSE)

for(i in seq_len(nrow(tmaps))){
  download.file(
    tmaps$file[i],
    destfile=paste0("~/neurovault_2447/tmap_", i, ".nii.gz"),
    mode="wb"
  )
}




# --------------------------------------------
# CREATE MASK
# --------------------------------------------

maps <- lapply(1:10, function(i) {
  RNifti::readNifti(paste0("~/neurovault_2447/tmap_", i, ".nii.gz"))
})

mask <- Reduce("+", lapply(maps, function(x) abs(x) > 0)) > 0
sum(mask)

RNifti::writeNifti(mask, "~/neurovault_2447/mask.nii.gz")

