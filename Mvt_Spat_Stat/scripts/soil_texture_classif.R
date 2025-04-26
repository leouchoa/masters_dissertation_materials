library(soiltexture)

alr_transformed_with_locations_UNIQUE <- read.csv("~/Documents/git/Mvt_Spat_Stat/scripts/discussion_misc/data/alr_transformed_with_locations_UNIQUE.csv", row.names=NULL)

soil_dts <- alr_transformed_with_locations_UNIQUE[,1:4]


aux_tbl <- setNames(BivMaternEstim:::alr_inv(soil_dts[,1:2]) * 100,c("SAND","CLAY","SILT"))

class_label <- as.data.frame(TT.points.in.classes(aux_tbl,class.sys= "SiBCS13.TT"))


# class_label$classif <- 

apply(class_label,1,function(x){
  aux_term <- colnames(class_label)[which(as.logical(x))];
  ifelse(length(aux_term) != 1,paste0(aux_term,collapse = "/"),aux_term)
})

class_label$classif <- aux_res

TT.plot( class.sys = "SiBCS13.TT" , lang = "pt")
