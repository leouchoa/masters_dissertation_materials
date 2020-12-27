library(class)
library(soiltexture)

necessary_cols <- c("prop1_areia","prop2_argila","prop3_silte","coord_x","coord_y")

# ----- CTB0562 LOOC Estimation -----

ctb0562 <- read.csv("data/alr_transformed_with_locations_UNIQUE.csv", row.names=NULL)

soil_dts_0562 <- ctb0562[,necessary_cols]


aux_tbl <- setNames(soil_dts_0562[,1:3] * 100,c("SAND","CLAY","SILT"))

class_label_0562 <- as.data.frame(TT.points.in.classes(aux_tbl,class.sys= "SiBCS13.TT"))


class_label_0562$classif <- apply(class_label_0562,1,function(x){
 aux_term <- colnames(class_label_0562)[which(as.logical(x))];
 ifelse(length(aux_term) != 1,paste0(aux_term,collapse = "/"),aux_term)
})


ctb0562$label <- as.factor(class_label_0562$classif)

train_cols <- c("label","coord_x","coord_y")

knn_train_dts <- ctb0562[,train_cols]


soildts_knn_results <- vector("list",length = 10)

for(i in seq_along(soildts_knn_results)){
  soildts_knn_results[[i]] <-
    knn.cv(
      train = knn_train_dts[,-1],
      cl = knn_train_dts$label,
      k = i
      )
}

knn_train_conf_mat <- lapply(
  soildts_knn_results,
  function(x){
    table(x,ctb0562$label)
  }
)

data.frame(
  K = seq_along(knn_train_conf_mat),
  "Acurácia" = unlist(
    lapply(knn_train_conf_mat,function(x){
      100 * round(
        sum(diag(x))/nrow(ctb0562),
        3
      )
    })
  )
)



# ----- CTB0809 Validation Test -----


ctb0809 <- read.csv("data/NEW_alr_transformed_with_locations_UNIQUE.csv", row.names=NULL)

soil_dts_0809 <- ctb0809[,necessary_cols]


aux_tbl <- setNames(soil_dts_0809[,1:3] * 100,c("SAND","CLAY","SILT"))

class_label_0809 <- as.data.frame(TT.points.in.classes(aux_tbl,class.sys= "SiBCS13.TT"))


class_label_0809$classif <- apply(class_label_0809,1,function(x){
  aux_term <- colnames(class_label_0809)[which(as.logical(x))];
  ifelse(length(aux_term) != 1,paste0(aux_term,collapse = "/"),aux_term)
})


ctb0809$label <- as.factor(class_label_0809$classif)

test_cols <- c("coord_x","coord_y")

knn_test_dts <- ctb0809[,test_cols]

soildts_knn_final <- vector("list",length = 10)

for(i in seq_along(soildts_knn_final)){
  soildts_knn_final[[i]] <-
    knn(
      train = knn_train_dts[,-1],
      test = knn_test_dts,
      cl = knn_train_dts$label,
      k = i
    )
}


test_set_conf_mat <- 
  lapply(
    soildts_knn_final,
    function(x){
      table(ctb0809$label,x)
    }
  )



test_set_acc <- 
  data.frame(
    K = seq_along(test_set_conf_mat),
    "Acurácia" = 
      unlist(
        lapply(
          test_set_conf_mat,
          function(x){
            sum(diag(x))/nrow(ctb0809)
          }
        )
      )
  )

# ---- Class Probability ----

test_prob <- vector("list",length = 10)

for(j in seq_along(test_prob)){
  test_prob[[j]] <-
    knn(
      train = knn_train_dts[,-1],
      test = knn_test_dts,
      cl = knn_train_dts$label,
      k = j,
      prob = TRUE
    )
}

lapply(
  test_prob,
  function(x){
    # as.character(x)
    # attr(x,"prob")
    data.frame(
      "Observado" = ctb0809$label,
      "Predito" = as.character(x),
      "Probabilidade" = attr(x,"prob")
    )
  }
)

