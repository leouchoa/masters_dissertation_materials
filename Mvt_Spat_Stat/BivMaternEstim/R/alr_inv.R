alr_inv <- function(gp_pred_with_no_coords){


  alr_inv_low_level <- function(gp_pred){

    aux <- exp(gp_pred) / (1 + sum(exp(gp_pred)))

    return(
      c(aux,1 - aux[1] - aux[2])
    )

  }

  return(
    setNames(
      as.data.frame(
        t(
          apply(
            gp_pred_with_no_coords,
            1,
            alr_inv_low_level
          )
        )
      ),
      paste0("comp_0",1:3)
    )
  )

}
