alr <- function(compositional_no_coords){


  alr_low_level <- function(x){

    return(
      c( log(x[1]/x[3]), log(x[2]/x[3]) )
    )

  }

  return(
    setNames(
      as.data.frame(
        t(
          apply(
            compositional_no_coords,
            1,
            alr_low_level
          )
        )
      ),
      paste0("alr_comp_0",1:2)
    )
  )

}
