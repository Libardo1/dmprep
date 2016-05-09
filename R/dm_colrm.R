
# test data:
Test.Mat <- data.frame(matrix(data = c(1:4,NA,6:10,11:14,NA,16:20,rep(0,5),rep(NA,5)), nrow = 5))

# remove columns that contain only NA values:

Test.Mat[,!apply(X = Test.Mat, MARGIN = 2, FUN = anyNA)]


# remove columns that contain entries of exactly 0 in all rows:

all.zero <- function(x){
    x.naom <- na.omit(x)
    if(length(x.naom) == 0){
      return(FALSE)  
    }
    x.u <- unique(x, na.rm = TRUE)
    if(length(x.u) > 1){
      return(FALSE)
    }
    if((length(x.u) == 1) & x.u == 0){
      return(TRUE)
    }
    if((length(x.u) == 1) & !(x.u == 0)){
      print('non zero unique value in all of this column')  
      return(FALSE)
    }

}



Test.Mat[,!apply(X = Test.Mat, MARGIN = 2, FUN = all.zero)]
