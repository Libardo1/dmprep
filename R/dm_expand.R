# aiapt = Add Interaction And Polynomial Terms

# note this doesn't recenter and rescale created polynomial and interaction terms either your modelling function will have to do that or you should do it yourself before supplying the output of this function to your modelling function

dm_expand <- function(LME = DM,
                  mpo = 4){
                   
  # LME = Linear Main Effects design matrix (as dataframe)
  # this function will only create interactions of linear main effects
  mpo = 4 # note the facility to alter this value would require editing of polynomial term filtering code below (currently a value of 4 is hard coded in) e.g. |2|3|4
  if( (max(abs((colMeans(LME)))) < 1e-10) | (max(colSums(LME^2)) > 1) ){
    print('You have not recentred and rescaled your design matrix prior to requesting the calculation of polynomial and interactions terms. I\'ll do this for you now, recentring to column means of zero and column L1 norms of 1.')
    X <- LME
    X.R <- scale(x = X, center = colMeans(X), scale = FALSE)
    X.RR <- scale(x = X.R, center = FALSE, scale = sqrt(colSums(X.R^2)))
    LME <- X.RR
  }
  Data = LME
  n.LME = ncol(Data)
  # Expand Desing Matrix to include Polynomial Terms up to order mpo for all covariates in LME
  P.ng = expand.grid(colnames(Data), as.character(2:mpo))
  colnames(P.ng) = c('Var','P')
  PDM = data.frame(matrix(0, nrow = nrow(Data), ncol = nrow(P.ng))) # Polynomial Terms Dataframe
  P.n = character(nrow(P.ng))
  for(i in 1:nrow(P.ng)){
    P.n[i] = paste(P.ng[i,'Var'], P.ng[i,'P'], sep = '.')
  }
  colnames(PDM) = P.n
  for(i in 1:nrow(P.ng)){
    PDM[,paste(P.ng[i,'Var'], P.ng[i,'P'], sep = '.')] = (Data[,paste(P.ng[i,'Var'])])^(as.numeric(paste(P.ng[i,'P'])))
  }
  # Expand Design Matrix to include Columns for Interactions to Order = 1/2 Max Polynomial Order
  # Order 2 Interactions: Products of two distinct linear terms
  IO2.2L.i = t(combn(x = 1:n.LME, m = 2))
  IO2.2L.n = character(nrow(IO2.2L.i))
  IO2.2L.DM = data.frame(matrix(0, nrow = nrow(Data), ncol = length(IO2.2L.n)))
  for(i in 1:nrow(IO2.2L.i)){
    IO2.2L.n[i] = paste(colnames(Data)[IO2.2L.i[i,1]], colnames(Data)[IO2.2L.i[i,2]], sep = '_X_')
  }
  colnames(IO2.2L.DM) = IO2.2L.n
  for(i in 1:ncol(IO2.2L.DM)){
    cols = strsplit(x = colnames(IO2.2L.DM)[i], split = '_X_')[[1]]
    IO2.2L.DM[,i] = Data[,cols[1]]*Data[,cols[2]]
  }
  Data.P4I2.DM = data.frame(Data, IO2.2L.DM, PDM)
  return(Data.P4I2.DM)
}
