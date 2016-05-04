# pacakge name: prepdm (Prepare Design Matrix)
# ccaaf = covariate collinearity assessment and filtering

# assess the collinearity between covariates
# filter to enforce a maxium permitted correlation coefficient magnitude between covariate pairs

dm_filter <- function(DM = DM,
                  prefer = c('band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'NDVI'),
                  mpccm = 0.8,
                  na.action = c('Ignore', 'Halt')){
   # DM = design matrix  
   # prefer = character vector of covariates to preferentially keep, order is order of preference 
   # mpccm = maximum permitted correlation coefficient magnitude
   # na.action = how to handle NAs
  require(dplyr) # these need to go, use 
  require(tidyr) # package::function() 
  require(ggplot2) # instead of require
  if(na.action == 'Halt'){    
    if( !(nrow(na.omit(DM)) == nrow(DM)) ) {
      stop('DM contains NAs')
    }
  }  
  cmat <- data.frame(cor(DM, use = 'complete.obs'))
  cmat$Var1 <- row.names(cmat)
  C.df <- gather(data = cmat, value = Cor, key = Var2, -Var1)
  # filer out all pairs which we couldnt calculate a correlation for due to missing values
  if(na.action == 'Ignore'){
    C.df <- na.omit(C.df)
  }
  if(na.action == 'Halt'){
    if(!(nrow(na.omit(C.df)) == nrow(C.df))){
      stop('NAs in correlation dataframe')
    }
  }
  # filter out the rows for the correlation of covariates with themselves             
  VCC.df <- filter(.data = C.df, !(Var1 == Var2) & abs(Cor) > mpccm) # Very Correlated Covariates
  # keep <- vector(mode = 'character')
  drop <- vector(mode = 'character')
  # keep[1] <- prefer[1]
  pref.l <- length(prefer)
  for(i in 1:pref.l){
    
    drop.i.v2 <- filter(.data = VCC.df, Var1 == prefer[i]) %>% select(Var2)

    drop.i.v1 <- filter(.data = VCC.df, Var2 == prefer[i]) %>% select(Var1)
    
    drop <- unique(c(drop, drop.i.v2[,1], drop.i.v1[,1]))

    prefer <- prefer[!(prefer %in% drop)]

    VCC.df <- filter(.data = VCC.df, !(Var1 %in% drop) & !(Var2 %in% drop))
    
  }

  DM.d1 <- select(.data = DM, -one_of(drop))
  
  cmat.d1 <- data.frame(cor(DM.d1, use = 'complete.obs'))
  cmat.d1$Var1 <- row.names(cmat.d1)
  C.d1.df <- gather(data = cmat.d1, value = Cor, key = Var2, -Var1)
  C.d1.df.F <- filter(.data = C.d1.df, !(Var1 == Var2) & abs(Cor) > mpccm)
  # final drop at random (if necessary)
  if( nrow(C.d1.df.F) > 0 ){
    while( (max(abs(C.d1.df.F$Cor)) > mpccm)){
      col.i <- sample(x = 1:2, size = 1) 
      drop.i <- select(.data = C.d1.df.F, col.i) %>% slice(.,1)
      drop <- unique(c(drop, drop.i[,1]))
      C.d1.df.F <- filter(.data = C.d1.df.F, !(Var1 %in% drop) & !(Var2 %in% drop))
      if( (nrow(C.d1.df.F) == 0)  ){
          break()
      }
    }
  }
  DM.d2 <- select(.data = DM, -one_of(drop))
  Output <- list(drop = drop, FDM = DM.d2)
  return(Output)
}


## Test:

# example data:

# Data <- read.csv(file = '~/unwgbdos/Data/CSV/Ground_Truthed_Crop_plus_raw_bands_VIs.csv')

# colnames(Data)

# DM <- select(.data = Data, -grep(pattern = '^X$|^crop$', x = colnames(Data)))

# Test <- ccaaf(DM = DM, prefer = c('band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'NDVI'), mpccm = 0.8, na.action = c( 'Halt'))

# colnames(DM)[!colnames(DM) %in% Test$drop]

# head(DM[,colnames(DM)[!colnames(DM) %in% Test$drop]])

# image(abs(cor(DM[,colnames(DM)[!colnames(DM) %in% Test$drop]])))
