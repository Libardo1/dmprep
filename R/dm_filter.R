# Design Matrix Filtering Function from the Design Matrix Preparation Package
#    Copyright (C) 2016  Benjamin R. Fitzpatrick
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# The package author may be contacted via electronic mail at <ben.r.fitzpatrick@gmail.com>

#' Filter a design matrix to reduce its collinearity.
#'
#' \code{dm_filter} Filters the provided design matrix to reduce its collinearity.  All pairs of covariates that have correlation coefficient magnitudes greater than the value specified to \code{mpccm} are extracted.  These overtly correlated pairs of covariates are then searched to identify which pairs contain the covariates named in the first element of the character vector \code{prefer}.  All covariates that are overtly correlated with the first element of the list of covariates to prefer are then filtered from the design matrix.  This proceedure is then repeated for each subsequent member of the list of covariates to preferentially retain.  Once this proceedure has been repeated for each covariate on the list of covariates to preferentially retain a simplified version of the filtering is conducted complete the filtering of the design matrix.  A covariate is choosen at random from the first pair remaining on the list of overtly correlated pairs of covariates and any covariates overtly correlated with this covariate are filtered from the design matrix.  This proceedure is then repeated until there are no longer any pairs of covariates among those included in the filtered design matrix which have correlation coefficient magnitudes greater than the value supplied to \code{mpccm}.
#'
#' @param DM The design matrix, each row is an observation, each column is a covariate.
#' @param prefer An ordered character vector of covariates names.  The first element will be retained by the filtering above all other covariates with which it is correlated. The second element will be retained after the first etc.
#' @param mpccm The maximum permitted correlation coefficient magnitude between remaining covariate pairs post filtering.
#' @param na.action What to do if there are NAs among the elements of the design matrix supplied to `DM'.
#' @return The output is a two element list consisting of the named elements \code{drop} and \code{FDM}.  The first element of this list, \code{drop}, is a character vector containing the names of all the covariates that were filtered from the design matrix due having correlation coefficient magnitude greater than the cut off value with covariates that have been retained in the filtered design matrix.  The second element of this list, \code{FDM}, is the design matrix that has been filtered to enforce the specified maximum permitted correlation coefficient between remaining covariates.
#'
#' @examples
#' # Filter the design matrix DM:
#' X.F <- dm_filter(DM = DM, prefer = c('band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'NDVI'), mpccm = 0.8, na.action = 'Ignore')
#' 
#' The covariates that have been filtered from the design matrix:
#' X.F[['drop']]
#'
#' The filtered design matrix:
#' X.F[['FDM']]

dm_filter <- function(DM = DM,
                  prefer = c('band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'NDVI'),
                  mpccm = 0.8,
                      na.action = c('Ignore', 'Halt')){
  # it appear 'raster' also has a select function so we will need to use package::function syntax to be safe
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
  VCC.df <- dplyr::filter(.data = C.df, !(Var1 == Var2) & abs(Cor) > mpccm) # Very Correlated Covariates
  drop <- vector(mode = 'character')
  pref.l <- length(prefer)
  for(i in 1:pref.l){
    drop.i.v2 <- dplyr::filter(.data = VCC.df, Var1 == prefer[i]) %>% dplyr::select(Var2)
    drop.i.v1 <- dplyr::filter(.data = VCC.df, Var2 == prefer[i]) %>% dplyr::select(Var1)
    drop <- unique(c(drop, drop.i.v2[,1], drop.i.v1[,1]))
    prefer <- prefer[!(prefer %in% drop)]
    VCC.df <- dplyr::filter(.data = VCC.df, !(Var1 %in% drop) & !(Var2 %in% drop))
  }
  DM.d1 <- dplyr::select(.data = DM, -one_of(drop))
  cmat.d1 <- data.frame(cor(DM.d1, use = 'complete.obs'))
  cmat.d1$Var1 <- row.names(cmat.d1)
  C.d1.df <- gather(data = cmat.d1, value = Cor, key = Var2, -Var1)
  C.d1.df.F <- dplyr::filter(.data = C.d1.df, !(Var1 == Var2) & abs(Cor) > mpccm)
  # final drop at random (if necessary)
  if( nrow(C.d1.df.F) > 0 ){
    while( (max(abs(C.d1.df.F$Cor)) > mpccm)){
      col.i <- sample(x = 1:2, size = 1) 
      drop.i <- dplyr::select(.data = C.d1.df.F, col.i) %>% slice(.,1)
      drop <- unique(c(drop, drop.i[,1]))
      C.d1.df.F <- dplyr::filter(.data = C.d1.df.F, !(Var1 %in% drop) & !(Var2 %in% drop))
      if( (nrow(C.d1.df.F) == 0)  ){
          break()
      }
    }
  }
  DM.d2 <- dplyr::select(.data = DM, -one_of(drop))
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
