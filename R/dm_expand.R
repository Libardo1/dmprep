# Design Matrix Expansion Function from the Design Matrix Preparation Package
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

#' Expand a design matrix to include pairwise interactions terms and polynomial terms.
#'
#' \code{dm_expand} expands the supplied design matrix of linear main effects terms for covariates to include columns for all possible pairwise interactions of covariates and columns for polynomial terms for all covariates up to the polynomial order specified to \code{mpo}.  The additional columns for interaction terms are named following the convention where the column named 'first.covariate.name_X_second.covariate.name' contains the products of the observations in 'first.covariate.name' and 'second.covariate.name'.  The additional columns for polynomial terms are named following the convention where the polynomial order is appended to the name of the consituent covariate such that the column named 'covariate.name.2' contains the values for the quadratic (polynomial order 2) terms for the covariate called 'covariate.name'. Note: to avoid one of the consitutent covariates dominating the values of the interaction term all linear main effects terms are recentered to have means of zero and L1 norms of one prior to calculating the pairwise interactions.  Similarly to avoid extremeley large magnitude values in the column for higher order polynomial terms these are also calculate from linear main effects terms that have been recentered to have means of zero and L1 norms of one prior to calculating the pairwise interactions.
#'
#' @param LME The design matrix, each row is an observation, each column is a covariate.  The covariates should include linear main effects only.
#' @param mpo the maximum polynomial order to include.
#' @return The expanded design matrix.
#'
#' @examples
#' # Expand the design matrix DM: 
#' X.P4I2 <- dm_expand(LME = DM, mpo = 4)

dm_expand <- function(LME = DM,
                  mpo = 4){
  # this function will only create interactions of linear main effects
  if( (max(abs((colMeans(LME)))) > 1e-10) | ((max(colSums(LME^2)) - 1) > 1e-10) ){
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
  # Currrently all that is implemented is the facility to include Order 2 Interactions:
  # Products of two distinct linear terms
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
