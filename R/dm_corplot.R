

# function to produce plots to depict the strengths of correlations between all possible pairs of covariates in a design matrix:

# raster plot is the easy option

# option to turn covariate labels on or off (if you do lots of covariates small enough font on all labels to fit along axis will not be legible)


dm_corplot(

DM = X.5km.WS.LSO.F$FDM

cmat <- data.frame(abs(cor(DM, use = 'complete.obs')))
cmat$Var1 <- row.names(cmat)
C.df <- gather(data = cmat, value = Cor, key = Var2, -Var1)

nrow(C.df)    

p <- ggplot(aes(x = Var1, y = Var2, fill = Cor), data = C.df)

p + geom_raster() + coord_equal() + scale_fill_distiller(palette = 'Spectral', limits = c(0,0.8), na.value = 'white') + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(vjust = 0.5), panel.background = element_rect(fill = 'white')) + labs(x = '', y = '')

# chord diagram is harder - however chord diagram easier to interpret

# option to switch between abs cor and cor

# mpccm cor scale cut off
