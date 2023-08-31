data_initialization=function()
{
data(MARKERS_GSI)
data(TRAIT_GSI)
data(WEIGHTS_GSI)
data(WEIGHTS_GSI)
data(QTL_GSI)
write.table(MARKERS_GSI,"00MARKERS_GSI.csv",sep =",",quote =FALSE)
write.table(TRAIT_GSI,"00TRAIT_GSI.csv",sep =",",quote =FALSE)
write.table(WEIGHTS_GSI,"00WEIGHTS_GSI.csv",sep =",",quote =FALSE)
write.table(QTL_GSI,"00QTL_GSI.csv",sep =",",quote =FALSE)

}