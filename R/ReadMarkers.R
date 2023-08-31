ReadMarkers <-
function(markers=NULL){
cat("\n","Choose markers' file (*.csv)","\n")
fname="markers"
if(is.null(markers)) 
{
fname <- file.choose()
markers <- read.csv(file=fname,na.strings = c(".","NA","na","N/A"))
}else
{
if(is.character(markers))
{
fname=markers
markers <- read.csv(file=fname,na.strings =  c(".","NA","na","N/A"))
}else
{
  for (u in 1:ncol(markers))
  {
  nna=(1:length(markers[,u]))[!is.na(markers[,u])]
  si_na=nna[markers[nna,u]=="."]
  if(length(si_na)>0) markers[si_na,u]=NA	
  if(!any(is.na(as.numeric(as.character(markers[!is.na(markers[,u]),u])))))
  markers[,u]=as.numeric(as.character(markers[,u]))
  }
}
}
recode <- function(i){
cod <- ifelse(markers[,i]==0,-1,
                                ifelse(markers[,i]==1,0,
                                        ifelse(markers[,i]==2,1,NA)))
    return(cod)
}
#if(markers[3,1]==3) markers <- markers[,-1]
if(nrow(as.matrix(table(markers[,1])))>4) markers <- markers[,-1]
nmark <- ncol(markers)
markers_vec=as.numeric(as.matrix(markers))
if(any(markers_vec[!is.na(markers_vec)]==2)) codes <- sapply(1:nmark,recode)  else codes <- markers
nmark <- ncol(markers)
colnames(codes) <- paste("M",1:nmark,sep="")
nas <- which(is.na(codes),arr.ind=T)
i=1
while(i <= nrow(nas)){
  codes[nas[i,1],nas[i,2]] <- 0 
  i=i+1
}
cat("\n",paste(fname, "has been read","\n"))
return(as.matrix(codes))
}

