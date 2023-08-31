ReadQTL <-
function(qtl=NULL){
if(is.null(qtl)) 
{
qtl0 <- choose.files(,"QTL file")
    qtl0<- read.csv(file=qtl0,skip=1)
    qtl <- qtl0[,seq(1,ncol(qtl0)-1,by=2)]
    link <- qtl0[,seq(2,ncol(qtl0),by=2)]
    
}else
{
if(is.character(qtl))
{
qtl0 <- qtl
    qtl0<- read.csv(file=qtl0,skip=1)
    qtl <- qtl0[,seq(1,ncol(qtl0)-1,by=2)]
    link <- qtl0[,seq(2,ncol(qtl0),by=2)]

}else
{

colnames(qtl)=as.matrix(qtl[1,])
qtl=qtl[-1,]

  for (u in 1:ncol(qtl))
  {
  nna=(1:length(qtl[,u]))[!is.na(qtl[,u])]
  si_na=nna[qtl[nna,u]==""]
  if(length(si_na)>0) qtl[si_na,u]=NA	
  if(!any(is.na(as.numeric(as.character(qtl[!is.na(qtl[,u]),u])))))
  qtl[,u]=as.numeric(as.character(qtl[,u]))
  }
qtl0=qtl
qtl <- qtl0[,seq(1,ncol(qtl0)-1,by=2)]
link <- qtl0[,seq(2,ncol(qtl0),by=2)]

}
}

    return(list(QTL=qtl,Link=link))
}

