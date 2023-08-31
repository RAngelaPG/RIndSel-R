Read <- function (file.dat = NULL, file.wgt = NULL, design=design,mis =  c(".","NA","na","N/A"), one.env=one.env, genomic=FALSE, markers=markers, block.ex, softR){ 
  if( design=="lattice"){
	if(one.env) fcol <<- 4 else fcol <<- 5
  }
  
  if( design=="rcb"){
	if(one.env) fcol <<- 3 else fcol <<- 4
  }
  
  if(design=="ACBD"){
	if(block.ex==TRUE){
		if(softR=="AsremlR"){
			if(one.env) fcol <<- 7 else fcol <<- 8
		}else{
			if(one.env) fcol <<- 11 else fcol <<- 12
		}
	}else{
		if(softR=="AsremlR"){
			if(one.env) fcol <<- 6 else fcol <<- 7
		}else{
			if(one.env) fcol <<- 10 else fcol <<- 11
		}
	}
  }
  
  traits1=file.dat
  
  nom <- colnames(traits1[, -(1:(fcol - 1))])
  n <- length(nom)
  
  if(!is.null(file.wgt)){
	wgts <<- read.csv(file = file.wgt, header = T)  
	wgts <<-as.data.frame(wgts)
  }
  ntraits <<- ncol(traits1) - (fcol - 1)  
  if(genomic){
    markers_actualizados<<- ReadMarkers(markers)
    markers =markers_actualizados
  }
  
  return(traits1)
}

