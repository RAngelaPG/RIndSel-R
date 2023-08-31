IndexName <-
function (op = NULL) 
{
    if (is.null(op)) {
        cat("Select an index", "\n")
        cat("(1) Smith                                                            ", 
            "\n")
        cat("(2) Restrictive Kempthorne & Nordskog Selection Index                ", 
            "\n")
        cat("(3) Restrictive Eigen Selection Index Method                         ", 
            "\n")
        cat("(4) Eigen Selection Index Method                                     ", 
            "\n")            
        cat("(5) Lande and Thompson Index                                         ",
            "\n")
        cat("(6) Molecular Eigen Selection Index                                  ",
            "\n") 

        cat("(7) Genomic Eiegen Selection Index                                   ",
            "\n") 
        cat("(8) Restrictive Kempthorne & Nordskog Genomic Selection Index        ",
            "\n") 
        cat("(9) Restrictive Genomic Eiegen Selection Index                       ",
            "\n") 
        cat("(10) Smith Genomic Selection Index                                   ",
            "\n") 
        cat("(11) Tallis predetermined proportional gains Eiegen Selection Index  ",
            "\n") 
        cat("(12) Tallis predetermined proportional gains Genomic Selection Index ",
            "\n") 

        cat("Selection index option?:", "\n")
        op <- scan(what = numeric(1), nmax = 1, quiet = TRUE)
    }
    SelInd <- switch(op, SmithIndex(),KNIndex(),RESIMIndex(),ESIMIndex(),LTIndex(),MESIMIndex(),GESIMIndex(),KNGSIIndex(),RGESIMIndex(),SGSIIndex(),TGESIMIndex(),TGSIIndex())
    return(SelInd)
}

