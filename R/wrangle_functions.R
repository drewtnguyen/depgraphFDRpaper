tp_ratio_wrt_BH <- function(confarr){
    ratio = confarr[,"tp"] / confarr["BH", "tp"]
    # is NA either if BH made no rejections, or there are no non-nulls. 
    # Leave the NAs as is
    ratio[is.nan(ratio)] = NA
    ratio
}

rej_ratio_wrt_BH <- function(confarr){
    ratio = (confarr[,"tp"] + confarr[,"fp"]) / (confarr["BH","tp"] + confarr["BH","fp"])
    # is NA either if BH made no rejections, or there are no non-nulls. 
    # Leave the NAs as is
    ratio[is.nan(ratio)] = NA
    ratio
}


fdp <- function(confarr){
    fdp = confarr[,"fp"]/(confarr[, "tp"] + confarr[, "fp"])
    # If no rejections, call the fdp 0
    fdp[is.nan(fdp)] = 0
    fdp
}

pow <- function(confarr){
    pow = confarr[,"tp"] / (confarr[, "tp"] + confarr[, "fn"])
    # is NA if there are no non-nulls. 
    # Leave the NAs as is
    pow[is.nan(pow)] = NA
    pow
}

runtime <- function(confarr){
    if("runtime" %in% names(confarr[1,])){
        confarr[, "runtime"]
    }
    else{
        0*confarr[,"fp"]
    }
}



