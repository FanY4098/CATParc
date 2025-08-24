##' MSA Data Preprocessing
##'
##' This function one hot encodes the original protein MSA sequence data.
##'
##' @param data The MSA data should be a N*m matrix, with N sequences and m positions.
##' Each amino acid is encoded by a number from 1 to 20, while 0 represents the gap.
##' @return Returns an list with \item{data1hot}{The one hot encoded data. Note that the columns
##' representing gaps are removed, so are columns with only 0 entries. The standardized
##' data are reported in the end.}
##' \item{Gsz}{A vector of length m. The binary vectors belonging to the same position form
##' a group, and each entry of Gsz is the group size for the associated position. }
##' @export


OneHot_Protein=function(data) {
  data=apply(data,2, as.numeric)
  data1hot=NA
  Gsz=c()
  for (i in 1:dim(data)[2]) {
    mtn=matrix(rep(0,21*dim(data)[1]),nrow=dim(data)[1],ncol=21)
    rows <- 1:dim(data)[1]
    cols <- data[,i]+1
    mtn[cbind(rows, cols)] <- 1
    check_mtn=apply(mtn,2,sum) #need to get rid of the columns with all 0's
    if (length(which(check_mtn==0))>0) mtn=mtn[,-which(check_mtn==0)]
    if (dim(mtn)[2]>2) {
    mtn=mtn[,-1] #remove gaps
    Gsz[i]=dim(mtn)[2] #record the group size of each position
    } else {
      mtn=mtn[,-1] #remove gaps
      Gsz[i]=1
    }
    data1hot=cbind(data1hot,mtn)
  }

  data1hot=data1hot[,-1] #remove NA
  data1hot=scale(data1hot, center = TRUE, scale = TRUE)

  return(list(data1hot=data1hot,Gsz=Gsz))
}


