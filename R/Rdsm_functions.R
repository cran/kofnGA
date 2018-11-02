# These functions are taken fromthe Rdsm package by N. Matloff, with unneeded parts
# removed.

# To avoid getting notes about global variables in CRAN check
if(getRversion() >= "2.15.1") utils::globalVariables("myinfo")

#' @keywords internal
getidxs <- function(m) {
  splitIndices(m,myinfo$nwrkrs)[[myinfo$id]]
}


#' @keywords internal
mgrinit <- function(cls) {
  # set up so that each worker node will have a global variable myinfo
  # that contains the thread ID and number of threads
  setmyinfo <- function(i,n) {
    assign("myinfo",list(id = i,nwrkrs = n),pos=tmpenv)
  }
  ncls <- length(cls)
  clusterEvalQ(cls,tmpenv <- new.env())
  clusterApply(cls,1:ncls,setmyinfo,ncls)
  clusterEvalQ(cls,myinfo <- get("myinfo",tmpenv))
  # we create global variables only at the workers, thus OK for CRAN,
  # but CRAN check complains anyway, so here is a workaround
  clusterEvalQ(cls,gbl <- globalenv())
  # send the threads needed Rdsm functions
  clusterExport(cls,"getidxs",envir=environment())
}


#' @keywords internal
mgrmakevar <- function(cls,varname,nr,nc,vartype="double") {
  tmp <- big.matrix(nrow=nr,ncol=nc,type=vartype)
  # make accessible to manager
  assign(varname,tmp,pos=parent.frame())
  # get the descriptor for this big.matrix object, to send to the
  # worker nodes
  clusterExport(cls,"varname",envir=environment())
  desc <- describe(tmp)
  clusterExport(cls,"desc",envir=environment())
  clusterEvalQ(cls, tmp <- attach.big.matrix(desc))
  clusterEvalQ(cls,assign(varname,tmp))
  invisible(0)
}


