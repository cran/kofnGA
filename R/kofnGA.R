# To avoid getting notes about global variables in CRAN check
if(getRversion() >= "2.15.1") utils::globalVariables("OFVALS")


#*****************************************************************************************
#*** kofnGA ******************************************************************************
#*****************************************************************************************

#' Search for the best subset of size k from n choices.
#'
#' \code{kofnGA} implements a genetic algorithm for subset selection.  The function searches
#' for a subset of a fixed size, k, from the integers 1:n, such that user-supplied function
#' \code{OF} is minimized at that subset.  The selection step is done by tournament selection
#' based on ranks, and elitism may be used to retain the best solutions from one generation
#' to the next. Population objective function values can be evaluated in parallel.
#'
#' \itemize{
#' \item Tournament selection involves creating mating "tournaments" where two groups of
#'   \code{tourneysize} solutions are selected at random without regard to fitness.  Within each
#'   tournament, victors are chosen by weighted sampling based on within-tournament fitness ranks
#'   (larger ranks given to more fit individuals). The two victors become parents of an offspring.
#'   This process is carried out \code{popsize} times to produce the new population.
#' \item Crossover (reproduction) is carried out by combining the unique elements of both parents and
#'   keeping \code{k} of them, chosen at random.
#' \item Increasing \code{tourneysize} will put more "selection pressure" on the choice of mating pairs,
#'   and will speed up convergence (to a local optimum) accordingly. Smaller \code{tourneysize} values will
#'   conversely promote better searching of the solution space.
#' \item Increasing the size of the elite group (\code{keepbest}) also promotes more homogeneity in the
#'  population, thereby speeding up convergence.
#' }
#'
#' @section Notes on parallel evaluation:
#'
#' Specifying a \code{cluster} allows \code{OF} to be evaluated over the population in parallel.
#' The population of solutions will be distributed among the workers in \code{cluster} using
#' static dispatching.  Any cluster produced by \code{\link[parallel]{makeCluster}}
#' should work, though the \code{sharedmemory} option is only appropriate for a cluster of
#' workers on the same multicore processor.
#'
#' Solutions must be sent to workers (and results gathered back) once per generation. This
#' introduces communication overhead.  Overhead can be reduced, but not eliminated, by setting
#' \code{sharedmemory=TRUE}. The impact of parallelization on run time will depend on
#' how the run time cost of evaluationg \code{OF} compares to the communication overhead.  Test runs
#' are recommended to determine if parallel execution is beneficial in a given situation.
#'
#' Note that only the objective function evaluations are distributed to the workers.  Other parts of
#' the algorithm (mutation, crossover) are computed serially. As long as \code{OF} is deterministic,
#' reproducibility of the results from a given random seed should not be affected by the use of
#' parallel computation.
#'
#' @section References:
#'
#' Mark A. Wolters (2015), \dQuote{A Genetic Algorithm for Selection of Fixed-Size Subsets,
#' with Application to Design Problems,} \emph{Journal of Statistical Software}, volume 68,
#' Code Snippet 1, \href{http://www.jstatsoft.org/article/view/v068c01}{available online}.
#'
#' @param n The maximum permissible index (i.e., the size of the set we are doing subset
#' selection from).  The algorithm chooses a subset of integers from 1 to \code{n}.
#' @param k The number of indices to choose (i.e., the size of the subset).
#' @param OF The objective function.  The first argument of OF should be an index vector of
#' length \code{k} containing integers in the range [1, \code{n}].  Additional arguments can
#' be passed to \code{OF} through \code{\dots}.
#' @param popsize The size of the population; equivalently, the number of offspring produced
#' each generation.
#' @param keepbest The \code{keepbest} least fit offspring each generation are replaced by the
#' \code{keepbest} most fit members of the previous generation. Used to implement elitism.
#' @param ngen The number of generations to run.
#' @param tourneysize The number of individuals involved in each tournament at the selection stage.
#' @param mutprob The probability of mutation for each of the \code{k} chosen indices in each
#' individual. An index chosen for mutation jumps to any other unused index, uniformly at random.
#' This probability can be set indirectly through \code{mutfrac}.
#' @param mutfrac The average fraction of offspring that will experience at least one mutation.
#' Equivalent to setting \code{mutprob} to \code{1 - (1 - mutfrac)^(1/k)}.  Only used if
#' \code{mutprob} is not supplied. This method of controlling mutation may be preferable
#' if the algorithm is being run at different values of \code{k}.
#' @param initpop A \code{popsize}-by-\code{k} matrix of starting solutions.  The final populations
#' from one GA search can be passed as the starting point of the next search.  Possibly useful
#' if using this function in an adaptive, iterative, or parallel scheme (see examples).
#' @param verbose An integer controlling the display of progress during search.  If \code{verbose}
#' takes positive value \code{v}, then the iteration number and best objective function value are
#' displayed at the console every \code{v} generations. Otherwise nothing is displayed. Default is
#' zero (no display).
#' @param cluster If non-null, the objective function evaluations for each generation are done in
#' parallel.  \code{cluster} can be either a cluster as produced by \code{\link[parallel]{makeCluster}},
#' or an integer number of parallel workers to use. If an integer, \code{makeCluster(cluster)} will
#' be called to create a cluster, which will be stopped on function exit.
#' @param sharedmemory If \code{cluster} is non-null and \code{sharedmemory} is \code{TRUE}, the
#' parallel computation will employ shared-memory techniques to reduce the overhead of repeatedly
#' passing the population matrix to worker threads. Uses code from the \code{Rdsm} package, which
#' depends on \code{bigmemory}.
#' @param \dots Additional arguments passed to \code{OF}.
#'
#' @return A list of S3 class "GAsearch" with the following elements:
#'   \item{bestsol}{A vector of length \code{k} holding the best solution found.}
#'   \item{bestobj}{The objective function value for the best solution found.}
#'   \item{pop}{A \code{popsize}-by-\code{k} matrix holding the final population, row-sorted
#'     in order of increasing objective function. Each row is an index vector representing
#'     one solution.}
#'   \item{obj}{The objective function values corresponding to each row of pop.}
#'   \item{old}{A list holding information about the search progress. Its elements are:}
#'   \item{old$best}{The sequence of best solutions known over the course of the search
#'     (an \code{(ngen+1)}-by-\code{k} matrix)}
#'   \item{old$obj}{The sequence of objective function values corresponding to the solutions
#'     in \code{old$best}}
#'   \item{old$avg}{The average population objective function value over the course of the
#'     search (a vector of length \code{ngen+1}).  Useful to give a rough indication of
#'     population diversity over the search.  If the average fitness is close to the best
#'     fitness in the population, most individuals are likely quite similar to each other.}
#'
#' @importFrom stats rbinom
#' @importFrom utils tail
#' @importFrom parallel clusterExport makeCluster stopCluster clusterEvalQ clusterApply splitIndices parRapply
#' @importFrom bigmemory big.matrix describe attach.big.matrix as.matrix
#'
#' @export
#' @keywords optimize design
#' @seealso \code{\link{plot.GAsearch}} plot method, \code{\link{print.GAsearch}} print
#' method, and \code{\link{summary.GAsearch}} summary method.
#' @examples
#' #---Find the four smallest numbers in a random vector of 100 uniforms---
#' # Generate the numbers and sort them so the best solution is (1,2,3,4).
#' Numbers <- sort(runif(100))
#' Numbers[1:6]                                              #-View the smallest numbers.
#' ObjFun <- function(v, some_numbers) sum(some_numbers[v])  #-The objective function.
#' ObjFun(1:4, Numbers)                                      #-The global minimum.
#' out <- kofnGA(n = 100, k = 4, OF = ObjFun, ngen = 50, some_numbers = Numbers)  #-Run the GA.
#' summary(out)
#' plot(out)
#'
#' \dontrun{
#' # Note: the following two examples take tens of seconds to run (on a 2018 laptop).
#'
#' #---Harder: find the 50x50 principal submatrix of a 500x500 matrix s.t. determinant is max---
#' # Use eigenvalue decomposition and QR decomposition to make a matrix with known eigenvalues.
#' n <- 500                                         #-Dimension of the matrix.
#' k <- 50                                          #-Size of subset to sample.
#' eigenvalues <- seq(10, 1, length.out=n)          #-Choose the eigenvalues (all positive).
#' L <- diag(eigenvalues)
#' RandMat <- matrix(rnorm(n^2), nrow=n)
#' Q <- qr.Q(qr(RandMat))
#' M <- Q %*% L %*% t(Q)
#' M <- (M+t(M))/2                                  #-Enusre symmetry (fix round-off errors).
#' ObjFun <- function(v,Mat)  -(determinant(Mat[v,v],log=TRUE)$modulus)
#' out <- kofnGA(n=n, k=k, OF=ObjFun, Mat=M)
#' print(out)
#' summary(out)
#' plot(out)
#'
#' #---For interest: run GA searches iteratively (use initpop argument to pass results)---
#' # Alternate running with mutation probability 0.05 and 0.005, 50 generations each time.
#' # Use the same problem as just above (need to run that first).
#' mutprob <- 0.05
#' result <- kofnGA(n=n, k=k, OF=ObjFun, ngen=50, mutprob=mutprob, Mat=M)  #-First run (random start)
#' allavg <- result$old$avg                       #-For holding population average OF values
#' allbest <- result$old$obj                      #-For holding population best OF values
#' for(i in 2:10) {
#'   if(mutprob==0.05) mutprob <- 0.005 else mutprob <- 0.05
#'   result <- kofnGA(n=n, k=k, OF=ObjFun, ngen=50, mutprob=mutprob, initpop=result$pop, Mat=M)
#'   allavg <- c(allavg, result$old$avg[2:51])
#'   allbest <- c(allbest, result$old$obj[2:51])
#' }
#' plot(0:500, allavg, type="l", col="blue", ylim=c(min(allbest), max(allavg)))
#' lines(0:500, allbest, col="red")
#' legend("topright", legend=c("Pop average", "Pop best"), col=c("blue", "red"), bty="n",
#'        lty=1, cex=0.8)
#' summary(result)
#' }
#'
kofnGA <- function(n,k,OF,popsize=200,keepbest=floor(popsize/10),ngen=500,
                   tourneysize=max(ceiling(popsize/10),2),mutprob=0.01, mutfrac=NULL,
                   initpop=NULL, verbose=0, cluster=NULL, sharedmemory=FALSE, ...) {


##=== Basic input checking =====================================================================##
if (!is.null(mutfrac)) {
  if (missing(mutprob)) {
    stopifnot(mutfrac >= 0, mutfrac <= 1)
    mutprob <- 1 - (1 - mutfrac)^(1/k)
  } else {
    warning("Both mutprob and mutfrac supplied. Ignoring mutfrac.")
  }
}
stopifnot(n %% 1 == 0, n > 0, n >= k,
          k %% 1 == 0, k > 0,
          popsize %% 1 == 0, popsize >= 2, popsize > keepbest, popsize > tourneysize,
          keepbest %% 1 == 0, keepbest >= 0,
          ngen %% 1 == 0, ngen >= 1,
          tourneysize %% 1 == 0, tourneysize >= 2,
          mutprob >= 0, mutprob <= 1,
          all(dim(initpop) == c(popsize,k)),
          verbose %% 1 == 0)


##=== Create needed objects ====================================================================##
indices <- 1:n                                   #-The set of possible objects to choose from.
if (keepbest>0) {
    elitespots <- 1:keepbest                     #-Portion of pop rows reserved for elites.
    newspots <- (keepbest+1):popsize             #-Portion of pop rows reserved for offspring.
}
fitness.old <- vector(mode="numeric",length=popsize)
fitness.new <- vector(mode="numeric",length=popsize)
offspring <- matrix(0,nrow=popsize,ncol=k)
old <- list()
old$best <- matrix(0,nrow=ngen+1,ncol=k)
old$obj <- vector(mode="numeric",length=ngen+1)
old$avg <- vector(mode="numeric",length=ngen+1)
Tourneys1 <- matrix(0,nrow=popsize,ncol=tourneysize)
Tourneys2 <- matrix(0,nrow=popsize,ncol=tourneysize)
if (is.null(cluster)) {
  useparallel = FALSE
  if (sharedmemory) {
    warning("cluster is NULL. Ignoring sharedmemory=TRUE.")
    sharedmemory <- FALSE
  }
} else {
  useparallel = TRUE
  if (!any(class(cluster) == "cluster")) {
    cluster <- makeCluster(cluster)
    on.exit(stopCluster(cluster))
  }
}

##=== Initialize the population ================================================================##
# When no population supplied, initialize the population using randomly-chosen indices.
if (is.null(initpop))
  pop <- t(replicate(popsize,sample(indices,k)))
else
  pop <- initpop


##=== Create function for evaluating population fitness ========================================##
# If not using parallel execution, just loop over the population, calling OF.
# For parallel case without shared memory, use the standard parRapply available in parallel pkg.
# For shared memory case, variables POP and OFVALS are created and subsequently modified in-place.
if (!useparallel) {
  getfitness <- function(P) {
    f <- vector(mode="numeric", length=popsize)
    for (i in 1:popsize) f[i] <- OF(P[i,], ...)
    f
  }
}
if (useparallel) {
  list(...)   #Evaluate any promises (see vignette("parallel"), sec. 10.4)
  if (!sharedmemory) {
    fcn <- function(v) OF(v, ...)
    clusterExport(cluster, "OF", envir=environment())
    getfitness <- function(P) {
      parRapply(cluster, P, fcn)
    }
  } else {
    mgrinit(cluster)
    mgrmakevar(cluster, "POP", popsize, k)
    mgrmakevar(cluster, "OFVALS", popsize, 1)
    POP[,] <- pop
    getOFvals <- function(P,result) {
      ix <- getidxs(length(result))
      for (i in ix) {
        result[i] <- OF(P[i,], ...)
      }
      0  # dummy return value
    }
    clusterExport(cluster, c("OF","getOFvals"), envir=environment())
    getfitness <- function(P) {
      POP[,] <- P
      clusterEvalQ(cluster, getOFvals(POP,OFVALS))
      0  # dummy return value
    }
  }
}


##=== Evaluate initial fitness ================================================================##
if (sharedmemory) {
  getfitness(pop)
  fitness.old <- as.matrix(OFVALS)
} else {
  fitness.old = getfitness(pop)
}
old$best[1,] = sort(pop[rank(fitness.old,ties.method="random")==1,])
old$obj[1] = min(fitness.old)
old$avg[1] = mean(fitness.old)


##=== Loop through generations =================================================================##
for (gen in 1:ngen) {

    # Choose mating pairs (selection)-------------------------------------------------------------
    # Matrix Tourneys1 is popsize-by-tourneysize.  Row i contains indices of the individuals
    # chosen to be in the tournament to become parent 1 of the ith offspring.  Matrix Tourneys2 is
    # created in the same way.  In each tournament, the probability of a competing individual
    # being chosen is proportional to its fitness rank in the tournament (larger ranks given to
    # more fit individuals).  This is done using pickfun().  The results, Parents1 and Parents2,
    # are vectors of length popsize containing the row-indices of parents to take from pop.
    Tourneys1[,] <- t(replicate(popsize,sample(1:popsize,tourneysize)))
    Tourneys2[,] <- t(replicate(popsize,sample(1:popsize,tourneysize)))
    pickfun <- function(v) sample(v,1,prob=rank(-fitness.old[v]))
    Parents1 <- apply(Tourneys1,1,pickfun)
    Parents2 <- apply(Tourneys2,1,pickfun)

    # Create offspring (crossover)----------------------------------------------------------------
    # For crossover, just combine the unique elements of both parents, then keep k of them, chosen
    # at random.
    for (i in 1:popsize) {
        combo <- unique(c(pop[Parents1[i],],pop[Parents2[i],]))
        offspring[i,] <- sample(combo,k)
    }

    # Perform mutation----------------------------------------------------------------------------
    # Logical matrix chosen has the same dimensions as pop, with TRUE elements indicating
    # locations chosen for mutation.  Then each offspring in turn is mutated by replacing its
    # chosen elements with a randomly-selected unused index.
    chosen <- matrix(as.logical(rbinom(popsize*k,1,mutprob)),nrow=popsize)
    nchosen <- apply(chosen,1,sum)
    for (i in 1:popsize) {
        if (nchosen[i]>0) {
            toadd <- sample(indices[-offspring[i,]],nchosen[i])
            offspring[i,chosen[i,]] = toadd
        }
    }

    # Evaluate fitness----------------------------------------------------------------------------
    # Put the fitness values for the offspring in a separate vector fitness.new.  The next
    # generation will be chosen based on both fitness.old and fitness.new, since the keepbest
    # elites from the previous generation will stay in.
    if (sharedmemory) {
      getfitness(offspring)
      fitness.new <- as.matrix(OFVALS)
    } else {
      fitness.new <- getfitness(offspring)
    }

    # Form the next generation, keeping the elites------------------------------------------------
    # If keepbest=0, then no elitism.  Just set pop and fitness.old to the offspring values.
    # If keepbest>0, use logical vectors old.keep and new.keep to identify the keepbest most fit
    # members of the previous generation and the (popsize-keepbest) most fit members of the
    # offspring.  Then overwrite pop and fitness.old appropriately.
    if (keepbest==0) {
        pop <- offspring
        fitness.old <- fitness.new
    }
    else {
        old.keep <- rank(fitness.old,ties.method="random") <= keepbest
        new.keep <- rank(fitness.new,ties.method="random") <= popsize-keepbest
        pop[elitespots,] = pop[old.keep,]
        pop[newspots,] = offspring[new.keep,]
        fitness.old[elitespots] = fitness.old[old.keep]
        fitness.old[newspots] = fitness.new[new.keep]
    }

    # Record progress of the search---------------------------------------------------------------
    old$best[gen+1,] <- sort(pop[rank(fitness.old,ties.method="random")==1,])
    old$obj[gen+1] <- min(fitness.old)
    old$avg[gen+1] <- mean(fitness.old)

    # Give output to console, if requested--------------------------------------------------------
    if (verbose>0 && gen%%verbose==0) {
      cat("Finished iteration ", gen, ". Best OF value = ", old$obj[gen+1], "\n")
    }
}


##=== Package the outputs ======================================================================##
# Note, when returning the overall best, need to look at the history, since elitism might be
# turned off (then final population might not contain best solution found).
out <- list()
# Return the record of all runs.
out$old <- old
# Return the final population, sorted.
ord <- order(fitness.old)
out$pop <- matrix(0,nrow=popsize,ncol=k)
for (i in 1:popsize) {
    out$pop[i,] = sort(pop[ord[i],])
}
# Return the objfun values for the final population.
out$obj <- fitness.old[ord]
# Return the best solution found, and its objfun value.
alltimebest <- which(old$obj==min(old$obj,na.rm=TRUE))
alltimebest <- tail(alltimebest,1)               #-In case multiple bests, take last one.
out$bestsol <- out$old$best[alltimebest,]
out$bestobj <- out$old$obj[alltimebest]

class(out) <- "GAsearch"
out

}
