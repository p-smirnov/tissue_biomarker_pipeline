library(boot)
library(MASS)


test <- mvrnorm(50, c(0,0), matrix(c(1,.4,.4,1),2,2))

bootCor <- function(mat, myx){

	return(cor(mat[myx,])[1,2])

}

boot.out <- boot(data=test, statistic=bootCor, R=1000)




sampleTwoLevelMultinom <- function(nSample = 1000, N=40,M=3,Nm =c(5,15,20)){


	sapply(seq_len(nSample), function(X) {
		outArrayList <- list()
	
	
		for(ii in seq_len(M)){
	
			outArrayList[[ii]] <- numeric(Nm[ii])
		}
	
		sampledBuckets <- sample.int(M, replace=TRUE)
	
		for(bucket in sampledBuckets){
			outArrayList[[bucket]] <- outArrayList[[bucket]] + tabulate(sample.int(Nm[bucket], replace=TRUE), Nm[bucket])
		}
	
		return(do.call(c, outArrayList))
	})
}


sampleMultinomWithFlatProbs <- function(nSample = 1000, N=40,M=3,Nm =c(5,15,20)){


	probVector <- do.call(c,sapply(Nm, function(Nbucket){
		rep(1/M*1/Nbucket, times=Nbucket)
	}))

	rmultinom(nSample, N, probVector)
}

sampleMultinomWithTwoLevelFlatProbs <- function(nSample = 1000, N=40,M=3,Nm =c(5,15,20)){

	sapply(seq_len(nSample), function(X) {
	sampledBuckets <- sample.int(M, replace=TRUE)



	probVector <- list()
	probVector <- sapply(Nm, function(N)return(numeric(N)))
	for(bucket in sampledBuckets){
		probVector[[bucket]] <- 1/M*1/Nm[bucket] + probVector[[bucket]]
	}
	probVector <- do.call(c, probVector)

	rmultinom(1, sum(Nm[sampledBuckets]), probVector)})
}
