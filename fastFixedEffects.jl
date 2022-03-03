using LinearAlgebra, Statistics, DataFrames, StatsBase, InteractiveUtils



# # # ## Current timing, 7 seconds per 5e3, N=320
# permFixedEffect <- function(model.data, R){

# 	# t <- foreach(i =icount(R), .combine=c, .inorder=FALSE) %dorng% {
# 	# 	library(RhpcBLASctl)
# 	# 	RhpcBLASctl::blas_set_num_threads(1)
#  #        x <- model.data[,"x"]
#  #        y <- sample(model.data[,"y"])
#  #        # coop::pcor(x,y)
#  #        cor(x,y)
#  #    }
# 	datasets <- unique(model.data[,"dataset"])
# 	for(ds in datasets){
# 		myx <- model.data[,"dataset"] == ds
# 		model.data[myx,"x"] <- scale(model.data[myx,"x"])
# 	}
# 	ds.vec <- lapply(datasets, function(ds) return(model.data[,"dataset"] == ds))
# 	denom  <- sum(model.data[,"x"]^2)
# 	xt <- model.data[,"x"]/denom
# 	t <- unlist(mclapply(seq_len(R), function(i, model.data, xt, ds.vec) {
#         # model.data[,"x"] <- model.data[,"x"]
#         iny <- sample(model.data[,"y"])
#         # coop::pcor(x,y)
#         y <- standardizeByDatasetInPerm(iny, ds.vec)
#         t0 <- crossprod(xt,y)
# 	}, model.data = model.data, xt = xt, ds.vec = ds.vec, mc.cores=nthread))

#     y <- standardizeByDatasetInPerm(model.data[,"y"], ds.vec)
#     t0 <- crossprod(xt,y)[1]
#     # t0 <- cor(x,y)
#     return(list(t0 = t0, t = t, R = R, pvalue = (sum(abs(t) > abs(t0)) + 1)/(length(t) + 1)))
# }

function scale(x::Array{Float64,1})::Array{Float64,1}
    return (x .- mean(x))/std(x)
end



function scaleWithinDataset(modelData::DataFrame, dataset::String)::DataFrame
    myDS = findall(modelData[!,:dataset].==dataset);
    datasetData = modelData[myDS,:];
    datasetData[!,:x] = scale(datasetData[!,:x]);
    datasetData[!,:y] = scale(datasetData[!,:y]);
    return datasetData
end

## probably the place where I can optimize the most
function standardizeByDatasetInPerm(iny, dsIndx)
	for myDS = dsIndx
		iny[myDS] = scale(iny[myDS]);
	end
	return iny
end





## need to fix the type intability of DataFrame

function permFixedEffectWrapper(modelData::DataFrame, R)

	x = modelData[!,:x]::Array{Float64,1}
	y = modelData[!,:y]::Array{Float64,1}
	dataset = modelData[!,:dataset]
	return permFixedEffect(x,y,dataset,R)
end 

function permFixedEffect(x, y, dataset, R)
	unDatasets = unique(dataset)
    nR = length(y);
    dsIndx = map(x -> findall(dataset.==ds), unDatasets);


	for myDS = dsIndx
		x[myDS] = scale(x[myDS]);
	end

	denom  = sum(x .^2)
 	xt = x ./denom

 	t = zeros(R);



 	for i=1:R
 		yt = sample(y,nR);
 		yt = standardizeByDatasetInPerm(yt, dsIndx);
 		t[i] = xt⋅yt;
 	end 

 	yt = standardizeByDatasetInPerm(y,dsIndx)

 	t0 = xt⋅yt;

 	pvalue = (sum(abs.(t) .> abs(t0))+1)/(R+1)

    return [t0, t, R, pvalue]
end
