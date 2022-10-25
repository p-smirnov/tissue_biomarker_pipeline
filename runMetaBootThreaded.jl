using MixedModels, CSV, DataFrames, LinearAlgebra, Statistics, GLM, Suppressor

LinearAlgebra.BLAS.set_num_threads(1)




drug = ARGS[1]

tissue = ARGS[2]

gene = ARGS[3]

filePath = ARGS[4]

outFileName = ARGS[5]


# nthread = 1::Int64  #40 threads faster than 80 on niagara

modelData = DataFrame(CSV.File(filePath, pool=false));

select!(modelData, Not(:Column1));


R = min(modelData[1,:R], 10000000);

select!(modelData, Not(:R));


function scale(x::Array{Float64,1})::Array{Float64,1}
    return (x .- mean(x))/std(x)
end


function sampleWithinDataset(modelData::DataFrame, dataset)::DataFrame
    myDS = findall(modelData[!,:dataset].==dataset);
    nDS = length(myDS);
    myx = rand(1:nDS, nDS);
    myDS = myDS[myx];
    datasetData = modelData[myDS,:];
    if any(names(datasetData).=="tissueid")
        tissues = unique(modelData[!,:tissueid])
        for tissue = tissues
            datasetData=removeTissueMean!(datasetData, tissue);
        end
    end
    datasetData[!,:x] = scale(datasetData[!,:x]);
    datasetData[!,:y] = scale(datasetData[!,:y]);
    return datasetData
end


function removeTissueMean!(datasetData::DataFrame, tissue)::DataFrame
    myTissue = findall(datasetData[!,:tissueid].==tissue);
    tissueMeanX = mean(datasetData[myTissue,:x]::Array{Float64,1});
    tissueMeanY = mean(datasetData[myTissue,:y]::Array{Float64,1});
    datasetData[myTissue,:x] = datasetData[myTissue,:x] .- tissueMeanX;
    datasetData[myTissue,:y] = datasetData[myTissue, :y] .- tissueMeanY;
    return datasetData
end



function scaleWithinDataset(modelData::DataFrame, dataset)::DataFrame
    myDS = findall(modelData[!,:dataset].==dataset);
    datasetData = modelData[myDS,:];
    if any(names(datasetData).=="tissueid")
        tissues = unique(modelData[!,:tissueid])
        for tissue = tissues
            datasetData=removeTissueMean!(datasetData, tissue);
        end
    end
    datasetData[!,:x] = scale(datasetData[!,:x]);
    datasetData[!,:y] = scale(datasetData[!,:y]);
    return datasetData
end


function getBootSample(modelData::DataFrame)::DataFrame
    sampledDatasets = rand(unique(modelData[!,:dataset]), 
                           length(unique(modelData[!,:dataset])));
    resampled = map(x -> sampleWithinDataset(modelData,x), sampledDatasets);
    resampled = reduce(append!, resampled);
    return resampled
end

function standardizeByDataset(modelData::DataFrame)
    sampledDatasets = unique(modelData[!,:dataset]);
    standardized = map(x -> scaleWithinDataset(modelData,x), sampledDatasets);
    standardized = reduce(append!, standardized);
    return standardized
end


t = zeros(R);

Threads.@threads for i = 1:R::Int64
    resampled = getBootSample(modelData);
    nDSS = length(unique(resampled[!,:dataset]));
    if nDSS > one(nDSS)
        # @suppress begin
            m1 = fit(LinearMixedModel, @formula(y ~ (x + 0| dataset) + x + 0), resampled);
            t[i] = coef(m1)[1];
        # end
    else
        # @suppress begin
            m1 = fit(LinearModel, @formula(y ~ x + 0), resampled);
            t[i] = coef(m1)[1];
        # end
    end
end

# this takes 8 seconds for 1e4, seems to scale linearly from here. 20x improvement!
# threaded on 40 threads, performance pan-cancer is 100000 in 94 seconds.
modelData2 = standardizeByDataset(modelData);

m0 = fit(LinearMixedModel, @formula(y ~ (x + 0| dataset) + x + 0), modelData2);
t0 = coef(m0)[1];

# badchars = r"[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\^]|[/]|[\\]|[ ]|[(]|[)]"

# tissueClean = replace(tissue, badchars => s".")

# drugClean = replace(drug, badchars => s".")

outfile = outFileName;

open(outfile, "w") do f
    println(f, "t0:");
    println(f, t0);
    println(f, "N:");
    println(f, nrow(modelData));
    println(f, "R:");
    println(f, R);
    println(f, "t:");
    for i in t
        println(f, i)
    end
end


# PD.0325901 Lung ENSG00000130477


