using MixedModels, CSV, DataFrames, LinearAlgebra, Statistics, GLM, Suppressor

LinearAlgebra.BLAS.set_num_threads(1)




drug = ARGS[1]::String

tissue = ARGS[2]::String

gene = ARGS[3]::String

filePath = ARGS[4]::String

outPath = ARGS[5]::String


nthread = 1::Int64  #40 threads faster than 80 on niagara

modelData = DataFrame(CSV.File(filePath, pool=false));

select!(modelData, Not(:Column1));


R = min(modelData[1,:R], 10000000);
select!(modelData, Not(:R));


function scale(x::Array{Float64,1})::Array{Float64,1}
    return (x .- mean(x))/std(x)
end


function sampleWithinDataset(modelData::DataFrame, dataset::String)::DataFrame
    myDS = findall(modelData[!,:dataset]::Array{String,1}.==dataset);
    nDS = length(myDS);
    myx = rand(1:nDS, nDS);
    myDS = myDS[myx];
    datasetData = modelData[myDS,:];
    datasetData[!,:x] = scale(datasetData[!,:x]);
    datasetData[!,:y] = scale(datasetData[!,:y]);
    return datasetData
end


function scaleWithinDataset(modelData::DataFrame, dataset::String)::DataFrame
    myDS = findall(modelData[!,:dataset]::Array{String,1}.==dataset);
    datasetData = modelData[myDS,:];
    datasetData[!,:x] = scale(datasetData[!,:x]);
    datasetData[!,:y] = scale(datasetData[!,:y]);
    return datasetData
end


function getBootSample(modelData::DataFrame)
    sampledDatasets = rand(unique(modelData[!,:dataset]::Array{String,1}), 
                           length(unique(modelData[!,:dataset]::Array{String,1})));
    resampled = map(x -> sampleWithinDataset(modelData,x), sampledDatasets);
    resampled = reduce(append!, resampled);
    return resampled
end

function standardizeByDataset(modelData::DataFrame)
    sampledDatasets = unique(modelData[!,:dataset]::Array{String,1});
    standardized = map(x -> scaleWithinDataset(modelData,x), sampledDatasets);
    standardized = reduce(append!, standardized);
    return standardized
end


t = zeros(R);

for i = 1:R::Int64
    resampled = getBootSample(modelData);
    nDSS = length(unique(resampled[!,:dataset]));
    if nDSS > one(nDSS)
        @suppress begin
            m1 = fit(LinearMixedModel, @formula(y ~ (x + 0| dataset) + x + 0), resampled);
            t[i] = coef(m1)[1];
        end
    else
        @suppress begin
            m1 = fit(LinearModel, @formula(y ~ x + 0), resampled);
            t[i] = coef(m1)[1];
        end
    end
end

# this takes 8 seconds for 1e4, seems to scale linearly from here. 20x improvement!
modelData2 = standardizeByDataset(modelData);

m0 = fit(LinearMixedModel, @formula(y ~ (x + 0| dataset) + x + 0), modelData2);
t0 = coef(m0)[1];

badchars = r"[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\^]|[/]|[\\]|[ ]"

tissueClean = replace(tissue, badchars => s".")

drugClean = replace(drug, badchars => s".")

outfile = outPath * "/metaBootRes_" * gene  *"_"* drugClean *"_"* tissueClean *"_out.txt" ;

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


