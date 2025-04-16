using MixedModels, CSV, DataFrames, LinearAlgebra, Statistics, GLM, Suppressor

LinearAlgebra.BLAS.set_num_threads(1)



drug = ARGS[1]::String

tissue = ARGS[2]::String

gene = ARGS[3]::String

filePath = ARGS[4]::String

outPath = ARGS[5]::String


nthread = 1::Int64  #40 threads faster than 80 on niagara

modelData = DataFrame(CSV.File(filePath, pool=false, types=Dict("x" => Float64, "dataset" => String)));

select!(modelData, Not(:Column1));


R = min(modelData[1, :R], 10000000);

select!(modelData, Not(:R));

if any(names(modelData) .== "tissueid") && length(unique(modelData[!, :tissueid])) < 2
    select!(modelData, Not(:tissueid))
end


function scale(x::Array{Float64,1})::Array{Float64,1}
    return (x .- mean(x)) / std(x)
end


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

modelData2 = standardizeByDataset(modelData);

m0 = fit(LinearMixedModel, @formula(y ~ (x + 0| dataset) + x + 0), modelData2);
t0 = coef(m0)[1];

numPerLoop = 1000;

t = zeros(R);

i = 1::Int64;
while i <= R
    # jj = 1::Int64;
    @suppress begin
        Threads.@threads for jj = 0:(numPerLoop-1)::Int64
            resampled = getBootSample(modelData);
            nDSS = length(unique(resampled[!,:dataset]));
            if nDSS > one(nDSS)
                    m1 = fit(LinearMixedModel, @formula(y ~ (x + 0| dataset) + x + 0), resampled);
                    t[i+jj] = coef(m1)[1];
                # end
            else
                # @suppress begin
                    m1 = fit(LinearModel, @formula(y ~ x + 0), resampled);
                    t[i+jj] = coef(m1)[1];
                # end
            end
            # global i = i + 1;
            # jj = jj + 1;
        end
    end
    global i = i + numPerLoop;
    if t0 > 0
        numOverZero = sum(t[1:(i-1)] .< 0);
    else 
        numOverZero = sum(t[1:(i-1)] .> 0);
    end
    p_hat = numOverZero / (i-1);
    margin = 2.576 * sqrt(p_hat * (1 - p_hat) / (i - 1)) ## 99% confidence interval
    if p_hat - margin > 0.05
        print("Early stopping at $i out of $R using 99% CI bounding p value larger than 0.05\n")
        break
    end
end


t = t[1:(i-1)];
R = i - 1;

# this takes 8 seconds for 1e4, seems to scale linearly from here. 20x improvement!


badchars = r"[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\^]|[/]|[\\]|[ ]|[(]|[)]"

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


