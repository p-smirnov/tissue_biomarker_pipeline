#/bin/bash


# -- Authenticate to Azure
az login


# -- Load configuration files
source .azure
source .aks
stgkey=$(
    az storage account keys list \
        -g $resgroup \
        -n $stgacct \
        --query "[0].value"
)



# -- Configure Blob storage access for
export AZ_BLOB_ACCOUNT_URL="https://${stgacct}.blob.core.windows.net"
export AZ_BLOB_CREDENTIAL="$stgkey"
export containername=$containername

snakemake --kubernetes --no-shared-fs \
    --default-remote-provider AzBlob \
    --default-remote-prefix $containername \
    --envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_CREDENTIAL containername \
    --container-image "bhklab/aks-petr" \
    --jobs 1 getSigGenes