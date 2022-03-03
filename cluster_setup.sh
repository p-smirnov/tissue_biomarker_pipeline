#!/bin/bash

## -- Authenticate with Azure (requires manual user input)
az login

# -- Read in configuration files
source .azure
source .aks

# -- Create a Blob storage container

# Get the first key for the specified storage account
stgkey=$(
    az storage account keys list \
        -g $resgroup \
        -n $stgacct \
        --query "[0].value"
)

# Create the container, fails if it exists
az storage container create \
    --resource-group $resgroup \
    --account-name $stgacct \
    --account-key $stgkey \
    --name $containername

# -- Upload data to blob container

# Create tmp directory, if it doesn't exist
[ ! -d tmp ] && mkdir tmp

# Fetch the data
git clone \
    https://github.com/snakemake/snakemake-tutorial-data.git \
    tmp/snakemake-tutorial-data

# Upload to blob storage
az storage blob upload-batch \
    -d $containername \
    --account-name $stgacct \
    --account-key $stgkey \
    -s tmp/snakemake-tutorial-data/data/ \
    --destination-path data

# -- Configure and create AKS cluster

# Provision the cluster
az aks create \
    --resource-group $resgroup \
    --name $clustername \
    --vm-set-type VirtualMachineScaleSets \
    --load-balancer-sku standard \
    --enable-cluster-autoscaler \
    --node-count $startnodes \
    --min-count $minnodes \
    --max-count $maxnodes \
    --node-vm-size $vmsize

# Get the login credentials, saved in .kube/config
az aks get-credentials \
    --resource-group $resgroup \
    --name $clustername

# Show basic cluster information
kubectl cluster-info