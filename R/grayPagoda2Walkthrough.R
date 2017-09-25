##############################
# Load PAGODA
##############################

library(GO.db)
library(Rtsne)
library(org.Hs.eg.db)
library(pcaMethods)
library(igraph)
library(irlba)
require(mgcv)
library(WGCNA)
library(Cairo)
library(parallel)
library(Matrix)
## custom installations
library("Rcpp", lib.loc="/usr/local/lib/R/site-library")
library("dbscan", lib.loc="/usr/local/lib/R/site-library")
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
## load from Nick's library for latest
library("pagoda2", lib.loc="/home/barkasn/R/x86_64-pc-linux-gnu-library/3.4/")

##############################
# Data preparation
##############################

# Read data from your file, rows as genes colums as cells
load('../data/grayScde.RData')
oldAnnotation <- groups

load('../data/grayNeuronNPCs.RData')
myCountMatrix <- cd
myCountMatrix <- as(myCountMatrix, 'dgCMatrix') # convert to sparse dgCMatrix
myAnnotation <- sg
head(myCountMatrix)
dim(myCountMatrix)
length(myAnnotation)

##############################
# Process the data
##############################

# Generate a new pagoda2 object
myPagoda2Object <- Pagoda2$new(x = myCountMatrix, n.cores = 10)

# Adjust the variance
myPagoda2Object$adjustVariance(plot = TRUE, gam.k = 10)

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes -- in this case 2000
myPagoda2Object$calculatePcaReduction(nPcs = 20, n.odgenes = 2.e3, maxit=2000)

# Generate K-nearest neighbour graph
myPagoda2Object$makeKnnGraph(k = 30, type = 'PCA', center = TRUE, distance = 'cosine')

##############################
# Identify clusters
##############################

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
myPagoda2Object$getKnnClusters(method = infomap.community, type = 'PCA', name = 'infomap')

##############################
# Generate embeddings
##############################

# Generate an embedding with tSNE on the basis of the PCA reduction
myPagoda2Object$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity=20)

##############################
# Plot the generated embeddings
##############################
par(mfrow=c(1,3))
myPagoda2Object$plotEmbedding(type = 'PCA',
                              clusterType = 'infomap',
                              embeddingType = 'tSNE',
                              mark.clusters = TRUE,
                              alpha = 0.5)
myPagoda2Object$plotEmbedding(type = 'PCA',
                              groups = myAnnotation,
                              embeddingType = 'tSNE',
                              mark.clusters = TRUE,
                              alpha = 0.5)
myPagoda2Object$plotEmbedding(type = 'PCA',
                              groups = oldAnnotation,
                              embeddingType = 'tSNE',
                              mark.clusters = TRUE,
                              alpha = 0.5)

################################
# Pathway overdispersion
################################
go.env <- p2.generate.human.go(myPagoda2Object)
myPagoda2Object$testPathwayOverdispersion(setenv = go.env, verbose = TRUE,
                                          correlation.distance.threshold = 0.1,
                                          recalculate.pca = FALSE,
                                          min.pathway.size = 100, max.pathway.size = 1000)

################################
# Generate the web application
################################

# Metadata to be displayed
additionalMetadata <- list()
additionalMetadata$newCluster <- p2.metadata.from.factor(myPagoda2Object$clusters$PCA$infomap, displayname = 'Infomap', s = 0.7, v = 0.8,start = 0, end = 0.5)
additionalMetadata$scdeCluster <- p2.metadata.from.factor(factor(oldAnnotation), displayname = 'SCDE Annotation', s = 0.7, v = 0.8,start = 0, end = 0.5)
additionalMetadata$originalCluster <- p2.metadata.from.factor(factor(myAnnotation), displayname = 'Original Annotation', s = 0.7, v = 0.8,start = 0, end = 0.5)

# Generate GO genesets for the web app
p2.generate.go.web <- function(gene.names, egALIAS2EG = NULL, egGO2ALLEGS = NULL) {
    require(GO.db)
    require(BiocGenerics)
    require(AnnotationDbi)
    require(parallel)

    if (is.null(egALIAS2EG)) {
        stop('egALIAS2EG cannot be null, it has to be an object like org.Hs.egALIAS2EG');
    }

    if (is.null(org.Hs.egGO2ALLEGS)) {
        stop('org.Hs.egGO2ALLEGS cannot be null it has to be an object like org.Hs.egGO2ALLEGS');
    }

    if (!is.character(gene.names)) {
        stop("gene.names needs to be a character vector of gene names");
    }

    ids <- unlist(mclapply(BiocGenerics::mget(gene.names, egALIAS2EG, ifnotfound = NA), function(x) x[1]))

    ## Swap names and ids
    rids <- names(ids)
    names(rids) <- ids

    ## Get go environment
    go.env <- AnnotationDbi::eapply(egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))

    ## Filter for go terms with at least 5 genes
    go.env <- go.env[unlist(lapply(go.env, length)) > 5]

                                        # TODO make this parallel
    geneSets <- lapply(names(go.env), function(x) {
        list(
            properties = list(
                locked = T,
                genesetname = x,
                shortdescription = GO.db::GOTERM[[x]]@Term
            ),
            genes = c(go.env[[x]])
        )
    })

    ## Name the genesets
    names(geneSets) <- names(go.env)

    ## return geneSets
    return(geneSets)
}
myGeneNames <- colnames(myPagoda2Object$counts)
goSets <- p2.generate.go.web(gene.names = myGeneNames, egALIAS2EG = org.Hs.egALIAS2EG,
                   egGO2ALLEGS = org.Hs.egGO2ALLEGS)


# Generate and display web app
myPagoda2WebObject <-
    make.p2.app(
        myPagoda2Object,
        dendrogramCellGroups = myPagoda2Object$clusters$PCA$infomap,
        additionalMetadata = additionalMetadata,
        geneSets = goSets,
        show.clusters = FALSE, # Hide the clusters that were used for the dendrogram from the metadata
        );

# Show the app
show.app(app=myPagoda2WebObject,name='pagoda_gray')
