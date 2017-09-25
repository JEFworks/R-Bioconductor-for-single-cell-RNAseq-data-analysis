## Install SCDE
require(devtools)
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")
install_github("hms-dbmi/scde")

## Load SCDE
library(scde)

## Load Neuron/NPC dataset from Gray Camp et al.
load('../data/grayNeuronNPCs.RData')

## Filter out lowly expressed genes and poor cells
cd <- clean.counts(cd)
dim(cd)
length(sg)

## Error model
knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = 10, min.count.threshold = 2, min.nonfailed = 5, verbose=1, save.model.plots=FALSE)
valid.cells <- knn$corr.a > 0
table(valid.cells)

## Variance normalize
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 10, plot = TRUE)
## Normalize out depth bias
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

## Get GO annotations
library(liger)
library(GO.db)
go.env <- org.Hs.GO2Symbol.list
desc <- AnnotationDbi::select(GO.db, keys = names(go.env), columns = c("TERM"), multiVals = "CharacterList")
names(go.env) <- paste(names(go.env), desc$TERM)
go.env <- list2env(go.env)
head(ls(go.env))

## Get overdispersed pathways
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 10)
## Get de novo gene sets
clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 10, plot = TRUE)

## Get full info on the top aspects
tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
hc <- pagoda.cluster.cells(tam, varinfo, method='ward.D2')
## determine overall cell clustering
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.6, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0, top=4)

## View
l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(sg)]
names(l2cols) <- names(sg)
groups <- cutree(hc, 5)
col.cols <- rbind(l2cols[hc$labels], rainbow(5)[groups])
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = col.cols)

## Make app
# compile a browsable app, showing top three clusters with the top color bar
app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = col.cols, cell.clustering = hc, title = "gray")
# show app in the browser (port 1468)
show.app(app, "gray", browse = TRUE, port = 1468)

## Plot few markers
markers <- c(
    "SCN2A","GRIK3","CDH6","NRCAM",
    "SOX11",
    "SLC24A2", "SOX4", "DCX", "TUBB3","MAPT",
    "RBFOX3",
    "PTBP2", "PTBP1",  "ZFP36L2",
    "HMGN2", "PAX6", "SFRP1",
    "SOX2", "HES1", "NOTCH2", "CLU","HOPX",
    "MKI67","TPX2",
    "EOMES", "NEUROD4"
    )
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
mat <- mat[markers,]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))

save.image(file="../data/grayScde.RData")
