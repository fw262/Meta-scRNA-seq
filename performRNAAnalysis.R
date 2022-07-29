#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("DropletUtils"))

 # create parser object
parser <- ArgumentParser()

parser$add_argument("-m", "--expressionMatrix",
                    help="Folder containing filtered_feature_bc_matrix generated from cellranger count.",
                    metavar="directory")
parser$add_argument("-k", "--krakenFile",
                    help="File containing single cell tags with Kraken output",
                    metavar="file")
# variable arguments
parser$add_argument("--minCells", default=3, type="integer",
                    help="Keep features that appear in at least this many cells [default %(default)s].")
parser$add_argument("--minFeatures", default=100, type="integer",
                    help="Keep cells with at least this many features [default %(default)s].")
parser$add_argument("--maxFeatures", default=4000, type="integer",
                    help="Keep cells with at most this many features [default %(default)s].")
parser$add_argument("--maxMito", default=10, type="integer",
                    help="Keep cells with at most this much percentage mitochondrial content [default %(default)s].")
parser$add_argument("--resParam", default=0.2, type="double",
                    help="FindClusters resolution parameter [default %(default)s].")
parser$add_argument("--nPCA", default=10, type="integer",
                    help="Number of PCAs to use in dimensional reduction [default %(default)s].")
parser$add_argument("--nDiffGenes", default=5, type="integer",
                    help="Number of differentially expressed genes to show for each cluster [default %(default)s].")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
print("Loaded libraries and arguments.")

# create Seurat Object with unmapped content
createSeuWithMeta<-function(inputDGE,inputMetaBarcodes,minCellsVal){
  geneSeu<-Read10X(data.dir = dirname(inputDGE))
  # load in metagenomic counts
  krakenBarcode<-read.delim(inputMetaBarcodes, header=FALSE)
  krakenBarcode<-krakenBarcode[!duplicated(krakenBarcode$V8),]
  metaGenome<-as.data.frame.matrix(with(krakenBarcode,table(V3,V6)))
  metaCells<-colnames(metaGenome)[colnames(metaGenome) %in% colnames(geneSeu)]
  metaGenome<-metaGenome[,metaCells]
  colnames(metaGenome)<-gsub(pattern = "-",replacement = ".",colnames(metaGenome))
  seuratData<-rbind.fill(data.frame(geneSeu),metaGenome)
  rownames(seuratData)<-c(rownames(geneSeu),rownames(metaGenome))
  seuratData[is.na(seuratData)]<-0
  # remove empty rownames
  seuratData<-seuratData[rownames(seuratData)!="",]
  seuratData<-seuratData[,colnames(seuratData)!=""]
  # Initialize the Seurat object with the raw (non-normalized data).
  seuObj<-CreateSeuratObject(counts = seuratData, min.cells = minCellsVal)
  # add percent data info
  seuObj[["percent.mt"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")
  validMetaFeat<-rownames(metaGenome)[rownames(metaGenome) %in% rownames(seuObj)]
  seuObj[["percent.unmapped"]] <- PercentageFeatureSet(seuObj, features = validMetaFeat)
  return(seuObj)
}
seuratObj<-createSeuWithMeta(args$expressionMatrix,
                             args$krakenFile,
                             args$minCells)
print("Initiated seurat object.")

# run clustering using all features
runClustering<-function(inputSeu,minFeat=200,maxFeat=4000,maxMit=10,resParam=0.2,maxPCA=10){
  inputSeu <- subset(inputSeu,subset = nFeature_RNA > minFeat & nFeature_RNA < maxFeat & percent.mt < maxMit)
  inputSeu <- NormalizeData(inputSeu)
  inputSeu <- ScaleData(inputSeu, features = rownames(inputSeu))
  inputSeu <- RunPCA(object = inputSeu, features = rownames(inputSeu)) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  inputSeu <- FindNeighbors(object=inputSeu,dims=1:maxPCA)
  inputSeu <- FindClusters(object=inputSeu,resolution=resParam)
  inputSeu <- RunUMAP(inputSeu,dims=1:maxPCA,check_duplicates = F)
  return(inputSeu)
}
seuratObj<-runClustering(seuratObj,
                         minFeat = args$minFeatures,
                         maxFeat = args$maxFeatures,
                         maxMit = args$maxMito,
                         resParam = args$resParam,
                         maxPCA = args$nPCA)
print("Finished dimensional reduction and clustering.")

# find differentially expressed unmapped features
metaMarkers<-FindAllMarkers(seuratObj,features = rownames(seuratObj)[grep("taxid", rownames(seuratObj))])

# find differentially expressed gene markers
geneMarkers<-FindAllMarkers(seuratObj,features = rownames(seuratObj)[-grep("taxid", rownames(seuratObj))])

makeDotPlotMeta<-function(seuratIn,featuresIn,titleIn){
  pOut<-DotPlot(seuratIn,
                features = unique(featuresIn),
                dot.scale = 4)+
    ylab("cell type")+
    theme(axis.title=element_blank(),
          plot.title = element_text(size=10))+
    ggtitle(paste0(titleIn))+
    theme(axis.text.x=element_text(size=8,angle = 45,hjust=1),
          axis.text.y=element_text(size=8))+
    coord_flip()
  return(pOut)
}

# make directory
dir.create(paste0(dirname(args$krakenFile),"/plots"))

# plot percent unmapped
VlnPlot(seuratObj,features="percent.unmapped")
ggsave(filename = "percentUnmapped.pdf",
       device = "pdf",
       path=paste0(dirname(args$krakenFile),"/plots/"))
print("Saved percent unmapped.")

# plot umap
DimPlot(seuratObj)
ggsave(filename = "umap.pdf",
       device = "pdf",
       path=paste0(dirname(args$krakenFile),"/plots/"))
print("Saved umap.")

# plot top N differentially expressed gene markers
if(!is.null(geneMarkers$gene)){
  topNgenes <- geneMarkers %>% group_by(cluster) %>% top_n(n = args$nDiffGenes, wt = avg_log2FC)
  makeDotPlotMeta(seuratObj,unique(topNgenes$gene),titleIn="differentially expressed genes")
  ggsave(filename = "diffGeneDot.pdf",
         device = "pdf",
         path=paste0(dirname(args$krakenFile),"/plots/"))
  print("Saved differentially expressed genes.")
}

metaToexclude<-c("Homo sapiens (taxid 9606)","root (taxid 1)","unclassified (taxid 0)")
# plot differentially expressed unmapped markers
if(!is.null(metaMarkers$gene)){
  potentialMeta<-metaMarkers$gene
  markersToPlot<-potentialMeta[!potentialMeta %in% metaToexclude]
  if(length(markersToPlot)>0){
    makeDotPlotMeta(seuratObj,markersToPlot,titleIn="differentially expressed unmapped features")
    ggsave(filename = "diffMetaDot.pdf",
           device = "pdf",
           path=paste0(dirname(args$krakenFile),"/plots/"))
    print("Saved differentially expressed unmapped features.")
  }
}

# ALL expressed metagenomic features
allMetaFeat<-rownames(seuratObj)[grep("taxid", rownames(seuratObj))]
markersToPlotAll<-allMetaFeat[!allMetaFeat %in% metaToexclude]
if(length(markersToPlotAll)>0){
  makeDotPlotMeta(seuratObj,markersToPlotAll,titleIn="all unmapped features")
  ggsave(filename = "allMetaDot.pdf",
         device = "pdf",
         path=paste0(dirname(args$krakenFile),"/plots/"))
  print("Saved all unmapped features.")
}

# save merged transcriptomic and metagenomic data
mergedMatPath<-paste0(dirname(args$krakenFile),"/Solo.out/merged")
write10xCounts(x = seuratObj@assays$RNA@counts, path = mergedMatPath)

