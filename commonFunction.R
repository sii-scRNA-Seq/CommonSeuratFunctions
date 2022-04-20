#This is a collection of functions useful for Seurat objects.

ProportionPlot <- function(object, name1, name2){
	#Tris function creates Proportion plot in % for two variables in the seurat object 
	#Example use: ProportionPlot(PBMC, "cluster", "Celltype") 
	
	#Get number of cells per cluster
	slot.metadata = slot(object, "meta.data")
	clust.table.group <- table(getElement(slot.metadata,name1), getElement(slot.metadata,name2))

	#Convert to proportion/percentage of total cells in each sample
	clust.prop.group <- as.data.frame(prop.table(x = clust.table.group, margin = 2))

	#Plot as stacked bar plot
	ggplot(data=clust.prop.group, aes(x=Var2, y=Freq, fill=Var1)) + theme_linedraw() +
	geom_bar(stat="identity", color="black") + labs(x=name2, y="Proportion of Cells", fill=name1)

	#if you want to change colours use:
	# + scale_fill_manual("legend", values = c("CLEC10Ahigh"="dodgerblue1","ICAM1+"="limegreen","ID2+"="deepskyblue","ISG15+"="mediumorchid","LYVE1+"="aquamarine3", "S100A12+"="deeppink","SPP1+"="firebrick1","TREM2high"="darkolivegreen3","TREM2low"="orange2"))
	}

#This is part of VlnPlot.median function
median.stat <- function(x){
	out <- quantile(x, probs = c(0.5))
	names(out) <- c("ymed")
	return(out) 
	}

VlnPlot.median <- function(object, features, colour = "black", pt.size = 0, ncol = NULL, legend = TRUE) {
	#Tris function creates VlnPlot plus a median 
	#VlnPlot.median(object, c("IL1B","CLEC10A"))
	
	myplots <- vector("list")
	
	#Create a plot for each gene and add median
	for (gene in features) {
		#print(gene)
		if (legend == TRUE) {
		myplots[[gene]] <- VlnPlot(object = object, features = gene, pt.size = pt.size) +
		stat_summary(fun = median.stat, geom='point', size = 1, colour = colour)
		} else {
		myplots[[gene]] <- VlnPlot(object = object, features = gene, pt.size = pt.size) +
		stat_summary(fun = median.stat, geom='point', size = 1, colour = colour) + NoLegend()
		}
	}
	#patchwork function to combine multiple plots
	wrap_plots(myplots, ncol=ncol)
	}

Add.SNPs.HT <- function(Seurat.Object, souporcell.file, verbose=FALSE) {
  #Tris function creates will add your souporcell cluster.tsv file as metadata to a Seurat object 
  #PBMC = Add.SNPs.HT(PBMC,"Mapped/clusters.tsv")
  
  SNPs = read.csv(souporcell.file, sep = "\t")
  SNPs
  common_barcode = intersect(colnames(Seurat.Object), SNPs$barcode)
  print(paste("Number of cells barcoded: ",length(common_barcode)))
  
  if (length(common_barcode)/length(colnames(Seurat.Object))<0.5) {
    warning("Less than 50% barcode found")
  }
  rownames(SNPs) <- SNPs$barcode
  SNPs <- SNPs[, c('status', 'assignment')]
  SNPs$SNP_cluster <- SNPs$assignment
  SNPs$assignment <- NULL
  if (verbose){
    head(SNPs)  
  }
  Seurat.Object <- SeuratObject::AddMetaData(Seurat.Object, SNPs)
  return(Seurat.Object)
}

#Generic form
Add.ADT <- function(Seurat.object, ADT.folder.path, verbose=FALSE){
  #Add ADT table to your Seurat object
  #Example use: Seurat.Object =  Add.ADT(Seurat.object, ADT.folder.path="your/path/umi_count/") 
  
  #Read in the ADT library
  ADT.data <- Read10X(ADT.folder.path, gene.column=1)
  #Add "-1"
  colnames(ADT.data) <- paste0(colnames(ADT.data),"-1")
  
  #Just select the common barcodes
  common.barcodes <- intersect(colnames(Seurat.object), colnames(ADT.data))
  print(paste0("Number of cells with ADT: ", length(common.barcodes)))

  no.ADT <- setdiff(colnames(Seurat.object), colnames(ADT.data))
  
  if (verbose==TRUE) {
    print(length(no.ADT))
    print(ncol(Seurat.object))
    print(ncol(ADT.data))
    print(head(ADT.data))
    print(typeof(ADT.data))
  }
  
  #create matrix with 4 columns
  no.ADT.table <- matrix(rep(0, times=nrow(ADT.data)*length(no.ADT)), ncol=length(no.ADT), byrow=TRUE)
  
  #define column names and row names of matrix
  colnames(no.ADT.table) <- no.ADT
  rownames(no.ADT.table) <- rownames(ADT.data)

  if (verbose==TRUE) {
    print(ncol(no.ADT.table))
    print(head(no.ADT.table))
    print(typeof(no.ADT.table))
  }
  
  ADT.final <- cbind(ADT.data, no.ADT.table)
  
  ## ADT preparation
  #For ADT remove unmapped
  ADT.data.filtered = ADT.final[c(1:(length(rownames(ADT.final))-1)) ,]
  #Remove barcode from the name (they are separated with an "-")
  rownames(ADT.data.filtered) = sapply(strsplit(rownames(ADT.data.filtered),"-"), '[', 1)
  
  # #Replace a ADT name with another (in this example with "MHCII")
  # rownames(ADT1.data.b.filtered)
  # if (rownames(ADT1.data.b.filtered)[23]=="HLA") {
  #   rownames(ADT1.data.b.filtered)[23] <- "MHCII"
  # } else {
  #   stop("MHCII is not 23rd position")
  # }
  # rownames(ADT1.data.b.filtered)
  
  print(barplot(rowSums(ADT.data.filtered), main = "ADT total reads", las=2))
  print(head(ADT.data.filtered))
  
  #Add the ADT as an independent assay, before to subset
  Seurat.object[["ADT"]] <- CreateAssayObject(counts = ADT.data.filtered)
  # Normalize ADT data, here we use centered log-ratio (CLR) transformation
  Seurat.object <- NormalizeData(Seurat.object, assay = "ADT", normalization.method = "CLR", margin = 2)
  
  return(Seurat.object)
}

#Generic form
Unique.Name.Of.Function <- function(Seurat.object, var1, varN){
	#Describe your function
	#Example use: Unique.Name.Of.Function(Seurat.object, "var1", "var2") 
	
	#Describe the code
	}









