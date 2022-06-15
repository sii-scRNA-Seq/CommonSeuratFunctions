#This is a collection of functions useful for Seurat objects.
#To load it run:
#source("commonFunctions.R")

ProportionPlot <- function(object, name1, name2, cols=NULL){
	#Tris function creates Proportion plot in % for two variables in the seurat object 
	#Example use: ProportionPlot(PBMC, "cluster", "Celltype") 
  ##if you want to change colours use:
  # ...cols = c("CLEC10Ahigh"="dodgerblue1","ICAM1+"="limegreen","ID2+"="deepskyblue","ISG15+"="mediumorchid","LYVE1+"="aquamarine3", "S100A12+"="deeppink","SPP1+"="firebrick1","TREM2high"="darkolivegreen3","TREM2low"="orange2"))
	
	#Get number of cells per cluster
	slot.metadata = slot(object, "meta.data")
	clust.table.group <- table(getElement(slot.metadata,name1), getElement(slot.metadata,name2))

	#Convert to proportion/percentage of total cells in each sample
	clust.prop.group <- as.data.frame(prop.table(x = clust.table.group, margin = 2))
  
	if (is.null(cols)==TRUE) {
  	#Plot as stacked bar plot
  	ggplot(data=clust.prop.group, aes(x=Var2, y=Freq, fill=Var1)) + theme_linedraw() +
  	geom_bar(stat="identity", color="black") + labs(x=name2, y="Proportion of Cells", fill=name1)
	} else {
	  ggplot(data=clust.prop.group, aes(x=Var2, y=Freq, fill=Var1)) + theme_linedraw() +
	    geom_bar(stat="identity", color="black") + labs(x=name2, y="Proportion of Cells", fill=name1) + 
	    scale_fill_manual("legend", values = cols)
	  }
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
		myplots[[gene]] <- Seurat::VlnPlot(object = object, features = gene, pt.size = pt.size) +
		stat_summary(fun = median.stat, geom='point', size = 1, colour = colour)
		} else {
		myplots[[gene]] <- Seurat::VlnPlot(object = object, features = gene, pt.size = pt.size) +
		stat_summary(fun = median.stat, geom='point', size = 1, colour = colour) + NoLegend()
		}
	}
	#patchwork function to combine multiple plots
	patchwork::wrap_plots(myplots, ncol=ncol)
	}

Add.SNPs.HT <- function(Seurat.Object, souporcell.file, verbose=FALSE) {
  #Tris function creates will add your souporcell cluster.tsv file as metadata to a Seurat object 
  #PBMC = Add.SNPs.HT(PBMC,"Mapped/clusters.tsv")
  
  SNPs = read.csv(souporcell.file, sep = "\t")
  common_barcode = intersect(colnames(Seurat.Object), SNPs$barcode)
  print(paste("Number of cells barcoded: ",length(common_barcode)))
  
  if (length(common_barcode)/length(colnames(Seurat.Object))<0.5) {
    warning("Less than 50% barcode found")
  }
  rownames(SNPs) <- SNPs$barcode
  SNPs <- SNPs[, c('status', 'assignment')]
  SNPs$SNP_cluster <- SNPs$assignment
  SNPs$SNP_status <- SNPs$status
  SNPs$assignment <- NULL
  SNPs$status <- NULL
  if (verbose){
    head(SNPs)  
  }
  Seurat.Object <- SeuratObject::AddMetaData(Seurat.Object, SNPs)
  return(Seurat.Object)
}

#Generic form
Add.ADT <- function(Seurat.Object, ADT.folder.path, replace.any=FALSE, verbose=FALSE){
  #Add ADT table to your Seurat object
  #Example use: Seurat.Object =  Add.ADT(Seurat.Object, ADT.folder.path="your/path/umi_count/",replace.any=c("HLA"="MHCII")) 
  
  #Read in the ADT library
  ADT.data <- Seurat::Read10X(ADT.folder.path, gene.column=1)
  #Add "-1"
  colnames(ADT.data) <- paste0(colnames(ADT.data),"-1")
  
  #Just select the common barcodes
  common.barcodes <- intersect(colnames(Seurat.Object), colnames(ADT.data))
  print(paste0("Number of cells in Seurat with ADT: ", length(common.barcodes)))
  
  if (length(common.barcodes)==0) {
    stop("Stop! No cells have ADT")
  }

  no.ADT <- setdiff(colnames(Seurat.Object), colnames(ADT.data))
  extra.ADT <- setdiff(colnames(ADT.data), colnames(Seurat.Object))
  
  if (length(extra.ADT)!=0) {
    ADT.data <- ADT.data[, !colnames(ADT.data) %in% extra.ADT]
  }
  
  if (verbose==TRUE) {
    print(paste0("Number of cells in Seurat: ",ncol(Seurat.Object)))
    print(paste0("Number of cells in ADT: ",ncol(ADT.data)))
    print(head(ADT.data))
    print(paste0("Number of cells in without any ADT: ",length(no.ADT)))
  }
  
  if (length(no.ADT)!=0) {
  
    #create matrix with 4 columns
    no.ADT.table <- matrix(rep(0, times=nrow(ADT.data)*length(no.ADT)), ncol=length(no.ADT), byrow=TRUE)
    
    if (verbose==TRUE) {
      print(head(no.ADT.table))
    }
    
    #define column names and row names of matrix
    colnames(no.ADT.table) <- no.ADT
    rownames(no.ADT.table) <- rownames(ADT.data)
  
    if (verbose==TRUE) {
      print(ncol(no.ADT.table))
      print(head(no.ADT.table))
      print(typeof(no.ADT.table))
    }
    
    ADT.final <- cbind(ADT.data, no.ADT.table)
    
    if (verbose==TRUE) {
      print(ncol(ADT.final))
    }
  } else {
    ADT.final=ADT.data
  }
  
  ## ADT preparation
  #For ADT remove unmapped
  ADT.data.filtered = ADT.final[c(1:(length(rownames(ADT.final))-1)) ,]
  #Remove barcode from the name (they are separated with an "-")
  rownames(ADT.data.filtered) = sapply(strsplit(rownames(ADT.data.filtered),"-"), '[', 1)
  
  #Replace a ADT name with another
  if (replace.any!=FALSE) {
    print("Replacing ADT names...")
    for (i in 1:length(replace.any))
    if (any(rownames(ADT.data.filtered)==names(replace.any[i]))) {
      pos = grep(names(replace.any[i]), rownames(ADT.data.filtered))
      rownames(ADT.data.filtered)[pos] <- replace.any[[i]]
    }
  }
  
  print(barplot(rowSums(ADT.data.filtered), main = "ADT total reads", las=2))
  print(head(ADT.data.filtered))
  
  #Add the ADT as an independent assay, before to subset
  Seurat.Object[["ADT"]] <- SeuratObject::CreateAssayObject(counts = ADT.data.filtered)
  # Normalize ADT data, here we use centered log-ratio (CLR) transformation
  Seurat.Object <- Seurat::NormalizeData(Seurat.Object, assay = "ADT", normalization.method = "CLR", margin = 2)
  
  return(Seurat.Object)
}

extract.HTO <- function(path, barcodeWhitelist, minCountPerCell = 5, methods = c("bff_cluster", "multiseq","dropletutils")) {
	#Use:
	#> my.HTO.table <- extract.HTO("/mypath/", c("HTO1","HTO6")))
  
  print("Process Count Matrix...")
	barcodeData <- cellhashR::ProcessCountMatrix(rawCountData = path, minCountPerCell = minCountPerCell, barcodeWhitelist = barcodeWhitelist)
	print("Generate Cell Hashing Calls...")
	calls.HTO <- cellhashR::GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = methods)

	HTOtable <- data.frame(row.names = calls.HTO$cellbarcode)
	HTOtable$HTO <- calls.HTO$consensuscall
	HTOtable$HTO_status <- calls.HTO$consensuscall.global
	rownames(HTOtable) <- paste0(rownames(HTOtable),"-1")

	# Inspect negative cells:
	print("Summarize Cells By Classification...")
	cellhashR::SummarizeCellsByClassification(calls = calls.HTO, barcodeMatrix = barcodeData)
	return(HTOtable)
}

Add.HTO <- function(Seurat.Object, path, barcodeWhitelist, minCountPerCell = 5, methods = c("bff_cluster", "multiseq","dropletutils")) {
	#Add HTO table to your Seurat object
	#Use:
	#> Seurat.Object <- (Seurat.Object, "/mypath/HTO_folder/", c("HTO1","HTO6")))
	HTOtable <- extract.HTO(path, barcodeWhitelist, minCountPerCell = minCountPerCell, methods = methods)
	Seurat.Object <- SeuratObject::AddMetaData(Seurat.Object, HTOtable)
	return(Seurat.Object)
}

SoupX.clean.from.CellRanger <- function(cellranger.folder) {
  #Clean your scRNAseq data from ambient contamination
  #Use:
  #clean.data = SoupX.clean(cellranger.folder="/my/cellranger/folder/")
  
  if (dir.exists(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))==FALSE){
    stop("filtered_feature_bc_matrix folder didn't find")
  }
  
  if (dir.exists(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))==FALSE){
    stop("raw_feature_bc_matrix folder didn't find")
  }
  
  if (dir.exists(paste(cellranger.folder,"analysis",sep = "/"))==FALSE){
    stop("analysis folder didn't find")
  }
  
  #Load data
  sc = SoupX::load10X(cellranger.folder)
  #Estimante the contamination
  sc = SoupX::autoEstCont(sc)
  
  print("Genes with highest expression in background:")
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  
  #Filterout the contamination
  out = SoupX::adjustCounts(sc)
  return(out)
}

SoupX.on.Seurat <- function(Seurat.object, cellranger.folder, min.cells = 5, min.features = 50) {
  #Clean your scRNAseq data from ambient contamination
  #Seurat object should be UNFILTERED (MT genes and doublets), and seurat_clusters should be present
  #Use:
  #clean.data = SoupX.clean(cellranger.folder="/my/cellranger/folder/")
  
  if (dir.exists(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))==FALSE){
    stop("filtered_feature_bc_matrix folder didn't find")
  }
  
  if (dir.exists(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))==FALSE){
    stop("raw_feature_bc_matrix folder didn't find")
  }
  
  if (any(colnames(Seurat.object@meta.data) == "seurat_clusters")) {
    print("seurat_clusters present...")
  } else {
    stop("seurat_clusters not present")
  }
  
  raw.matrix = Seurat::Read10X(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))
  filt.matrix = Seurat::Read10X(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))
  
  sc  <- SoupX::SoupChannel(raw.matrix, filt.matrix)
  
  meta    <- Seurat.object@meta.data
  umap    <- Seurat.object@reductions$umap@cell.embeddings
  sc  <- SoupX::setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- SoupX::setDR(sc, umap)
  head(meta)
  
  #With defined clusters, run the main SoupX function, calculating ambient RNA profile.
  sc  <- SoupX::autoEstCont(sc)
  
  print("Genes with highest expression in background:")
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  
  adj.matrix = SoupX::adjustCounts(sc, roundToInt = T)
  
  New.Seurat.object = Seurat::CreateSeuratObject(counts = adj.matrix, min.cells = 5, min.features = 50, meta.data = meta)
  return(New.Seurat.object)
}

make.add.meta <- function(Seurat.Object, metadata) {
  #From a metadata table and a Seurat.Object, 
  #it creates a proper metadata table with cell barcodeID and add it to Seurat.object 
  # Use:
  # Seurat.object = make.add.meta(Seurat.Object, metadata)
  
  #if there is only 1 row, add it to the whole Seurat.Object
  if (nrow(metadata)==1) {
    df.cells <- data.frame(row.names = colnames(Seurat.Object))
    for (name in colnames(metadata)) {
      #print(name)
      df.cells[name]=metadata[name]
    }
  } else {
    #check if length of idents is different of lenght metadata
    if (length(setdiff(Idents(Seurat.Object), rownames(metadata)))!=0 || length(setdiff(rownames(metadata),Idents(Seurat.Object)))!=0) {
      stop("Seurat object Idents and metadata rows are not matching.")
    }
    
    df.cells <- data.frame()
    
    #Select the Idents matching “Name” in metadata
    #Idents(Seurat.Object) <- "Condition_name"
    
    #For each idents within the seurat object (cluster or sample you want)
    for (name in unique(Idents(Seurat.Object)))
    {
      print(name)
      
      #Select the cells and put them in another df
      new_df <- data.frame(row.names = WhichCells(Seurat.Object, idents = name))
      
      #Select the row of interest from metadata, corresponding to the metadata to add to those cells
      meta_row=metadata[rownames(metadata) == name,]
      for (col_name in colnames(meta_row)){
        #add the specific metadata you need
        new_df[col_name]=meta_row[col_name]
      }
      #merge in a big df
      df.cells=dplyr::bind_rows(df.cells, new_df)
    }
  }
  
  #Finally add to the object
  Seurat.Object <- AddMetaData(Seurat.Object, metadata = df.cells)
  return(Seurat.Object)
}

Read.BD.Rhap <- function(path, project.name="BD Rhapsody", remove.multiplets = TRUE, remove.undetermined=TRUE,  min.cells = 3, min.features = 20, show.plots=TRUE, verbose=FALSE) {
  # path = folder containing DBEC_MolsPerCell and SampleTag files from 7B
  # TO DO: implement RSEC
  # See https://scomix.bd.com/hc/en-us/articles/360044971032-Bioinformatics
  
  fileList <- list.files(path=path)
  
  #Keep only files of interest (i.e. those which contain UMInormalized counts and are sample tagged)
  print("Reading files...")
  DBEC.file <- fileList[grepl("DBEC_MolsPerCell", fileList)]
  #RSEC.file <- fileList[grepl("RSEC_MolsPerCell", fileList)]  
  sampleTag.file<- fileList[grepl("Sample_Tag", fileList)]
  
  count <- read.csv(file = paste(path, DBEC.file, sep="/"), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, comment.char = "#")
  
  #traspose gene/cells
  count <- t(count)
  if (verbose) {
    print("Counts rownames:")
    print(head(rownames(count)))
    print("Counts colnames:")
    print(head(colnames(count)))
  }
  
  #Read sampleTag
  ST <- read.csv(file = paste(path, sampleTag.file, sep="/"), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, comment.char = "#")
  if (verbose) {
    print(sampleTag.file)
    print("SampleTag:")
    print(head(ST))
  }
  
  #Create Seurat
  Seurat.object <- CreateSeuratObject(counts = count, min.cells = min.cells, min.features = min.features, project = project.name)
  Seurat.object = AddMetaData(Seurat.object, ST)
  Seurat.object@meta.data$file <- sampleTag.file
  print(Seurat.object)
  
  Idents(Seurat.object) <- "Sample_Name"
  if (show.plots==TRUE){
    print(qplot(Seurat.object$Sample_Name, main=project.name))
  }
  
  if (remove.multiplets==TRUE){
    print("Removing Multiplets...")
    Seurat.object = subset(Seurat.object, idents = "Multiplet", invert = TRUE)
  }
  if (remove.multiplets==TRUE){
    print("Removing Undetermined")
    Seurat.object = subset(Seurat.object, idents ="Undetermined", invert = TRUE)
  }
  if (show.plots==TRUE){
    hist(colSums(Seurat.object@assays$RNA@counts), main=project.name, breaks = 50)
  }
  return(Seurat.object)
}

Read.BD.Rhap.simple <- function(MolsPerCell.file, project.name="BD Rhapsody", min.cells = 3, min.features = 20, show.plots=TRUE, verbose=FALSE) {
  # MolsPerCell.file = file containing MolsPerCell and SampleTag files from 7B
  
  #Keep only files of interest (i.e. those which contain UMInormalized counts and are sample tagged)
  print("Reading file...")
  
  count <- read.csv(file = MolsPerCell.file, sep = ',', header = TRUE, row.names = 1, check.names = FALSE, comment.char = "#")
  
  #traspose gene/cells
  count <- t(count)
  if (verbose) {
    print("Counts rownames:")
    print(head(rownames(count)))
    print("Counts colnames:")
    print(head(colnames(count)))
  }
  
  #Create Seurat
  Seurat.object <- CreateSeuratObject(counts = count, min.cells = min.cells, min.features = min.features, project = project.name)
  print(Seurat.object)

  if (show.plots==TRUE){
    hist(colSums(Seurat.object@assays$RNA@counts), main=project.name, breaks = 50)
  }
  return(Seurat.object)
}

plot.depth.seq <- function(Seurat.object, metadata.col) {
  #metadata.col need to be a metadata column in Seurat.object@meta.data
  
  Idents(Seurat.object) <- metadata.col
  libs = unique(Idents(Seurat.object))
  
  for (i in 1:length(libs)) {
    lib.seurat = subset(Seurat.object, idents = libs[i])

    colSums.counts = as.data.frame(colSums(lib.seurat@assays$RNA@counts))
    colSums.counts$colSums <- colSums.counts$`colSums(lib.seurat@assays$RNA@counts)`
    colSums.counts$`colSums(lib.seurat@assays$RNA@counts)` <- NULL
    #print(colnames(colSums.counts))
    a=ggplot(colSums.counts, aes(x=colSums))+ geom_histogram(color="black", fill="white", bins=50) + ggtitle(libs[[i]])
    b=ggplot(colSums.counts, aes(x=colSums))+ geom_histogram(color="black", fill="white", bins=50, aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666") + ggtitle(libs[[i]])
    print(a+b)
    #print(ggplot(colSums.counts, aes(x=colSums)) + geom_histogram(color="black", fill="white", bins=50, aes(y=..density..)) + geom_density(alpha=.2, fill="#FF6666") + ggtitle(unique(Idents(libs[[i]]))))
  }
}

#Generic form
Unique.Name.Of.Function <- function(Seurat.Object, var1, varN){
	#Describe your function
	#Example use: Unique.Name.Of.Function(Seurat.Object, "var1", "var2") 
  
	#Describe the code
  #Code
	}









