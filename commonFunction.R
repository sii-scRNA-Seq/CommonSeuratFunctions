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


#Generic form
Unique.Name.Of.Function <- function(Seurat.object, var1, varN){
	#Describe your function
	#Example use: Unique.Name.Of.Function(Seurat.object, "var1", "var2") 
	
	#Describe the code
	}
