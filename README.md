# CommonSeuratFunctions
Common Functions useful to work/plot with Seurat objects

## How to use it
Just download the .R file, and load it at the beginning of your script with:

```
source("Path/of/your/file/CommonFunctions.R")
```

and start using the functions

## Functions present

### Proportion Plot
To create

![Proportion Plot example](https://i.stack.imgur.com/4fPY5.png)

use:
```
ProportionPlot(Seurat.Object, "var1", "va2")
```

"var1" and "var2" need to be in your meta.data of course

e.g
```
ProportionPlot(PBMC, "seurat_cluster", "Condition")
```

---

### Violin Plot + Median
To create

![VlnPlot and Median](https://user-images.githubusercontent.com/34346930/161740981-db3f078d-df2b-4841-a669-ade8b56c821e.png)

use:
```
VlnPlot.median(Seurat.Object, features, colour = "black", pt.size = 0, ncol = NULL, legend = TRUE)
```

colour = "black"    Is the color of the median dot

pt.size = 0         size of the cell dots

ncol = NULL         n of columns to use for multiple genes plot

legend = TRUE       whether or not show a legend

---

### Add souporcell clusters.tsv metadata to Seurat

Souporcell generates a clusters.tsv with the SNPs cluster for each sample.

use:
```
Seurat.Object = Add.SNPs.HT(Seurat.Object,"Path/of/your/file/clusters.tsv")
```

---

### Calculate and Add ADT (CITEseq) to Seurat

CITE-seq-Count tool map the ADT reads and generate a folder named "umi_count" with the data.
You can also replace some ADT names with others 

use:
```
Seurat.Object = Add.ADT(Seurat.Object,"Path/of/your/file/umi_count/", replace.any=c("OLD_NAME"="NEW_NAME"))
```

---

### Calculate and Add HTO (CITEseq) to Seurat (using cellhashR)

Use CITE-seq-Count tool to map the HTO reads and generate a folder named "umi_count" with the data. 
This function use cellhashR to deconvolute the hashtag information and add it to your seurat object

use:
```
Seurat.Object = Add.HTO(Seurat.Object,"Path/of/your/file/umi_count/")
```

---

### Use SoupX to clean cellranger matrix from Ambient RNA contaminat

You need "filtered_feature_bc_matrix", "raw_feature_bc_matrix", and "analysis" folders to use SoupX. You use this instead of Read10X function.

use:
```
clean.data = SoupX.clean(cellranger.folder="/my/cellranger/folder/")
```

---

### Create a Seurat.object cleaned with SoupX, starting from a Seurat.object

Starting from a Seurat.object, it will clean it (using SoupX) and generates another Seurat.object.

You need "filtered_feature_bc_matrix", "raw_feature_bc_matrix" folders (not "analysis")  from Cellranger. 

Seurat object need to have "seurat_clusters", and need to be UNFILTERED!! (don't remove doublets, high MT cells, and so on).

```
Seurat.object.clean <- SoupX.on.Seurat(Seurat.object, cellranger.folder)
```

---

### Add metadata to your seurat object

Starting from a seurat object and a metadata table, 

this function creates a proper metadata table with cell barcodeID and add it to Seurat.object.

Of course table rownames MUST match Seurat.object Idents. 

```
Seurat.object <- make.add.meta(Seurat.object, Seurat.object)
```

---

### Read BD Rhapsody folder

Starting from a BD Rhapsody mapping folder, this function read it and create a seurat object

At the moment it works using DBEC MolPerCell, to implement RSEC

use:

```
Read.BD.Rhap(path, project.name="BD Rhapsody", remove.multiplets = TRUE, remove.undetermined=TRUE,  min.cells = 3, min.features = 20, show.plots=TRUE, verbose=FALSE)
```

path 							folder path

project.name="BD Rhapsody"		Seurat object folder name

remove.multiplets = TRUE		whether remove or not multiplets

remove.undetermined = TRUE		whether remove or not undetermined

min.cells = 3					min number of cells to keep a feature

min.features = 20 				min number of features to keep a cell

show.plots = TRUE				whether show or not some readcounts plot
---

### Read BD Rhapsody file

Starting from a BD Rhapsody MolPerCell folder, this function read it and create a seurat object

This function does NOT read the SampleTag file, so it's not able to associate each cell with a sample tag

```
Read.BD.Rhap.simple <- function(MolsPerCell.file, project.name="BD Rhapsody", min.cells = 3, min.features = 20, show.plots=TRUE, verbose=FALSE)
```

MolsPerCell.file 				file path

project.name="BD Rhapsody"		Seurat object folder name

min.cells = 3					min number of cells to keep a feature

min.features = 20 				min number of features to keep a cell

show.plots = TRUE				whether show or not some readcounts plot

---

### Plot Depth of sequencing

Starting from a seurat object, plot the depth  of sequencing (splitted for a specific metadata column) 

```
plot.depth.seq(Seurat.object, metadata.col)
```