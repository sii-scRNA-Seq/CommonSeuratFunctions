# CommonSeuratFunctions
Common Functions useful to work/plot with Seurat objects

## How to use it
Just download the .R file, and load it at the beginning of your script with:

```
source("Path/of/your/file/CommonFunction.R")
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

### Add souporcell clusters.tsv metadata to Seurat

Souporcell generates a clusters.tsv with the SNPs cluster for each sample.

use:
```
Seurat.Object = Add.SNPs.HT(Seurat.Object,"Path/of/your/file/clusters.tsv")
```

### Calculate and Add ADT (CITEseq) to Seurat

CITE-seq-Count tool map the ADT reads and generate a folder named "umi_count" with the data.
You can also replace some ADT names with others 

use:
```
Seurat.Object = Add.ADT(Seurat.Object,"Path/of/your/file/umi_count/", replace.any=c("OLD_NAME"="NEW_NAME"))
```

### Calculate and Add HTO (CITEseq) to Seurat (using cellhashR)

Use CITE-seq-Count tool to map the HTO reads and generate a folder named "umi_count" with the data. 
This function use cellhashR to deconvolute the hashtag information and add it to your seurat object

use:
```
Seurat.Object = Add.HTO(Seurat.Object,"Path/of/your/file/umi_count/")
```

### Use SoupX to clean cellranger matrix from Ambient RNA contaminat

You need "filtered_feature_bc_matrix", "raw_feature_bc_matrix", and "analysis" folders to use SoupX. You use this instead of Read10X function.

use:
```
clean.data = SoupX.clean(cellranger.folder="/my/cellranger/folder/")
```

### Create a Seurat.object cleaned with SoupX, starting from a Seurat.object

Starting from a Seurat.object, it will clean it (using SoupX) and generates another Seurat.object.

You need "filtered_feature_bc_matrix", "raw_feature_bc_matrix" folders (not "analysis")  from Cellranger. 

Seurat object need to have "seurat_clusters", and need to be UNFILTERED!! (don't remove doublets, high MT cells, and so on).

```
Seurat.object.clean <- SoupX.on.Seurat(Seurat.object, cellranger.folder)
```