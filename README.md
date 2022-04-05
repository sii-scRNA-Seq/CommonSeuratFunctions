# CommonSeuratFunctions
Common Functions useful to work/plot with Seurat objects

## How to use it
Just download the .R file in one folder, and load it at the beginning of your script with:
source("Path/of/your/file/CommonFunction.R")

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
VlnPlot.median(object, features, colour = "black", pt.size = 0, ncol = NULL, legend = TRUE)
```


colour = "black"    Is the color of the median dot

pt.size = 0         size of the cell dots

ncol = NULL         n of columns to use for multiple genes plot

legend = TRUE       whether or not show a legend
