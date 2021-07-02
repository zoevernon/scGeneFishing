# GeneFishing
## Description
This is the GitHub repository for the GeneFishing R package.  [GeneFishing] (https://www.pnas.org/content/116/38/18943) (Liu et al., 2019) is a method identify genes involved in a biological using a set of genes previously known to be a part of the process.  The original paper looked at bulk RNA-seq data and found novel genes involved in the cholesterol biosynthesis pathway.   In this package we updated the method, allowing the user to define a large set of genes that can potenially be used as bait and utilize UMAP to improve method performance on scRNA-seq data.  

To do GeneFishing you need two pieces of data: 
1. A count matrix (appropriate cleaned and normalized, if appropriate) with genes as rows and cells or samples as columns. 
2. A set of bait related to the biological process being studied, where the genes in that set are highly coexpressed in the count matrix. 

In this we go over how one might select a bait set from a larger set of genes, such as a GO term or pathway, and subsequently perform GeneFishing.

## Installing
To install the package use the following comands.  Note, you should uncomment the first line ```# install.packages("devtools")``` if you have not installed the ```devtools``` package.  
```{r, eval = F}
# install.packages("devtools")
devtools::install_github("zoevernon/GeneFishing")
```

## Step 1: define a bait set 
The ```probeFishability``` function will automatically search through a set of potential bait genes to determine if there are any subsets which can be used in ```GeneFishing```.  
```{r, eval = F}
discovered_bait <- probeFishability(X, 
                                    potential_bait, 
                                    n_rounds = 100, 
                                    min_tightness = 0.5, 
                                    alpha = 5,
                                    umap = TRUE,
                                    ncores = 2)
```
The tuning parameters ```min_tightness``` and ```alpha``` control the level of coexpression which is required to say a set of genes can be used as bait.  

## Step 2: run ```geneFishing```
If the user already has a bait set they want to use, or has used ```probeFishability``` to find one, they can then run ```geneFishing``` to search for additional genes that are likely to related to the function of the bait set. 

```{r, eval = F}
discovered_bait <- geneFishing(X, 
                               bait_genes, 
                               k = 2, 
                               umap = TRUE,
                               alpha = 5,
                               n_rounds = 1000,
                               min_tightness = 0.5,
                               ncores = 2)
```

In this function ```k```  is the number of clusters to use in ```geneFishing``` and when ```umap = TRUE``` the algorithm will use UMAP coordinates to do the clustering.  Before doing GeneFishing it will first check if the supplied bait is tight enough to use, if not it will suggest the user try ```probeFishability```.  Otherwise, it will iterate through ```n_rounds``` of GeneFishing and returned a data.frame with the capture frequency rate (CFR) of all the genes in the data.  


