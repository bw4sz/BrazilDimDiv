Using data.table to read in raw betadiversity output
========================================================

The package data.table provides huge opportunities to use sql like syntax to make really efficient queries of large tables.

I'll show an example for how to compute mean and sd on each cell, but this is just the **unique** rows and not the final data. I want to pass you the raw data so that you are not constrained by the types of data output i compute.

We will all talk about it more but Gabriel and I initially brainstormed a few metrics
  1. Avg betadiversity for taxonomic, phylogenetic, trait for each cell
  2. Variance in betadiversity for each cell
  3. Quantile of the mean and standard deviation, this would entail:
    * Finding steps 1 and 2
    * Finding the quantile of the mean of each cell to the distribution of means for all cells
    * Same for variance,
    *This would be helpful for deliniating 'high' or 'low' betadiversity quantiles
  4. Compare combinations of betadiversity to env betadiversity
  5. Comparing combinations of betadiversity to distance
  6. Overlaying cells on a biome map and looking at within and between
  
There is alot more you can do with the raw data, which is why i wanted to show the power of the data.table packages

All of this data is available on the github repo https://github.com/bw4sz/BrazilDimDiv .


```r
setwd("C:/Users/Jorge/Documents/BrazilDimDiv/cluster")
require(data.table)
```

```
## Loading required package: data.table
```

```r
require(reshape)
```

```
## Loading required package: reshape
```

```r
# use fread to read in files, this is much much faster than read.csv
dat <- fread("FinalData.csv")
head(dat)
```

```
##    combo To.OriginalRow To.xcord To.ycord To.x From.OriginalRow From.xcord
## 1:  8_12            727  9913718  7594048   12            12052   12213718
## 2:  8_15           3097 -4786282  5734048   15            12052   12213718
## 3:  8_17          11943 10613718 -3937952   17            12052   12213718
## 4:  8_17          12001 10513718 -4061952   17            12052   12213718
## 5:  8_24           2932 -4186282  5858048   24            12052   12213718
## 6:  8_25           3439 -4286282  5486048   25            12052   12213718
##    From.ycord From.x To.y From.y BetaSim Sorenson  MNTD
## 1:   -4185952      8    8     12  0.4561      1.0 5.292
## 2:   -4185952      8    8     15  0.4561      1.0 4.365
## 3:   -4185952      8    8     17  0.0000      0.5 2.796
## 4:   -4185952      8    8     17  0.0000      0.5 2.796
## 5:   -4185952      8    8     24  0.4561      1.0 5.209
## 6:   -4185952      8    8     25  0.4561      1.0 4.640
```


There is alot of duplicitious info here, let's clean off some columns and just get what we want.


```r
dat[, `:=`(c("combo", "To.x", "From.x", "From.y", "To.y"), NULL)]
```

```
##        To.OriginalRow To.xcord To.ycord From.OriginalRow From.xcord
##     1:            727  9913718  7594048            12052   12213718
##     2:           3097 -4786282  5734048            12052   12213718
##     3:          11943 10613718 -3937952            12052   12213718
##     4:          12001 10513718 -4061952            12052   12213718
##     5:           2932 -4186282  5858048            12052   12213718
##    ---                                                             
## 15881:           3270 -4186282  5610048             3440   -4186282
## 15882:           2760 -4186282  5982048             3267   -4486282
## 15883:           2931 -4286282  5858048             3267   -4486282
## 15884:           3102 -4186282  5734048             3267   -4486282
## 15885:           3270 -4186282  5610048             3267   -4486282
##        From.ycord BetaSim Sorenson    MNTD
##     1:   -4185952 0.45611  1.00000 5.29150
##     2:   -4185952 0.45611  1.00000 4.36480
##     3:   -4185952 0.00000  0.50000 2.79591
##     4:   -4185952 0.00000  0.50000 2.79591
##     5:   -4185952 0.45611  1.00000 5.20910
##    ---                                    
## 15881:    5486048 0.00000  0.06667 0.02505
## 15882:    5610048 0.01855  0.12500 0.14251
## 15883:    5610048 0.01855  0.12500 0.14251
## 15884:    5610048 0.01855  0.12500 0.14251
## 15885:    5610048 0.01855  0.12500 0.14251
```

```r
head(dat)
```

```
##    To.OriginalRow To.xcord To.ycord From.OriginalRow From.xcord From.ycord
## 1:            727  9913718  7594048            12052   12213718   -4185952
## 2:           3097 -4786282  5734048            12052   12213718   -4185952
## 3:          11943 10613718 -3937952            12052   12213718   -4185952
## 4:          12001 10513718 -4061952            12052   12213718   -4185952
## 5:           2932 -4186282  5858048            12052   12213718   -4185952
## 6:           3439 -4286282  5486048            12052   12213718   -4185952
##    BetaSim Sorenson  MNTD
## 1:  0.4561      1.0 5.292
## 2:  0.4561      1.0 4.365
## 3:  0.0000      0.5 2.796
## 4:  0.0000      0.5 2.796
## 5:  0.4561      1.0 5.209
## 6:  0.4561      1.0 4.640
```


That's a bit better, deleting columns in data.table is the same as setting them to NULL

We can look at the stucture of the data.table

```r
str(dat)
```

```
## Classes 'data.table' and 'data.frame':	15885 obs. of  9 variables:
##  $ To.OriginalRow  : int  727 3097 11943 12001 2932 3439 9632 8684 8378 10297 ...
##  $ To.xcord        : num  9913718 -4786282 10613718 10513718 -4186282 ...
##  $ To.ycord        : num  7594048 5734048 -3937952 -4061952 5858048 ...
##  $ From.OriginalRow: int  12052 12052 12052 12052 12052 12052 12052 12052 12052 12052 ...
##  $ From.xcord      : num  12213718 12213718 12213718 12213718 12213718 ...
##  $ From.ycord      : num  -4185952 -4185952 -4185952 -4185952 -4185952 ...
##  $ BetaSim         : num  0.456 0.456 0 0 0.456 ...
##  $ Sorenson        : num  1 1 0.5 0.5 1 1 1 1 1 1 ...
##  $ MNTD            : num  5.29 4.36 2.8 2.8 5.21 ...
##  - attr(*, ".internal.selfref")=<externalptr>
```


Note that the Original.Row columns are integers


We can set columns in a similiar format to deleting them, here let's make the original row columns in a character so 1 becomes "1"


```r
dat[, `:=`(To.OriginalRow, as.character(To.OriginalRow))]
```

```
##        To.OriginalRow To.xcord To.ycord From.OriginalRow From.xcord
##     1:            727  9913718  7594048            12052   12213718
##     2:           3097 -4786282  5734048            12052   12213718
##     3:          11943 10613718 -3937952            12052   12213718
##     4:          12001 10513718 -4061952            12052   12213718
##     5:           2932 -4186282  5858048            12052   12213718
##    ---                                                             
## 15881:           3270 -4186282  5610048             3440   -4186282
## 15882:           2760 -4186282  5982048             3267   -4486282
## 15883:           2931 -4286282  5858048             3267   -4486282
## 15884:           3102 -4186282  5734048             3267   -4486282
## 15885:           3270 -4186282  5610048             3267   -4486282
##        From.ycord BetaSim Sorenson    MNTD
##     1:   -4185952 0.45611  1.00000 5.29150
##     2:   -4185952 0.45611  1.00000 4.36480
##     3:   -4185952 0.00000  0.50000 2.79591
##     4:   -4185952 0.00000  0.50000 2.79591
##     5:   -4185952 0.45611  1.00000 5.20910
##    ---                                    
## 15881:    5486048 0.00000  0.06667 0.02505
## 15882:    5610048 0.01855  0.12500 0.14251
## 15883:    5610048 0.01855  0.12500 0.14251
## 15884:    5610048 0.01855  0.12500 0.14251
## 15885:    5610048 0.01855  0.12500 0.14251
```

```r

dat[, `:=`(From.OriginalRow, as.character(From.OriginalRow))]
```

```
##        To.OriginalRow To.xcord To.ycord From.OriginalRow From.xcord
##     1:            727  9913718  7594048            12052   12213718
##     2:           3097 -4786282  5734048            12052   12213718
##     3:          11943 10613718 -3937952            12052   12213718
##     4:          12001 10513718 -4061952            12052   12213718
##     5:           2932 -4186282  5858048            12052   12213718
##    ---                                                             
## 15881:           3270 -4186282  5610048             3440   -4186282
## 15882:           2760 -4186282  5982048             3267   -4486282
## 15883:           2931 -4286282  5858048             3267   -4486282
## 15884:           3102 -4186282  5734048             3267   -4486282
## 15885:           3270 -4186282  5610048             3267   -4486282
##        From.ycord BetaSim Sorenson    MNTD
##     1:   -4185952 0.45611  1.00000 5.29150
##     2:   -4185952 0.45611  1.00000 4.36480
##     3:   -4185952 0.00000  0.50000 2.79591
##     4:   -4185952 0.00000  0.50000 2.79591
##     5:   -4185952 0.45611  1.00000 5.20910
##    ---                                    
## 15881:    5486048 0.00000  0.06667 0.02505
## 15882:    5610048 0.01855  0.12500 0.14251
## 15883:    5610048 0.01855  0.12500 0.14251
## 15884:    5610048 0.01855  0.12500 0.14251
## 15885:    5610048 0.01855  0.12500 0.14251
```

```r

# look at the str again
str(dat)
```

```
## Classes 'data.table' and 'data.frame':	15885 obs. of  9 variables:
##  $ To.OriginalRow  : chr  "727" "3097" "11943" "12001" ...
##  $ To.xcord        : num  9913718 -4786282 10613718 10513718 -4186282 ...
##  $ To.ycord        : num  7594048 5734048 -3937952 -4061952 5858048 ...
##  $ From.OriginalRow: chr  "12052" "12052" "12052" "12052" ...
##  $ From.xcord      : num  12213718 12213718 12213718 12213718 12213718 ...
##  $ From.ycord      : num  -4185952 -4185952 -4185952 -4185952 -4185952 ...
##  $ BetaSim         : num  0.456 0.456 0 0 0.456 ...
##  $ Sorenson        : num  1 1 0.5 0.5 1 1 1 1 1 1 ...
##  $ MNTD            : num  5.29 4.36 2.8 2.8 5.21 ...
##  - attr(*, ".internal.selfref")=<externalptr>
```


Set key
---

Data.tables have keys which sort the objects

```r
setkey(dat, To.OriginalRow)
```



Susbetting
-------

Data.table subsets tables very very fast. Think of the code inside of the [] as a logical query of the key NOT as a index! This is the important, and difficult to adapt to distinction.



```r
# Subsetting is faster when done directly on the key column
system.time(row1 <- dat["1"])
```

```
##    user  system elapsed 
##    0.02    0.00    0.00
```

```r

# Asking the column directly is slightly slower (still fast)
system.time(row1 <- dat[To.OriginalRow %in% "1"])
```

```
##    user  system elapsed 
##       0       0       0
```

```r

# Asking for two columns, since cell 8 could be in To or From
system.time(row1 <- dat[To.OriginalRow %in% "1" | From.OriginalRow %in% "1"])
```

```
##    user  system elapsed 
##       0       0       0
```

```r
head(row1)
```

```
##    To.OriginalRow To.xcord To.ycord From.OriginalRow From.xcord From.ycord
## 1:              1 -2386282  8710048            12052   12213718   -4185952
## 2:              1 -2386282  8710048            11877   10613718   -3813952
## 3:              1 -2386282  8710048            11878   10713718   -3813952
## 4:              1 -2386282  8710048            11998   10213718   -4061952
## 5:              1 -2386282  8710048            12045   10213718   -4185952
## 6:              1 -2386282  8710048              727    9913718    7594048
##    BetaSim Sorenson  MNTD
## 1:  0.4561        1 4.556
## 2:  0.6484        1 4.055
## 3:  0.6484        1 4.055
## 4:  0.6484        1 4.055
## 5:  0.6484        1 4.055
## 6:  0.2571        1 4.037
```


Check out how fast that was remember dat has `r{dim(dat)}` rows, and i'm just on my old desktop.

Columns as functions
-----
To get both columns was a little tricky, i had to go to stack overflow: http://stackoverflow.com/questions/23521323/r-data-table-for-computing-summary-stats-across-multiple-columns


```r
dat.B <- dat[, list(c(To.OriginalRow, From.OriginalRow), BetaSim, MNTD, Sorenson)]

# it made a column that combined both To and From into a new column V1, not
# sure why it names it V1
setkey(dat.B, V1)

# make a function to compute tests
stat_test <- function(x) {
    c(mean(x[is.finite(x)]), var(x[is.finite(x)]))
}

dat.stat <- dat.B[, c(list(Stat = c("mean", "var")), lapply(.SD, stat_test)), 
    by = V1]

m.stat <- melt(dat.stat, id.var = c("V1", "Stat"))
head(cdat <- cast(m.stat, V1 ~ Stat + variable))
```

```
##      V1 mean_BetaSim mean_MNTD mean_Sorenson var_BetaSim var_MNTD
## 1     1       0.4928     3.858        0.9889     0.03545   0.3302
## 2    10       0.4928     3.858        0.9889     0.03545   0.3302
## 3   100       0.4928     3.858        0.9889     0.03545   0.3302
## 4 10055       0.4871     4.358        0.9990     0.01966   0.3255
## 5 10062       0.5259     3.768        0.9812     0.02866   0.4216
## 6   101       0.4928     3.858        0.9889     0.03545   0.3302
##   var_Sorenson
## 1    0.0075479
## 2    0.0075479
## 3    0.0075479
## 4    0.0002062
## 5    0.0074714
## 6    0.0075479
```

```r

# reset to data.table, the above might be slow because we went back to
# data.frame for casting, maybe this could be made better.

cdat <- data.table(cdat)
setkey(cdat, V1)
```


Okay, you may have noticed we lost the x y coordinates, i've tried this in a couple ways and i find it easier just to remerge the data.table. it will also help to show merging in the new data.table syntax

If you need to index that is NOT a logical statement, but gives the name add a ,with=FALSE to get data.frame like indexing.


```r
# get spatial info from the beginning table
Todat <- dat[, c("To.OriginalRow", "To.xcord", "To.ycord"), with = F]
Fromdat <- dat[, c("From.OriginalRow", "From.xcord", "From.ycord"), with = F]

# name the same colames, data.table style, call the cell V1 to equal to the
# cdat column name
setnames(Todat, colnames(Todat), c("V1", "X", "Y"))
setnames(Fromdat, colnames(Fromdat), c("V1", "X", "Y"))

# bind data.tables together
spdat <- rbind(Todat, Fromdat)

# remove duplicates, nice data.table function
spd <- unique(spdat)

setkey(spd, V1)
head(spd)
```

```
##       V1        X        Y
## 1:     1 -2386282  8710048
## 2:    10 -1486282  8710048
## 3:   100 -4186282  8338048
## 4: 10055 14513718 -1085952
## 5: 10062 15713718 -1085952
## 6:   101 -4086282  8338048
```


Merge
-----


```r

# merge is just tomerge[rows]
mergeD <- spd[cdat]
```


Visualize spatially
-----------------


```r
# i think we need a data.frame here, just give it the x y coords and colum n
# you want to be the values

df <- data.frame(mergeD[, c("X", "Y", "mean_BetaSim"), with = F])

library(raster)
```

```
## Loading required package: sp
```

```r
# create spatial points data frame
spg <- df
coordinates(spg) <- ~X + Y

# coerce to SpatialPixelsDataFrame

gridded(spg) <- TRUE
```

```
## Warning: grid has empty column/rows in dimension 1
## Warning: grid topology may be corrupt in dimension 1
## Warning: grid has empty column/rows in dimension 2
```

```r
# coerce to raster
rasterDF <- raster(spg)
rasterDF
```

```
## class       : RasterLayer 
## dimensions  : 109, 291, 31719  (nrow, ncol, ncell)
## resolution  : 1e+05, 124000  (x, y)
## extent      : -11536282, 17563718, -4743952, 8772048  (xmin, xmax, ymin, ymax)
## coord. ref. : NA 
## data source : in memory
## names       : mean_BetaSim 
## values      : 0.2419, 0.6044  (min, max)
```

```r
plot(rasterDF)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 


Instead of directly plotting values, another idea would be to take the quantiles of these means and then project those into space simulataneously, using the code similiar to ana's risk/climate velocity map so show where taxonomic/phylogenetic/trait are all high/low in a mapped product. Gabriel has that code.

