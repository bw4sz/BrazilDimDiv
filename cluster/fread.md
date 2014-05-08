Data table crashes
========================================================

I am having a similar issue to this post: http://r.789695.n4.nabble.com/fread-crash-td4683394.html

The file is about 550MB, i'm unsure how many rows it actually is (several million).

When i try to run fread, Rstudio just crashes with no error. 
I can read in up to about 15 rows


```r

require(data.table)
```

```
## Loading required package: data.table
```

```r

# env dist table

env <- fread("EnvData.csv", nrows = 15, verbose = TRUE)
```

```
## Input contains no \n. Taking this to be a filename to open
## File opened, filesize is  0.543B
## File is opened and mapped ok
## Detected eol as \n only (no \r afterwards), the UNIX and Mac standard.
## Using line 30 to detect sep (the last non blank line in the first 'autostart') ... sep=','
## Found 4 columns
## First row with 4 fields occurs on line 2 (either column names or first row of data)
## Some fields on line 2 are not type character (or are empty). Treating as a data row and using default column names.
## Count of eol after first data row: 15989212
## Subtracted 0 for last eol and any trailing empty lines, leaving 15989212 data rows
## nrow limited to nrows passed in (15)
## Type codes: 4113 (first 5 rows)
## Type codes: 4113 (after applying colClasses and integer64)
## Type codes: 4113 (after applying drop or select (if supplied)
## Allocating 4 column slots (4 - 0 NULL)
##    0.000s (  0%) Memory map (rerun may be quicker)
##    0.000s (  0%) sep and header detection
##    0.721s (100%) Count rows (wc -l)
##    0.000s (  0%) Column type detection (first, middle and last 5 rows)
##    0.000s (  0%) Allocation of 15x4 result (xMB) in RAM
##    0.000s (  0%) Reading data
##    0.000s (  0%) Allocation for type bumps (if any), including gc time if triggered
##    0.000s (  0%) Coercing data already read in type bumps (if any)
##    0.000s (  0%) Changing na.strings to NA
##    0.721s        Total
```

```r
head(env)
```

```
##    V1 V2 V3     V4
## 1:  1  2  1  249.3
## 2:  2  3  1  536.9
## 3:  3  4  1 1161.8
## 4:  4  5  1 1234.0
## 5:  5  6  1 1513.4
## 6:  6  7  1 1757.1
```


However when i run fread with more than 20 rows, it crashes Rstudio. 


```r
# not run
env <- fread("EnvData.csv", nrows = 25, verbose = TRUE)
```


*verbose on the error output reads:*

Input contains no \n. Taking this to be a filename to open

File opened, filesize is  0.543B

File is opened and mapped ok

Detected eol as \n only (no \r afterwards), the UNIX and Mac standard.

Using line 30 to detect sep (the last non blank line in the first 'autostart') ... sep=','

Found 4 columns

First row with 4 fields occurs on line 2 (either column names or first row of data)

Some fields on line 2 are not type character (or are empty). Treating as a data row and using default column names.

Count of eol after first data row: 15989212

Subtracted 0 for last eol and any trailing empty lines, leaving 15989212 data rows

nrow limited to nrows passed in (25)

Type codes: 4113 (first 5 rows)

Type codes: 4113 (+middle 5 rows)

*Look at the file, nothing seems wrong*


```r

env <- read.csv("EnvData.csv", nrows = 25)

env
```

```
##    V1 V2     V3
## 1   2  1  249.3
## 2   3  1  536.9
## 3   4  1 1161.8
## 4   5  1 1234.0
## 5   6  1 1513.4
## 6   7  1 1757.1
## 7   8  1 2176.7
## 8   9  1 2644.0
## 9  10  1 3033.3
## 10 11  1 3721.2
## 11 12  1 4432.8
## 12 13  1 4609.6
## 13 14  1 5378.8
## 14 15  1 5953.6
## 15 16  1 5913.9
## 16 17  1 6281.3
## 17 18  1 6669.7
## 18 19  1 6449.7
## 19 20  1 6218.4
## 20 21  1 6493.4
## 21 22  1 6056.6
## 22 23  1 5275.8
## 23 24  1 4605.2
## 24 25  1 3153.9
## 25 26  1 2532.1
```


