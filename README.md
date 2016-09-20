# RUBIC v1.0.2#
## What is RUBIC? ##
RUBIC is an R package.

RUBIC detects recurrent copy number aberrations using copy number breaks,
rather than recurrently amplified or deleted regions. This allows for a vastly simplified
approach as recursive peak splitting procedures and repeated re-estimation of the background model
are avoided. Furthermore, the false discovery rate is controlled on the level of called regions,
rather than at the probe level.

## <a name="summary"></a>Summary ##
- [How do I get set up?](#setup)
- [Case studies](#case-studies)
    - [Simplest execution](#simplest-execution)
    - [Execution break-down with partial saves](#execution-breakdown)
    - [Choosing gene locations](#gene-locations)
    - [Plot only selected genes](#selected-genes)
    - [Plot line segments instead of steps](#line-segments)
- [Example](#example)
- [API reference](#api)
    - [The rubic function](#rubic-function)
    - [The RUBIC class](#rubic-class)
        - [The estimate.parameters method](#estimate-parameters)
        - [The segment method](#segment)
        - [The call.events method](#call-events)
        - [The call.focal.events method](#call-focal-events)
        - [The save.focal.gains method](#save-focal-gains)
        - [The save.focal.losses method](#save-focal-losses)
        - [The save.plots method](#save-plots)
        - [The save method](#save)
- [File formats](#file-formats)
    - [The seg.cna file](#seg-file)
    - [The markers file](#markers-file)
    - [The sample file](#sample-file)
    - [The gene file](#gene-file)
- [FAQ](#faq)

## <a name="setup"></a>How do I get set up? ##
RUBIC depends on the following libraries:

  - `data.table (= 1.9.4 or >= 1.9.8)` read the [FAQ](#faq) for more details
  - `pracma`
  - `digest`
  - `ggplot2`
  - `grid`
  - `gtable`
  - `grDevices`

If you don't have them, you need to install them with (for example) the following command:
```r
install.packages(c('data.table', 'pracma', 'digest', 'ggplot2', 'gtable'))
```
*NOTE*: `grid` and `grDevices` are already included in modern R distributions.

Once you have installed all the required libraries, the following code should get you up and running:
```r
install.packages('/path/to/RUBIC_1.0.2.tgz', repos=NULL)
library(RUBIC)
```

If you are using a windows machine, you should use the package `RUBIC_1.0.2.zip`.
Alternatively, you can install the source package as follows:
```r
install.packages('/path/to/RUBIC_1.0.2.tar.gz', repos=NULL, source=TRUE)
library(RUBIC)
```

## <a name="case-studies"></a>Case studies ##

### <a name="simplest-execution"></a>Simplest execution ###
The simplest way to use the package is to initilize RUBIC with a minimal set
of parameters and ask to save the results as follows:
```r
rbc <- rubic(0.25, '/path/to/segcna.tsv', '/path/to/markers.tsv')
rbc$save.focal.gains('/here/I/want/the/focal_gains.tsv')
```
*IMPORTANT*: By default, RUBIC expects that the segments and the marker file have
headers while the samples file do not have any. Please adjust these defaults using
the options `seg.cna.header`, `markers.header`, and `samples.header` respectively.

*IMPORTANT*: Be aware that since November 9, 2015, Biomart suspended its web services,
this means that it is necessary to specify the `gene` option as follow when creating
the RUBIC object:
```r
rbc <- rubic(0.25, '/path/to/segcna.tsv', '/path/to/markers.tsv', genes='/path/to/genes.tsv')
```

*NOTE*: The first call to `save.{focal.gains,focal.losses,plots}` need to compute the results first.
Therefore it might need long time to complete.
```r
rbc$save.focal.losses('/here/I/want/the/focal_losses.tsv')
rbc$save.plots('/here/I/want/the/plots')
```
*NOTE*: `/here/I/want/the/plots` should point to a directory (it will be created if doens't exists).

### <a name="execution-breakdown"></a>Execution break-down with partial saves ###
Some times a full RUBIC run might be a very long and memory intensive job.
Therefore, it might be wise to run RUBIC in steps and save each partial RUBIC result.
In this way, if the computation environment become unreliable,
it will be possible to restore the partial results already computed.
Here is an example:
```r
rbc <- rubic(0.25, '/path/to/segcna.tsv', '/path/to/markers.tsv')
rbc$save('/here/I/want/the/post_processed_rubic.rds')
```
*IMPORTANT*: Be aware that since November 9 2015 Biomart suspended its web services,
this means that it is necessary to specify the `gene` option as follow when creating
the RUBIC object:
```r
rbc <- rubic(0.25, '/path/to/segcna.tsv', '/path/to/markers.tsv', genes='/path/to/genes.tsv')
```

At this point RUBIC preproccesed all the input. The `rbc` object has been saved in `post_processed_rubic.rds` and can
be restored using the standard `rbc <- readRDS('/here/I/want/the/post_processed_rubic.rds')` command.
```r
rbc$estimate.parameters()
rbc$save('/here/I/want/the/post_estimate_rubic.rds')

rbc$call.events()
rbc$save('/here/I/want/the/post_events_rubic.rds')

rbc$call.focal.events()
rbc$save('/here/I/want/the/final_rubic.rds')

rbc$save.focal.gains('/here/I/want/the/focal.gains.tsv')
rbc$save.focal.losses('/here/I/want/the/focal_losses.tsv')
rbc$save.plots('/here/I/want/the/plots')
```

### <a name="gene-locations"></a>Choosing gene locations ###
If the `gene` argument is not specified when calling `rubic`,
RUBIC will automagically download and use the human reference from Biomart.
Otherwise it is possible to specify a custom built in Biomart format as follow:
```r
rbc <- rubic(0.25, '/path/to/segcna.tsv', '/path/to/markers.tsv',
             genes='/path/to/custom_build.tsv')
```

### <a name="selected-genes"></a>Plot only selected genes ###
If `save.plots` is called without specifying the optional argument `genes`,
the function will plot the position of the same genes used for calling the
focal regions. It is possible to specify a different (usually more restrictive)
list of genes as follow:
```r
rbc$save.plots('/here/I/want/the/plots', genes='/path/to/custom_build.tsv')
```
In this case the argument `genes` is in the same format as in the `rubic` function.

### <a name="line-segments"></a>Plot line segments instead of steps ###
It is possible to plot only line segments corresponding to the segments instead
of steps as follow:
```r
rbc$save.plots('/here/I/want/the/plots', steps=F)
```

## <a name="case-studies"></a>Example ##
An example dataset has been included in the RUBIC package.To run the example
is possible to do the following.

1. Load the RUBIC package.
```r
library(RUBIC)
```
2. Retrieve the path of the example files included in the package.
```r
markers <- system.file("extdata", "markers.tsv", package="RUBIC")
samples <- system.file("extdata", "samples.tsv", package="RUBIC")
seg.cna <- system.file("extdata", "segcna.tsv", package="RUBIC")
genes <- system.file("extdata", "genes.tsv", package="RUBIC")
```
Alternatively these same data have been preloaded by the RUBIC package
and can be loaded in your R session running the following command.
```r
data(seg.cna, markers, samples, genes)
```

3. Crate and initialise a RUBIC object.
```r
rbc <- rubic(0.25, seg.cna, markers, samples, genes)
```
4. Save focal gains and losses to files.
```r
rbc$save.focal.gains('focal_gains.tsv')
rbc$save.focal.losses('focal_losses.tsv')
```
5. Create and save gain and losses plots for each chromosome.
```r
rbc$save.plots('plots')
```

## <a name="api"></a>API Reference ##
In order to use the RUBIC method it is necessary to create e initialise a new RUBIC object
using data and parameters specific to your analysis. Your can do that using the `rubic` function.

### <a name="rubic-function"></a>The rubic function ###

#### Description ####
With this function, it is possible to create and initialise a new RUBIC object.

#### Usage ####
```r
rubic(fdr, seg.cna, markers, samples = NULL, genes = NULL,
  amp.level = 0.1, del.level = -0.1, min.seg.markers = 1,
  min.mean = NA_real_, max.mean = NA_real_, min.probes = 260000,
  focal.threshold = 1e+07, seg.cna.header = T, markers.header = T,
  samples.header = F, col.sample = 1, col.chromosome = 2, col.start = 3,
  col.end = 4, col.log.ratio = 6, ...)
```
#### Arguments ####

fdr

  : a number indicating the wanted event based false discovery rate. For example 0.25.

seg.cna

  : Either a character string naming a file or a data.table containing the segmented CNA data in appropriate format. In case a file name is provided the file will be open and read as a data.table.

markers

  : Either a character string naming a file or a data.table containing the markers information in appropriate format. In case a file name is provided the file will be opened and read as a data.table.
  
samples

  : Either a character string naming a file or a character vector containing the sample IDs used for the analysis. In case a file name is provided the file will be open and read as a data.table. In case this parameter is not specified, all the samples presents in the seg.cna file will be used.
  
genes

  : Either a character string naming a file or a data.table containing the gene locations in appropriate format. In case a file name is provided the file will be open and read as a data.table. In case this parameter is not specified, RUBIC will attempt to download the most recent annotations using Biomart.

amp.level

  : A positive number specifying the threshold used for calling amplifications. The default value is set to 0.1.
  
del.level

  : A negative number specifying the threshold used for calling deletions. The default value is set to -0.1.
  
min.seg.markers

  : A positive integer specifying the number of probes allowed in each segment. The default value is set to 1, which means that no segments will be merged. If this parameter is set to a number larger than 1, all segments with less than min.seg.markers will be joined with adjacent segments until a segment with at least min.seg.markers will be formed.
  
min.mean

  : A number specifying the minimum mean copy number allowed. By default segments will not be filtered based on their minimum mean copy number.
  
max.mean

  : A number specifying the maximum mean copy number allowed. By default segments will not be filtered based on their maximum mean copy number.
  
min.probes

  : The minimum number of probes to be considered in the analysis.
  
focal.threshold

  : The maximum length of a recurrent region to be called focal. Only regions smaller than focal.threshold bases will be called focal. By default 10000000 bases.
  
seg.cna.header

  : A logical value indicating whether the seg.cna file contains the names of the variables as its first line. The default is set to TRUE.
  
markers.header

  : A logical value indicating whether the markers file contains the names of the variables as its first line. The default is set to TRUE.
  
samples.header

  : A logical value indicating whether the samples file contains the names of the variables as its first line. Please notice that while the default value of seg.cna.header and markers.header is set to TRUE, RUBIC does not expect a header by default for the samples file (and, therefore, this value is set by default to FALSE).
  
col.sample

  : The number of the column containing the sample name in the seg.cna file. By default, RUBIC expect the sample name to be in the 1st column.
  
col.chromosome

  : The number of the column containing the chromosome name in the seg.cna file. By default, RUBIC expects the chromosome name to be in the 2nd column.
  
col.start

  : The number of the column containing the start position of each segment in the seg.cna file. By default, RUBIC expects the start position to be in the 3rd column.
  
col.end

  : The number of the column containing the end position of each segment in the seg.cna file. By default, RUBIC expects the end position to be in the 4th column.
  
col.log.ratio

  : The number of the column containing the log.ratio name in the seg.cna file. By default, RUBIC expects the log.ratio to be in the 6th column.

#### Details ####
In order to use the RUBIC method it is necessary to create e initialise a new RUBIC object using data and parameters specific to your analysis. This function reads and preprocesses all the data and returns an object of the RUBIC class. It is possible to start the analysis using the method estimate.parameters or any other method of the RUBIC class.

#### Value ####
A new RUBIC object.

#### See Also ####
[RUBIC-class](#rubic-class)

### <a name="rubic-class"></a>The RUBIC class ###
The RUBIC class is a ReferenceClass that contains all the data needed and generated by the RUBIC algorithm.

#### Methods ####

<a name="estimate-parameters"></a>
```r
estimate.parameters(quiet = T, test.env = NULL)
```
Estimates the parametersr necessary for segmentation and event calling.

_Arguments:_

quiet

  : A logical value that controls whether the parameters estimation should print on the console its progress. By default (`TRUE`), it runs in quiet mode suppressing all messages.

test.env

  : An environment which must contain a matrix named `random.matrix`. The number of columns in `random.matrix` should be equal to the number of samples to be analyzed. If this parameter is not NULL,
  each time RUBIC will need to perform a random permutation it will draw the random numbers from the matrix contained in this environment instead of generating new ones. This parameter is used internally to
  perform code testing.

_Details:_

This is one of the most computationally intensive tasks. For particularly large datasets, is good practice to [`save`](#save) the RUBIC object after completing this task. See the [execution break-down use case](#execution-breakdown) for more information.
<br><br>

<a name="segment"></a>
```r
segment()
```
Generate positive and negative segments.

_Details:_

This is a step necessary before calling [`call.events`](#call-events). If this method has not been called already, it will be called automatically before [`call.events`](#call-events). This method call automatically [`estimate.parameters`](#estimate.parameters) if it has not been previously called.
<br><br>

<a name="call-events"></a>
```r
call.events()
```
Call recurrent events.

_Details:_

This is a step necessary before calling [`call.focal.events`](#call-focal-events). If this method has not been called already, it will be called automatically before [`call.focal.events`](#call-events). This method call automatically [`segment`](#segment) if it has not been previously called.
<br><br>

<a name="call-focal-events"></a>
```r
call.focal.events(genes = NULL)
```
Call recurrent events.

_Arguments:_

genes

  : Either a character string naming a file or a data.table containing the gene locations in the [appropriate format](#gene-file). In case a file name is provided the file will be open and read as a data.table. In case this parameter is specified, it will override the `gene` argument passed to the [`rubic`](#rubic-function) function. In case neither this parameter nor the `gene` argument passed to [`rubic`](#rubic-function) was specified, RUBIC will attempt to download the most recent annotations using Biomart.

_Details:_

This is the last step in the method which computes the actual significant focal events. If this method has not been called already, it will be called automatically before calling any output function, such as [`save.focal.gains`](#save-focal-gains), [`save.focal.losses`](#save-focal-losses), and [`save.plots`](#save-plots). This method call automatically [`call.events`](#call-events) if it has not been previously called.
<br><br>

<a name="save-focal-gains"></a>
```r
save.focal.gains(file)
```
Save focal gains to a file in TSV format.

_Arguments:_

file
  : A [connection](https://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html) or the name of the file where RUBIC will save the focal gains in TSV format.

_Details:_ This method call automatically [`call.focal.events`](#call-focal-events) if it has not been previously called.
<br><br>

<a name="save-focal-losses"></a>
```r
save.focal.losses(file)
```
Save focal losses to a file in TSV format.

_Arguments:_

file
  : A [connection](https://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html) or the name of the file where RUBIC will save the focal losses in TSV format.

_Details:_ This method call automatically [`call.focal.events`](#call-focal-events) if it has not been previously called.
<br><br>

<a name="save-plots"></a>
```r
save.plots(dir, genes = NULL, steps = T, width = 11, height = 5)
```
Save gains and losses plots for each chromosome.

_Arguments:_

dir

  : The directory in which all plots will be created. If the directory doesn't
  exist, it will be created.

genes

  : Either a character string naming a file or a data.table containing the locations of the genes of interest in the [appropriate format](#gene-file). This parameter can be used to plot only a subset of genes
  in which the user is particularly interested, without changing the results of [`call.focal.events`](#call-focal-events).
  
steps

  : A logical value specifying whether to plot only the lines corresponding to the segments instead of the full steps (i.e. including the vertical lines).

width

  : The width of each plot in inches.
  
height

  : The height of each plot in inches.

_Details:_

This method creates and saves two plots for each chromosome; one plot showing the gains and one plot
showing the losses. In each plot is shown the location of the genes used to compute the focal events. However, it is possible to plot a different set of genes using the `genes` parameter. See the relative [use case](#selected-genes) for more information. This method call automatically [`call.focal.events`](#call-focal-events) if it has not been previously called.
<br><br>

<a name="save"></a>
```r
save(file)
```
Save the current state of the RUBIC object to file.

_Arguments:_

file
  : A [connection](https://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html) or the name of the file where the RUBIC object is going to be saved.

_Details:_

It is possible to restore a previously saved RUBIC object using the [readRDS](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) function. See the [execution break-down use case](#execution-breakdown) for more information.

## <a name="file-formats"></a>File formats ##
RUBIC expects either files in a specific format or data.tables (or object that can be coerced in data.tables) with specific columns. It is possible and recommended to instantiate a new RUBIC object by using the convenience function [`rubic`](#rubic-function).

### <a name="seg-file"></a>The seg.cna file ###
The seg.cna is in TSV format and must contain at least a column for sample name, chromosome, segment start, segment end, log ratio. By default,
it expects sample names to be the 1st column, chromosome the 2nd, segment start the 3rd, segment end  the 4th, and log ratio the *6th*. This is because
in the 5th column there is often the count of probes per segment.
If that is not the case it is possible to specify the positions using the optional arguments
`col.sample, col.chromosome, col.start, col.end, col.log.ratio` to the function `rubic()` as in this example:
```r
rbc <- rubic(0.25, '/path/to/segcna.tsv', '/path/to/markers.tsv', col.log.ratio=5)
```
In the example above the 6th column will be discarded (if present) and the log.ratio will be read from column number 5.
If the user chooses to provide a data.table instead of a file name, such data.table must contain at least five columns with the following names:
'Sample', 'Chromosome', 'Start', 'End', 'LogRatio'. Like in the following example:
```
       Sample Chromosome    Start       End  LogRatio
    1:  S0001          1  3218610   5838773  0.408563
    2:  S0001          1  5844802   9012691  0.845405
    3:  S0001          1  9012737  28712545  0.838700
    4:  S0001          1 28716679  33008483  2.682358
    5:  S0001          1 33019130  33032412  1.244822
   ---                                               
15728:  S0040         22 44534620  44538370 -0.364691
15729:  S0040         22 44545522  45407499  0.000000
15730:  S0040         22 45424360  45425535 -0.659315
15731:  S0040         22 45425671  49331012  0.000000
15732:  S0040          X  3157107 154905589  0.000000
```

### <a name="markers-file"></a>The markers file ###
The markers file indicates the exact locations of measurement probes (markers) for the given platform.  For sequencing data, copy number values are often estimated with fixed bin sizes (prior to segmentations). In this case each marker should be associated with a bin and the center genomic position of the bin.
The file must be in TSV format and must contain at least 3 columns. The 1st must contain the probe name, the 2nd the chromosome, the 3rd must contain the
location on the chromosome. Otherwise, if the user choose to provide a data.table instead of a file name, such data.table must contain at least
three columns with the following names: 'Name', 'Chromosome', 'Position'. Like in the following example:
```
          Name Chromosome  Position
    1: P000001          1   3252007
    2: P000002          1   3714033
    3: P000003          1   4083031
    4: P000004          1   4214828
    5: P000005          1   4309034
   ---                             
 9996: P009996          X 153269127
 9997: P009997          X 153376055
 9998: P009998          X 153392576
 9999: P009999          X 153641851
10000: P010000          X 153689045
```

### <a name="sample-file"></a>The sample file ###
The sample file must contain one column with sample names. If the file contains more column only the first will be used.

### <a name="gene-file"></a>The gene file ###
The gene file format is the same provided by Biomart (TSV). Therefore it must contains the following column headers: 'Ensembl Gene ID', 'Associated Gene Name',
'Chromosome Name', 'Gene Start (bp)', 'Gene End (bp)'. Column order doesn't matter as long as the column names are identical.
If the user chooses to provide a data.table instead of a file name, such data.table must contain at least columns with the following names,
'ID', 'Name', 'Chromosome', 'Start', 'End', like in the following example:
```
                    ID       Name               Chromosome     Start       End
    1: ENSG00000252830    5S_rRNA                        1 143439605 143439714
    2: ENSG00000263418    5S_rRNA CHR_HSCHR6_MHC_MANN_CTG1  32117921  32118041
    3: ENSG00000274097    5S_rRNA                       11 102057854 102057960
    4: ENSG00000201285    5S_rRNA                        X 148008100 148008215
    5: ENSG00000278457    5S_rRNA               KI270442.1    380608    380726
   ---                                                                        
66588: ENSG00000274075     uc_338                       12  27824207  27824382
66589: ENSG00000277035     uc_338                        2  28134408  28134585
66590: ENSG00000276673     uc_338                        8  76733302  76733476
66591: ENSG00000275939     uc_338                        5 131440031 131440196
66592: ENSG00000213076 yR211F11.2                        6 158921271 158922150
```
However, the 'ID' (or 'Ensembl Gene ID') column can be empty.


## <a name="faq"></a>FAQ ##

### I'm using Mac OSX and `save.plots()` doesn't seem to work. What's wrong? ###
RUBIC uses cairo_ps to produce high-quality EPS plots with transparencies. Most likely you have to (re-)install XQuartz (http://xquartz.macosforge.org).

### While I run RUBIC, I get the following error from Biomart: ###
```
Error in download.file(input, tt, mode = "wb") : 
  cannot open URL 'http://www.biomart.org/biomart/martservice?
```
Unfortunately, as of November 9 2015 [Biomart](http://www.biomart.org) suspended
its web services (click [here](http://www.biomart.org/notice.html) to know more).
This means that RUBIC cannot download automatically gene locations
from Biomart until they will restore their services. Meanwhile, we suggest to
manually specify the parameter `gene` in the [`rubic`](#rubic-function) function.

### While I run RUBIC, I get the following error: ###
```
Error in set(i, j = lc, value = newfactor) : 
  .SD is locked. Updating .SD by reference using := or set are reserved for future use. Use := in j directly. Or use copy(.SD) as a (slow) last resort, until shallow() is exported.
```
This problem is due to a bug in the package `data.table` described [here](https://github.com/Rdatatable/data.table/issues/1341), which affects at least `data.table` version 1.9.6. If you are using `data.table` version 1.9.6, we stuggest to downgrade your `data.table` to 1.9.4 or upgrade to version 1.9.8 (if possible), in which this bug should be corrected.

### I'm running RUBIC using `Rscript` and I get the following error: ###
```
Error in initFields(scales = scales) :
  could not find function "initRefFields"
Calls: <Anonymous> ... initialise -> initialise -> <Anonymous> -> initFields
Execution halted
```
This happens because `Rscript` doesn't attach the `methods` package (while the normal `R` environment does). Therefore, you need to include the following line at the beginning of your script.
```r
library('methods')
```