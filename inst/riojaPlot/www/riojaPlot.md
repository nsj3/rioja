---
title: "riojaPlot Version 1.0"
output:
  html_document: default
  pdf_document: default
---
# Steve Juggins Oct 2020

## Introduction

riojaPlot is a web app for plotting stratigraphic data found in micropalaeontological or palaeoecological studies. It can plot multiple biological and non-biological proxies against depth or age and allows trends to be visualized using a range of plot types.  

## Input data

Data should be in an Excel spreadsheet, with a single row of column names.  These names will be used to label the x-variables so choose them carefully.  Any variables sharting with the words "Depth" or "Age" will be used as y-variables.  For example, "Depth_in_m"" or "AgeBP" are fine, "Years BP" or "core depth" are not.  This allows one to plot the data using different age models by labelling the columns "Age01", "Age02" etc.  If you do not include a Depth or Age variable samples data be plotted against sample number.  A column starting with the word "Label" can contain text that will be used to label samples on the the y-axis. All other columns should be numeric and are used as x-variables.  Data are plotted 'as is' and in the order they appear in the file.  That is, if you want to plot a diagram of microfossil counts expressed as relative abundances you should supply percentage data with columns in correct the order.  An example spreadsheet containing pollen stratigraphic data from Birks & Mathews (1978) can be downloaded from the Help tab.  A screenshot of the data is listed at the end of this help file.

## Getting started

Upload a Excel file and select the worksheet to use, or check the box **Use example dataset**.  If your data look like biostratigraphic data with row sums in the range 50-150 riojaPlot will assume they are percentages and scale each curve appropriately, and de-select any variables with a maximum value less than 2%. If not, each curve will have equal width and be scaled from min-max data values (option **Scale for %** under the **Settings tab**).  If you have more than 50 variables on the first 50 will be selected.  You can change this selection on the **Variables** tab.  

## Options

### Variables

Select a variable for the **Y axis** (if you have included more than one *Depth*, *Age* or *Label* columns).  riojaPlot will try to guess the type of data you have.  If the data look like biological counts transformed to percentages then variables with a maximum value of less than 2.0 will be de-selected to avoid over-crowding.  Use **Select X variables** this option to remove or add variables to the diagram.  Use the **Scale for %** checkbox to force riojaPlot to plot the data as percentage (or not).  

**Auto Select X vars** can be used with % scaling to automatically select variables. riojaPlot calculates the maximum value of each variable and selects those variables whose maximum falls between the min / max cutoff.  This is useful for quickly de-selecting rare taxa in large datasets.

### Settings

- **Style**: Choose to plot each variable with lines, symbols or silhouettes (filled curves).  Some combinations do not make sense or look ugly so choose wisely.  Shilouettes may not plot properly if there are missing values in the data.

- **Show bar**: Show horizontal bars, either alone or superimposed on a line, symbol or silhouette plot.  **Curve** extends the bar to the data value, **Full** extends the bar the full width of the plotting area.  Try the latter with silhouettes and **Bars on top** unchecked.

- **Bars on top**: plot bars on top of silhouette or below. 

- **Settings**: **Reverse Y axis** does what it says.  Usually this option should be selected so depths / ages run from young at the top to old at the bottom.  De-select if your data are in years CE.  **Show min/max**: with **Scale for %** unchecked, shows either min / max or multiple values on x-axes (to prevent label crowding).  **Auto sort vars**:  sorts variables to highlight sequence from those with high values at base on left, to those with high values at top on right (can be useful to visualize trends in percentage biostratigraphic data).

- **Exaggeration**: Show an exaggerated shilouete, either in the **exaggeration colour** selected on the colour tab, or, if **Auto Col** is selected, a pale version of the shilouette colour.  The latter is useful if you groups of variables with different colour shilouettes. **Exag mult**. is the multiplier for the exaggerated curve.

### Colours

Select colour for lines, bars, silhouettes, symbols, zones and exaggeration.

### Sizes

Adjust axis, label font size and label rotation.

### Zonation

Add a zonation (constrained clustering) to the diagram using CONISS (Grimm 1987).  Optional show zones on the diagram with the number of zones determined automatically using a broken-stick model (Bennett 1996) or chosen manually. Data can be transformed prior to calculation of dissimilarities.  **Sqrt** is good for biological data expressed as percentage relative abundance (Legendre & Gallagher, 2001).  **Scale** scales each variable to zero mean and unit standard deviation and is good for non-biological data with different units of measurement.

### Groups

Assign each variable to a group and plot silhouettes in a different colour for each group.

## Save the plot

Save the plot as a pdf, png or svg file.  Scalable vector graphics (svg) format is good for importing into Powerpoint or Word.

## Details

**riojaPlot** is powered by the function `strat.plot` in the R package [rioja](https://cran.r-project.org/web/packages/rioja/index.html).  The web interface is built using [shiny](https://shiny.rstudio.com/) and [shinydashboard](https://rstudio.github.io/shinydashboard/index.html).  

## Example dataset

Pollen stratigraphic data from the Abernethy Forest, Scotland, spanning approximately 5500 - 12100 BP (from Birks & Mathews 1978).

![Abernethy Forest Pollen Stratigraphic data (from Birks & Matthews, 1978)](riojaPlot01.jpg)

## Contact

Bug reports and suggestions for improvement to Steve Juggins: Stephen.Juggins@ncl.ac.uk.

## References

Bennett, K (1996) Determination of the number of zones in a biostratigraphic sequence. *New Phytologist*, **132**, 155-170.

Birks, HH & Mathews, RW (1978) Studies in the vegetational history of Scotland V. Late Devensian and early Flandrian macrofossil stratigraphy at Abernethy Forest, Invernessshire. *New Phytologist*, **80**, 455-84.

Grimm, EC (1987) CONISS: A FORTRAN 77 program for stratigraphically constrained cluster analysis by the method of incremental sum of squares. *Computers & Geosciences*, **13**, 13-35.

Legendre, P & Gallagher E (2001) Ecologically meaningful transformations for ordination of species data. *Oecologia*, **129**,  271-280.
	
