suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(markdown))
suppressPackageStartupMessages(library(rhandsontable))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(rioja))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(graphics))

email <- tags$html( tags$body( a(href="mailto:Stephen.Juggins@ncl.ac.uk")))
riojaURL <- a("rioja", href="https://cran.r-project.org/web/packages/rioja/index.html")
SJURL <- a("Steve Juggins", href="mailto:Stephen.Juggins@ncl.ac.uk")
inputWidth <- "300px"
numericWidth <- "100px"
numericWidth2 <- "75px"
version <- "1.01"
colWidth <- "120px"
colList <- c(list('black' ),
    as.list(brewer.pal(9, "Blues")[-(1:3)]),
    as.list(brewer.pal(9, "Greens")[-(1:3)]),
    as.list(brewer.pal(9, "Reds")[-(1:3)]),
    as.list(brewer.pal(9, "Greys"))
  )

removeMsg <- function() { lapply(errorMsg, function(x) { removeNotification(x); errorMsg <- list(); }) }

addNotification <- function(msg, ...)
{
   id <- showNotification(msg, ...)
   nErr <- length(errorMsg)
   errorMsg[[nErr+1]] <<- id
}

N <- 0

errorMsg <- list()
currentSheet <- ""
currentFile <- ""
fn <- ""
fn2 <- ""
mydata <- NULL
maxGroups <- 6
messageDur <- 20

selectVars <- function(session, input)
{
   minC <- input$minCut
   maxC <- input$maxCut

   removeMsg()
   if (is.na(minC) | minC < 0 | minC > 20 | is.na(maxC) | maxC < 50) {
     addNotification("Please enter sensible values for Min cutoff and Max cutoff", 
                     duration=messageDur, type="message", closeButton=TRUE)
     return()
   }
   selTaxa <- colnames(mydata$spec)    
   mx <- apply(mydata$spec, 2, max, na.rm=TRUE)
   del <- mx > minC & mx < maxC
   
   if (sum(del) > 80) {
      removeMsg()
      msg <- paste0("This selection will result in ", sum(del), " variables. This is too many to plot.  
                    Change the cutoff values and try again.")
      addNotification(msg, duration=messageDur, type="message", closeButton=TRUE)
      return()
   }   
   
   selTaxa <- selTaxa[del]
   updatePickerInput(session, 'selTaxa', choices=colnames(mydata$spec), selected=selTaxa)
   if (sum(del) < length(del)) {
      msg <- paste0("Your data has ", length(mx), " variables. Only ", length(selTaxa), " with maximum abundance 
                    between ", minC, "% and ", maxC , "% have been plotted.  Use the selector under the 'Variables' tab to 
                    add or remove additional variables.")
      addNotification(msg, duration=messageDur, type="message", closeButton=TRUE)
   }
}

readIt <- function(fName, sheet, input, session, forceRead=TRUE) #, forceByCol=FALSE, forceAutoSelect=FALSE) 
{
#   print(paste0("File: ", fName))
#   print(paste0("Sheet: ", sheet))
#   print(paste0("CurrentFile: ", currentFile))
#   print(paste0("Currentsheet: ", currentSheet))
   if (is.character(sheet) & nchar(sheet)==0)
      return()
   if ((currentSheet == sheet & currentFile == fName) & !forceRead)
      return()
# reset UI elements    
   removeMsg()
   isolate(updateNumericInput(session, 'yMin', value=NA))
   isolate(updateNumericInput(session, 'yMax', value=NA))
   isolate(updateNumericInput(session, 'yInterval', value=NA))
   isolate(updateTextInput(session, 'yLabel', value=""))
   isolate(updateCheckboxGroupInput(session, 'ySettings', selected=1))
   isolate(updateCheckboxInput(session, 'doClust', value=FALSE))
    
   currentSheet <<- sheet
   currentFile <<- fName
   d <- read_excel(fName, sheet=sheet)
    
   sel1 <- toupper(substring(colnames(d), 1, 5)) %in% "DEPTH"
   sel2 <- toupper(substring(colnames(d), 1, 3)) %in% "AGE"
   sel <- sel1 | sel2
# check column types:
   cols <- sapply(d, mode)
   nonChar = which(cols != "numeric")
   if (length(nonChar) > 0) {
      ccols <- colnames(d)[nonChar]
      dropCol <- which (toupper(substring(ccols, 1, 5)) != "LABEL")
      if (length(dropCol)>0)
        nonChar <- nonChar[dropCol]
      msg <-  paste0("Spreadsheet has ", length(nonChar), " non-numeric column(s). They have been removed.")
      addNotification(msg, duration=messageDur, type="message", closeButton=TRUE)
      d <- d[, -nonChar, drop=FALSE]
   }
   
# delete cols with all missing values   
   del <- apply(is.na(d), 2, sum)
   
   del <- which(del == nrow(d))
   if (any(del)) {
      msg <-  paste0("Spreadsheet has ", sum(!del), " empty columns. They have been removed.")
      addNotification(msg, duration=messageDur, type="message", closeButton=TRUE)
      d <- d[, !del, drop=FALSE]
   }     

# check column headings
   colNames <- suppressWarnings(sum(!is.na(as.numeric(colnames(d)))))
   if (colNames > 0) {
      msg <-  paste0("Spreadsheet has ", colNames, " column names that are numbers. ", 
                     "Are you sure your data is formatted with variables in columns and observations in rows?")
      addNotification(msg, duration=30, type="error", closeButton=TRUE)
   }
   sel1 <- toupper(substring(colnames(d), 1, 5)) %in% "DEPTH"
   sel2 <- toupper(substring(colnames(d), 1, 3)) %in% "AGE"
   sel3 <- toupper(substring(colnames(d), 1, 5)) %in% "LABEL"
   sel <- sel1 | sel2 | sel3
   if (any(sel)) {
      dd <- NULL
      if (any(sel1 | sel2)) {
         dd <- d[, sel1 | sel2, drop=FALSE]
      }
      if (any(sel3)) {
        d3 <- sapply(d[, sel3, drop=FALSE], as.character)
        if (!is.null(dd))
           dd <- cbind(dd, d3)
        else 
           dd <- d3
      }
      mydata$chron <<- dd
   } else {
      mydata$chron <<- data.frame(SampleNo=1:nrow(d))
   }
   mydata$spec <<- d[, !sel]
   nMiss <- sum(is.na(d))
   if (nMiss > 0) {
      msg <- paste0("Spreadsheet has ", nMiss, " missing value(s). ",
                    "Zonation is disabled.")
      addNotification(msg, duration=messageDur, type="message", closeButton=TRUE)
      shinyjs::disable(selector = "#showZones")
      shinyjs::disable(selector = "#nZones")
      shinyjs::disable(selector = "#doClust")
      shinyjs::disable(selector = "#dataTransformation")
   } else {
      shinyjs::enable(selector = "#showZones")
      shinyjs::enable(selector = "#nZones")
      shinyjs::enable(selector = "#doClust")
      shinyjs::enable(selector = "#dataTransformation")
   }
   yvars <- colnames(mydata$chron)
   yvar <- yvars[1]
   isolate(updateSelectInput(session, 'yvar', choices=yvars))

# delete columns with no data
   
   mx <- apply(abs(mydata$spec), 2, sum, na.rm=TRUE)
   if (any(mx < 10e-6)) {
      del <- which(mx < 10e-6)
      mydata$spec <- mydata$spec[, -del]
      msg <- paste0("Your data has ", length(del), " variables with no data. They have been removed.")
      addNotification(msg, duration=messageDur, type="message", closeButton=TRUE)
   }

# determine if % data   

   mx <- apply(mydata$spec, 1, sum, na.rm=TRUE)
   r <- range(mx, na.rm=TRUE)
   scalePC <- TRUE
   if (r[1] < 50 | r[2] > 250) {
      scalePC <- FALSE
   } 

   selTaxa <- colnames(mydata$spec)
   
   if (scalePC) {
     shinyjs::enable("autoSelect")
     shinyjs::enable("minCut")
     shinyjs::enable("maxCut")
     isolate(updateCheckboxGroupInput(session, "xSettings", selected=c(1, 2)))
     isolate(updateCheckboxGroupInput(session, "exag", selected=c(1)))
     isolate(updateCheckboxGroupInput(session, "style", selected=c(3)))
     mx <- apply(mydata$spec, 2, max, na.rm=TRUE)
     minCut <- 2
     maxCut <- 150
     del <- mx > minCut & mx < maxCut
     if (sum(del) < 50) {
       selTaxa <- selTaxa[del]
     } else {
       minCut <- 5
       del <- mx > minCut & mx < maxCut
       selTaxa <- selTaxa[del]
     }
     if (length(mx) > length(selTaxa)) {
        msg <- paste0("Your data has ", length(mx), " variables. Only ", length(selTaxa), " with maximum abundance 
                      between ", minCut, "% and ", maxCut , "% have been plotted.  Use the selector under the 'Variables' tab to 
                      add or remove additional variables.")
        addNotification(msg, duration=messageDur, type="message", closeButton=TRUE)
     }
   } else {
     isolate(updateCheckboxGroupInput(session, "xSettings", selected=character(2)))
     isolate(updateCheckboxGroupInput(session, "style", selected=c(1)))
     isolate(updateCheckboxGroupInput(session, "exag", selected=character(0)))
     shinyjs::disable("autoSelect")
     shinyjs::disable("minCut")
     shinyjs::disable("maxCut")
     nSel <- length(selTaxa)
     if (nSel > 50) {
       selTaxa <- selTaxa[1:50]
       msg <- paste0("Your data has ", nSel, " variables. Only the first 50 variables have been plotted. ",  
                    "Use the selector under the 'Variables' tab to add or remove additional variables. ",  
                    "If the data are percentages use the Auto select X vars button.")
       addNotification(msg, duration=messageDur, type="message", closeButton=TRUE)
     }
   }
   sel2 <- colnames(mydata$spec) %in% selTaxa
   updatePickerInput(session, 'selTaxa', choices=colnames(mydata$spec), selected=selTaxa)
   group <- rep("1", length(sel2))
   mydata$group <<- data.frame(Name=colnames(mydata$spec), Group=group, Selected=sel2, stringsAsFactors = FALSE)
   set_data("varList", row=1:length(sel2), col=1:3, mydata$group, session)
}

plotIt <- function(mydata, input, session, fileFlag=FALSE) 
{
  if (is.null(mydata$spec) | is.null(mydata$chron) | is.null(input$selTaxa))
     return();
  
#  print(paste("Plot It", N))
#  N <<- N + 1

  selTaxa <- input$selTaxa
  yvar <- input$yvar
  style <- list()

  style$scalePC <- 1 %in% input$xSettings
  style$scaleMinMax <- 2 %in% input$xSettings
  
  style$yrev <- 1 %in% input$ySettings

  style$autoOrder <- "none"
  if (3 %in% input$xSettings)
    style$autoOrder <- "bottomleft"

  style$bar <- FALSE
  style$line <- FALSE
  style$poly <- FALSE
  style$symbol <- FALSE
  doZone <- FALSE
  
  if (1 %in% input$style)
    style$line <- TRUE
  if (2 %in% input$style)
      style$symbol <- TRUE
  if (3 %in% input$style)
    style$poly <- TRUE

  if (input$hLine == "None")
    style$bar <- FALSE
  else if (input$hLine == "Curve")
    style$bar <- TRUE
  else  
    style$bar <- "Full"
  
  style$yMin <- input$yMin
  style$yMax <- input$yMax
  style$yInterval <- input$yInterval
  style$exag <- 1 %in% input$exag
  
  doZone <- input$doClust
    
  style$lwdBar <- input$barSize;
  style$colBar <- input$hLineCol
  style$SymbSize <- input$symbSize

  clust <- NULL
  showZones <- 0
  
   d <- mydata$spec
   if (!is.null(selTaxa) & length(selTaxa) > 0) {
     if (any(!(selTaxa %in% colnames(d)))) {
        print("selTaxa names don't match")
        return()
     }
     d <- d[, selTaxa]
   }

  if (!is.null(yvar) & (nchar(yvar) > 1)) {
     yvar <- mydata$chron[, yvar, drop=TRUE]
  }
    else {
     yvar <- mydata$chron[, 1, drop=TRUE]
  }
  if (doZone) {
     d2 <- d
     if (input$dataTransformation == "Sqrt") {
        d2 <- sqrt(d)
     } 
     if (input$dataTransformation == "Scale") {
        d2 <- scale(d, TRUE, TRUE)
     }
     diss <- dist(d2)
     clust <- chclust(diss)
     if (input$showZones == "No")
        showZones <- 0
     else if (input$showZones == "Auto") {
        x <- bstick(clust)
        x2 <- x$dispersion <= x$bstick
        showZones <- which(x2)[1]
        if (showZones < 2) {
         addNotification("There are no significant zones in these data.", duration=30, type="message", closeButton=TRUE)
      }
    } else {
       showZones = input$nZones
    }
  } 

  groupColours <- rep(input$fillCol, nrow(d))
  if (!is.null(input$varList)) {
    newG <- hot_to_r(input$varList)
    newG <- newG[newG$Name %in% input$selTaxa, ]
    groups <- as.numeric(newG$Group)   
    cnames <- c(input$groupCol1, input$groupCol2, input$groupCol3, input$groupCol4, input$groupCol5, input$groupCol6)
    groupColours <- cnames[groups]
  }    

  nms <- colnames(d)
  yTop <- 0.8
  xRight = 1.0
  yMax <- max(sapply(nms, function(x) strwidth(x, units="figure", cex=input$nameSize))) 
  plotRatio<- session$clientData$output_myPlot_width / session$clientData$output_myPlot_height
  yTop <- 1.0 - (yMax * plotRatio * cos(pi/180 * (90-input$nameAngle))) - 0.018
  yTop <- min(yTop, 0.95)
  if (is.null(clust))
    xRight <- 1 - (strwidth(nms[length(nms)], units='figure', cex=input$nameSize) * 0.9 * sin(pi/180 * (90-input$nameAngle)))
 
  mx1 <- max(sapply(as.character(yvar), function(x) strwidth(x, units='figure', cex=input$axisSize)))
  mx2 <- max(sapply(as.character(yvar), function(x) strheight(x, units='figure')))
     
  ylabPos <- mx1/mx2*plotRatio*0.6 + 2.5
  xLeft <- mx1 + 0.03 / plotRatio
  if (nchar(input$yLabel) > 1)
    xLeft = xLeft + strheight(input$yLabel, units='figure', cex=input$axisSize) + 0.02 / plotRatio
  
  if (names(dev.cur())[1]=="pdf" | fileFlag) {
    xLeft = xLeft + .02
  }
  poly.line <- NA 
  if (style$line)
    poly.line <- input$outlineCol

  if (!is.null(dim(yvar)))
    yvar <- yvar[, 1, drop=TRUE]
  yLabels <- NULL
  style$ytks <- NULL
  ylim <- NULL
  if (is.character(yvar)) {
    yLabels <- yvar
    yvar <- 1:length(yvar)
    style$ytks <- yvar
  } else {
     ylim <- range(yvar, na.rm=TRUE)
     if (!is.na(style$yMin) & is.na(style$yMax)) {
         ylim[1] <- style$yMin
         ylim[2] <- max(yvar, na.rm=TRUE)
     } else if (!is.na(style$yMax) & is.na(style$yMin)) {
         ylim[2] <- style$yMax
         ylim[1] <- min(yvar, na.rm=TRUE)
     } else if (!is.na(style$yMax) & !is.na(style$yMin)) {
         ylim[1] <- style$yMin
         ylim[2] <- style$yMax
     }
     if (!is.null(ylim) & !is.na(input$yInterval)) {
       style$ytks <- seq(ylim[1], ylim[2], by=input$yInterval)
     }
  }
  exagCol <- input$exagCol
  if (2 %in% input$exag) {
     exagCol="auto"
  }

  exagMult <- input$exagMult
  if (is.na(input$exagMult)) {
    exagMult <- 1.0
    style$exag <- FALSE
  }
  
  xNames <- colnames(d)
  
  if (input$italicise) {
    xNames <- as.expression(sapply(xNames, function(x) bquote(italic(.(x))) ))

  }
  
  retVal <- tryCatch(x <- strat.plot(d, yvar = yvar, y.rev=style$yrev, scale.percent=style$scalePC, 
             plot.bar=style$bar, plot.line=style$line, plot.poly=style$poly, plot.symb=style$symbol, 
             col.poly=groupColours, col.bar=style$colBar, lwd.bar=style$lwdBar, col.symb=input$symbCol, 
             col.poly.line=poly.line, col.line=input$outlineCol, symb.cex=input$symbSize, exag=style$exag, 
             wa.order=style$autoOrder, bar.back=!input$barTop, clust=clust, cex.xlabel=input$nameSize, 
             srt.xlabel=input$nameAngle, yTop=yTop, xRight=xRight, ylabel=input$yLabel, 
             cex.yaxis=input$axisSize, cex.axis=0.8*input$axisSize, ylabPos=ylabPos,
             cex.ylabel=input$axisSize, tcl=-.4, mgp=c(3, input$axisSize/3, 0.3), xLeft=xLeft, 
             scale.minmax=style$scaleMinMax, ylim=ylim, y.tks=style$ytks, y.tks.labels=yLabels, 
             col.bg=NULL, col.exag=exagCol, exag.mult=exagMult, exag.alpha=0.15, x.names=xNames),
        error=function(e) return(e))
  if (inherits(retVal, "error")) {
    removeMsg()
    addNotification(retVal$message, duration=30, type="error", closeButton=TRUE)
  } 
  if (showZones > 1) {
     addClustZone(x, clust, showZones, col="red", yaxs="i")
  }
}
