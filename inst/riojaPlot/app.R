##
## Copyright (c) 2020, Steve Juggins
##
## License GPL-2
##
## Permission is hereby granted to use, copy, modify and distribute the software in accordance with
## the GPL-2 license and subject to the following condition:
##
## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
## WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(rioja))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(graphics))
suppressPackageStartupMessages(library(markdown))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(rhandsontable))

email <- tags$html( tags$body( a(href="mailto:Stephen.Juggins@ncl.ac.uk")))
riojaURL <- a("rioja", href="https://cran.r-project.org/web/packages/rioja/index.html")
SJURL <- a("Steve Juggins", href="mailto:Stephen.Juggins@ncl.ac.uk")
inputWidth <- "300px"
numericWidth <- "100px"
version <- "1.0"
colWidth <- "120px"
colList <- c(list('black' ),
    as.list(brewer.pal(9, "Blues")[-(1:3)]),
    as.list(brewer.pal(9, "Greens")[-(1:3)]),
    as.list(brewer.pal(9, "Reds")[-(1:3)]),
    as.list(brewer.pal(9, "Greys"))
  )

errorMsg <- ""
errorMsg2 <- ""
# errorMsg3 <- ""
currentSheet <- ""
currentFile <- ""
mydata <- NULL
fn <- ""
fn2 <- ""
DF <- NULL

plotIt <- function(fName, sheet, input, session, fileFlag=FALSE) {
#  sheet <- input$sheet
  
  errorMsg <<- ""
  if (is.character(sheet) & nchar(sheet)==0)
    return()
  selTaxa <- input$selTaxa
  yvar <- input$yvar
  style <- list()
  style$scalePC <- FALSE
  shinyjs::hideElement(id="errorBox2")
  if (currentSheet != sheet || currentFile != fName) {

    updateNumericInput(session, 'yMin', value=NA)
    updateNumericInput(session, 'yMax', value=NA)
    updateNumericInput(session, 'yInterval', value=NA)
    updateTextInput(session, 'yLabel', value="")

    sel <- input$misc
    sel <- setdiff(sel, c(2, 3, 4))
    updateCheckboxGroupInput(session, 'misc', selected=sel)
    updateCheckboxInput(session, 'doClust', value=FALSE)

    print(paste0("File: ", fName))
    print(paste0("Sheet: ", sheet))
    print(paste0("CurrentFile: ", currentFile))
    print(paste0("Currentsheet: ", currentSheet))
    
    currentSheet <<- sheet
    currentFile <<- fName
    d <- read_excel(fName, sheet=sheet)
    
# check column types:
    cols <- sapply(d, mode)
    nonChar = which(cols == "character")
    if (length(nonChar) > 0) {
      msg <-  paste0("Spreadsheet has ", length(nonChar), " non-numeric column(s). They have been removed.")
      showNotification(msg, duration=30, type="error", closeButton=TRUE)
      d <- d[, -nonChar]
    }
# check column headings
    colNames <- sum(!is.na(as.numeric(colnames(d))))
    if (colNames > 0) {
      msg <-  paste0("Spreadsheet has ", colNames, " column names that are numbers. 
                     Are you sure your data is formatted with variables in columns and observations in rows?")
      showNotification(msg, duration=30, type="error", closeButton=TRUE)
    }
    
    nMiss <- sum(is.na(d))
    if (nMiss > 0) {
      msg <- paste0("<p>Spreadsheet has ", nMiss, " missing value(s). ",
                     "Plotting silhouette diagrams and zonation is disabled.",  
                     "Check the help and example data for the correct data format.</p>")
      errorMsg0$mes1 <<- msg
      shinyjs::showElement(id="errorBox")
      shinyjs::disable(selector = "#showZones")
      shinyjs::disable(selector = "#nZones")
      shinyjs::disable(selector = "#doClust")
      shinyjs::disable(selector = "#style input[value='3']")
      updateCheckboxGroupInput(session, 'style', selected=as.numeric(1))
#        return("") 
    } else {
      shinyjs::hideElement(id="errorBox")
      errorMsg0$mes1 <<- ""
      shinyjs::enable(selector = "#showZones")
      shinyjs::enable(selector = "#nZones")
      shinyjs::enable(selector = "#doClust")
      shinyjs::enable(selector = "#style input[value='3']")
      updateCheckboxGroupInput(session, 'style', selected=as.numeric(3))
    }
    errorMsg2 <<- errorMsg0$mes1

    chron <- c("DEPTH", "AGE")
    sel1 <- toupper(substring(colnames(d), 1, 5)) %in% "DEPTH"
    sel2 <- toupper(substring(colnames(d), 1, 3)) %in% "AGE"
    sel <- sel1 | sel2
    if (any(sel)) {
      mydata$chron <<- d[, sel, drop=FALSE]
    } else {
      mydata$chron <<- data.frame(SampleNo=1:nrow(d))
    }
    mydata$spec <<- d[, !sel]
    yvars <- colnames(mydata$chron)
    yvar <- yvars[1]
    updateSelectInput(session, 'yvar', choices=yvars)
    
    selTaxa <- colnames(mydata$spec)
    mx <- apply(mydata$spec, 1, sum, na.rm=TRUE)
    r <- range(mx, na.rm=TRUE)
    sel <- input$misc
    if (r[1] > 50 & r[2] < 250) {
      sel <- apply(mydata$spec, 2, max, na.rm=TRUE) > 2
      selTaxa <- selTaxa[sel]
      style$scalePC <- TRUE
      sel <- input$misc
      sel <- unique(c(sel, "3"))
      updateCheckboxGroupInput(session, 'misc', selected=as.numeric(sel))
    } else {
        sel <- setdiff(sel, 3)      
        updateCheckboxGroupInput(session, 'misc', selected=as.numeric(sel))
        style$scalePC <- FALSE
    }
    sel2 <- colnames(mydata$spec) %in% selTaxa
    group <- "Group 1" # ifelse(sel2, "Group 1", NA)
    nGroups <- 6
    Groups <- factor(group, levels=paste0("Group ", 1:nGroups))
    DF <<- data.frame(Name=colnames(mydata$spec), Group="1", Selected=sel2, stringsAsFactors = FALSE)
    values$vDF <<- DF
    updatePickerInput(session, 'selTaxa', choices=colnames(mydata$spec), selected=selTaxa)
  }

  style$scaleMinMax <- FALSE
  style$exag <- FALSE
  style$yrev <- FALSE
  style$autoOrder <- "none"
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
  
  if (1 %in% input$misc)
    style$yrev <- TRUE
  if (2 %in% input$misc)
    style$exag <- TRUE
  if (3 %in% input$misc)
    style$scalePC <- TRUE
  if (4 %in% input$misc)
    style$scaleMinMax <- TRUE
  if (5 %in% input$misc)
    style$autoOrder <- "bottomleft"

  style$yMin <- input$yMin
  style$yMax <- input$yMax
  style$yInterval <- input$yInterval

  doZone <- input$doClust
    
  style$lwdBar <- input$barSize;
  style$colBar <- input$hLineCol
  style$SymbSize <- input$symbSize

  clust <- NULL
  showZones <- 0
  
   d <- mydata$spec
  
   if (!is.null(selTaxa) & length(selTaxa) > 0) {
     d <- d[, selTaxa]
   }

   nCol <- ncol(d)
    if (nCol > 50) {
       far <- ""
       if (nCol > 100) {
          errorMsg0$mes1 <<- paste0(errorMsg2, "<p>Your data has ", nCol, " variables. 
                                This is ", far, "too many to plot (the curves will be too small to visualise).  
                                All variables have been de-selected.  Use the selector under the 'Variables' tab 
                                to add variables to the diagram.</text></p>")
          updatePickerInput(session, 'selTaxa', choices=colnames(mydata$spec), selected=character(0))
          selTaxa <- NULL
          warning("too many to plot")
          return()
       } else {
         errorMsg0$mes1 <<- paste0(errorMsg2, "<p>Your data has ", nCol, " variables. 
                                This is too many to plot (the curves will be too small to visualise).  
                                Use the selector under the 'Variables' tab to de-select some and try again.  
                                80 is the about the maximum than can be plotted.  More than 50 makes no sense.</text></p>")
       }
   } else {
      errorMsg0$mes1 <<- errorMsg2
   }
 if (!is.null(yvar) & (nchar(yvar) > 1)) {
   yvar <- mydata$chron[, yvar, drop=TRUE]
 }
 else {
    yvar <- mydata$chron[, 1, drop=TRUE]
 }
 if (doZone) {
    diss <- dist(sqrt(d))
    clust <- chclust(diss)
    if (input$showZones == "No")
       showZones <- 0
    else if (input$showZones == "Auto") {
      x <- bstick(clust)
      x2 <- x$dispersion <= x$bstick
      showZones <- which(x2)[1]
      if (showZones < 2) {
        showNotification("There are no significant zones in these data.", duration=30, type="warning", closeButton=TRUE)
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
  yTop <- 1.0 - (yMax * plotRatio * cos(pi/180 * (90-input$nameAngle)))
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

  ylim <- NULL
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
  style$ytks <- NULL
  if (!is.null(ylim) & !is.na(input$yInterval)) {
    style$ytks <- seq(ylim[1], ylim[2], by=input$yInterval)
  }
  retVal <- tryCatch(x <- strat.plot(d, yvar = yvar, y.rev=style$yrev, scale.percent=style$scalePC, 
             plot.bar=style$bar, plot.line=style$line, plot.poly=style$poly, plot.symb=style$symbol, 
             col.poly=groupColours, col.bar=style$colBar, lwd.bar=style$lwdBar, col.symb=input$symbCol, 
             col.poly.line=poly.line, col.line=input$outlineCol, symb.cex=input$symbSize, exag=style$exag, 
             wa.order=style$autoOrder, bar.back=!input$barTop, clust=clust, cex.xlabel=input$nameSize, 
             srt.xlabel=input$nameAngle, yTop=yTop, xRight=xRight, ylabel=input$yLabel, 
             cex.yaxis=input$axisSize, cex.axis=0.8*input$axisSize, ylabPos=ylabPos,
             cex.ylabel=input$axisSize, tcl=-.4, mgp=c(3, input$axisSize/3, 0), xLeft=xLeft, 
             scale.minmax=style$scaleMinMax, ylim=ylim, y.tks=style$ytks, col.bg=NULL),
        error=function(e) return(e))
  if (inherits(retVal, "error")) {
    errorMsg0$mes3 <<- retVal$message
    shinyjs::showElement(id="errorBox2")
    warning(retVal$message)
  } else {
    errorMsg0$mes3 <<- ""
  }
  if (showZones > 1) {
     addClustZone(x, clust, showZones, col="red")
  }
}

"riojaPlot @ Newcastle University"

header <- dashboardHeader(title="riojaPlot", tags$li(a(href='riojaPlot.pdf', icon('question-circle-o'), title='Help', width=8, height=85), class='dropdown'))

D_ui <- dashboardPage(header, 
        dashboardSidebar(disable = TRUE),                       
        dashboardBody(
        tags$head(tags$style(HTML('
                           .skin-blue .main-header .logo { background-color: #3c8dbc; font-size: 2em; text-align: left; }
                              .skin-blue .main-header .logo:hover { background-color: #3c8dbc; }
                              .box {-webkit-box-shadow: none; -moz-box-shadow: none; box-shadow: none; }
                              .dropdown { padding: 0px; margin: 0px; font-size: 2em; border: 0; height: 50px; }
                              .errorMsg { color: red; padding-top: -20px; padding-bottom: -20; margin-top: -10px; 
                              margin-bottom: -20px; line-height: 1em; }
                              .main-header {max-height: 50px}
                              .main-header .logo {height: 50px;}
                              .sidebar-toggle {height: 50px; padding-top: 1px !important;}
                              .navbar {min-height:50px !important}
                              .error { color: red; padding: 0; margin-top: -20px;}
                      '))),
    # Boxes need to be put in a row (or column)
    fluidRow(shinyjs::useShinyjs(),
      column(width=2,
             box(fileInput("fn", "Upload Excel spreadsheet:", 
              accept=c("application/vnd.ms-excel",
                       "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        ".xlsx", ".xls")),
              div(style = "margin-top: -20px"),
              checkboxInput('exampleDF', 'Or use example data', value=FALSE), 
              div(style = "margin-top: -10px;"),
              width=NULL, solidHeader=TRUE),
        box(selectInput("sheet", "Select worksheet:", ""), width=NULL, solidHeader=TRUE),
        box(downloadButton("downloadResults", "Save plot"), 
            radioButtons('saveType', "", choices=c('pdf', 'png', 'svg'), selected='png', inline=TRUE), width=NULL, solidHeader=TRUE),
        box(p(tags$b(paste0("riojaPlot "), version)),
          paste0("Powered by "), 
              riojaURL, 
            p(paste0(" Version ", utils::packageDescription("rioja", fields="Version"), 
                  " (", utils::packageDescription("rioja", fields="Date"), ")")),
            "Comments & bug reports to ", SJURL, ".", width=NULL, solidHeader=TRUE)
      ),
      column(width=10,
        shinyjs::hidden(wellPanel(id="errorBox", htmlOutput("errorText", class="errorMsg"), width=NULL, solidHeader=TRUE)),
        shinyjs::hidden(wellPanel(id="errorBox2", htmlOutput("errorText2", class="errorMsg"), width=NULL, solidHeader=TRUE)),
        box(plotOutput("myPlot", height="600px"), width=NULL, solidHeader=TRUE),

        fluidRow(
          tabBox(title='',id='tab1', width=8, height="250px",
            tabPanel('Variables',      
              box( 
               selectInput('yvar', 'Y axis', choices='', selected='', width=inputWidth, multiple=FALSE),
              solidHeader=TRUE, width=4),
              box( 
                pickerInput('selTaxa', 'Select X variables', choices='', multiple=TRUE, options=pickerOptions(dropupAuto=FALSE)),
              solidHeader=TRUE, width=4),
            ),
            tabPanel(title='Settings', height="250px",
               box( 
                  checkboxGroupInput('style', 'Style', choices = c(Line=1, Symbols=2, Silhouette=3), selected=3),
               solidHeader=TRUE, width=3),
               box( 
                    radioButtons("hLine", "Show bar", choices=c("None", "Curve", "Full"), selected="Curve", inline=FALSE), 
                    checkboxInput('barTop', 'Bars on top', value=TRUE),
               solidHeader=TRUE, width=3),
               box(
                 numericInput('symbSize', 'Symbol size', value=0.5, min=0.2, max=4, step=0.1, width=numericWidth),
                 numericInput('barSize', 'Bar width', value=1, min=1, max=10, step=1, width=numericWidth), 
                 solidHeader=TRUE, width=3),
               box( 
                  checkboxGroupInput('misc', 'Settings', c('Reverse Y axis'=1, 'Show 5x exag'=2, 'Scale for %'=3, 
                                                           'Show min/max'=4, 'Auto sort vars'=5), 
                                                   selected=c(1, 2)),
               solidHeader=TRUE, width=3)
            ),
            tabPanel(title='Y axis', height="250px",
               box( 
                 numericInput('yMin', 'Y axis min', value=NA, width=numericWidth),
                 numericInput('yMax', 'Y axis max', value=NA, width=numericWidth),
              solidHeader=TRUE, width=3),
              box( 
                numericInput('yInterval', 'Y axis interval', value=NA, width=numericWidth),
                textInput('yLabel', 'Label for Y axis'),
                solidHeader=TRUE, width=3)
            ),
          tabPanel(title='Colours', 
            box(
              spectrumInput("outlineCol", "Line", "black", width=colWidth, choices = colList, update_on='dragstop'),
              spectrumInput("fillCol", "Silhouette", "darkgreen", width=colWidth, choices = colList, update_on='dragstop'),
              solidHeader=TRUE, width=2),
            box(
              spectrumInput("hLineCol", "Bar", "#D9D9D9", width=colWidth, choices = colList, update_on='dragstop'), 
             spectrumInput("zoneCol", "Zone", "red", width=colWidth, choices = colList, update_on='dragstop'),
           solidHeader=TRUE, width=2),
           box(
             spectrumInput("symbCol", "Symbol", "#000000", width=colWidth, choices=colList, update_on='dragstop'),
             solidHeader=TRUE, width=2),
          ),
          tabPanel(title='Sizes',  
            box(
              numericInput('nameSize', 'X-var font size', value=1, min=.4, max=2, step=0.05, width=numericWidth), 
              numericInput('nameAngle', 'Rotate names', value=90, min=0, max=90, step=5, width=numericWidth), 
            solidHeader=TRUE, width=3),
            box(
              numericInput('axisSize', 'Axis font size', value=1, min=.4, max=2, step=0.05, width=numericWidth), 
              solidHeader=TRUE, width=3)
          ),
          tabPanel(title='Zonation', 
            box(
             checkboxInput('doClust', 'Add zonation', value=FALSE),
             radioButtons("showZones", "Show zones", choices=c("No", "Auto", "Choose"), selected="Auto", inline=FALSE), 
            solidHeader=TRUE, width=3),
            box(
             numericInput('nZones', 'Number of zones', value=2, min=2, max=10, step=1, width="120px"),
            solidHeader=TRUE, width=4),
          ),
          tabPanel(title='Groups', 
                   box(
                     spectrumInput("groupCol1", "Group 1", "darkgreen", width=colWidth, choices = colList, update_on='dragstop'),
                     spectrumInput("groupCol2", "Group 2", "darkkhaki", width=colWidth, choices = colList, update_on='dragstop'),
                     solidHeader=TRUE, width=2),
                   box(
                     spectrumInput("groupCol3", "Group 3", "darkorange", width=colWidth, choices = colList, update_on='dragstop'),
                     spectrumInput("groupCol4", "Group 4", "darkred", width=colWidth, choices = colList, update_on='dragstop'),
                     solidHeader=TRUE, width=2),
                   box(
                     spectrumInput("groupCol5", "Group 5", "deepskyblue", width=colWidth, choices = colList, update_on='dragstop'),
                     spectrumInput("groupCol6", "Group 6", "firebrick", width=colWidth, choices = colList, update_on='dragstop'),
                     solidHeader=TRUE, width=2),
                   box(rHandsontableOutput("varList"),
                       solidHeader=TRUE, width=6)
                   ),
          tabPanel(title='Help', 
                   includeMarkdown("riojaHelp.md"),
                   solidHeader=TRUE, width=12)
          )                   
        )
      )
    )
  )
)

summarise_data <- function(fn, sheet, data) {
  nsam <- NULL; nsp <- NULL;
   if (!is.null(data)) {
      nsam <- nrow(data)
      nsp <- ncol(data)
   }
  paste("File name: ", fn, "\n\rSheet:", sheet, "\n\nNumber of samples: ", nsam, "\nNumber of taxa: ", nsp )
}

values <- reactiveValues(vDF=DF)
errorMsg0 <- reactiveValues(mes1="", mes3="")

D_server <- function(input, output, session) {
   shinyjs::disable("downloadResults")
   shinyjs::hide(id="errorBox")
   errorMsg <<- ""
   
   observeEvent(input$selTaxa, {
     if (!is.null(input$selTaxa) & !is.null(input$varList)) {
        DF <<- hot_to_r(input$varList)
        DF$Selected <<- colnames(mydata$spec) %in% input$selTaxa
        isolate(set_data('varList', 1:length(DF$Selected), 3, DF$Selected, session))
        values$vDF <<- DF
     }
   }, ignoreInit=TRUE)

   output$varList = renderRHandsontable({
     if (!is.null(values$vDF)) {
       sel <- which(colnames(mydata$spec) %in% input$selTaxa) - 1
#       newG <- values$vDF
#       sel <- which(newG$Selected) - 1 
       rht <- rhandsontable(values$vDF, selected=sel, readOnly=TRUE, height=160, width=220, contextMenu = FALSE, 
                   rowHeaders=FALSE) %>%              
       hot_validate_numeric("Group", min=1, max=6, choices=c("1", "2", "3", "4", "5", "6")) %>%
       hot_col("Group", readOnly=FALSE, format="0") %>%
       hot_col("Selected", width=0.1)  %>%
       hot_cols(renderer = " function (instance, td, row, col, prop, value, cellProperties) {
          Handsontable.renderers.TextRenderer.apply(this, arguments);
          if (instance.params) {
               hrows = instance.params.selected
               hrows = hrows instanceof Array ? hrows : [hrows]              
               if (instance.params && hrows.includes(row)) {
                 td.style.color = 'black'; 
               } else {
                 td.style.color = 'lightgrey';
               }
            }
       }"
       )
      }
   }) 
   
   observeEvent(input$fillCol, {
      updateSpectrumInput(session, 'groupCol1', selected=input$fillCol)
    }, ignoreInit=TRUE)

   observeEvent(input$groupCol1, {
     updateSpectrumInput(session, 'fillCol', selected=input$groupCol1)
   }, ignoreInit=TRUE)
   
   
   observeEvent(input$exampleDF, {
    output$myPlot <- renderPlot({
      if (input$exampleDF) {
        shinyjs::hideElement(id="errorBox")
        shinyjs::hideElement(id="errorBox2")
        fn <- system.file("riojaPlot/www/aber.xlsx", package="rioja")
        plotIt(fn, "Aber", input, session)
        shinyjs::enable("downloadResults")
      } else {
        if (nchar(fn2) > 5) {
           plotIt(fn2, input$sheet, input, session)
           shinyjs::enable("downloadResults")
        } else {
          updateSelectInput(session, 'yvar', choices=integer(0))
          updatePickerInput(session, 'selTaxa', choices=integer(0))
          mydata <<- NULL
          fn2 <- ""
          currentFile <<- ""
          currentSheet <<- ""
          values$vDF <<- NULL
        }         
      }
    })
    
  }, ignoreNULL = TRUE)

  observeEvent(errorMsg0, { 
    if (nchar(errorMsg0$mes3) > 1) {
       output$errorText2 <- renderText(errorMsg0$mes3)
       shinyjs::showElement(id="errorBox2")
    } else {
      shinyjs::hideElement(id="errorBox2")
    }      
    if (nchar(errorMsg0$mes1) > 1) {
      output$errorText <- renderText(errorMsg0$mes)
      shinyjs::showElement(id="errorBox")
      } else {
        shinyjs::hideElement(id="errorBox")
      }      
    })

#  observeEvent(input$errorText, { 
#    output$errorText <- renderText({
#      if (nchar(errorMsg) > 1) {
#        shinyjs::showElement(id="errorBox")
#        return(errorMsg)
#      } else {
#        shinyjs::hideElement(id="errorBox")
#      }      
#    })
#  })
  
   
  output$errorText2 <- renderText({
    if (nchar(errorMsg0$mes3) > 1) {
        shinyjs::showElement(id="errorBox2")
        return(errorMsg0$mes3)
    } else {
        shinyjs::hideElement(id="errorBox2")
    }      
  })  
  output$errorText <- renderText({
    if (nchar(errorMsg0$mes1) > 1) {
       shinyjs::showElement(id="errorBox")
       return(errorMsg0$mes1)
  } else {
       shinyjs::hideElement(id="errorBox")
  }      
  })  

  observeEvent(input$fn$name, {
    get_Sheets <- function(file) {
      sheets <- readxl::excel_sheets(file)
    }
    fn1 <- input$fn$name
    cat(paste("Input file = ", input$fn$name, "\n"))
    cat(paste("File type =", input$fn$type, "\n"))
    cat(paste("Local path =", input$fn$datapath, "\n"))
    if (is.null(fn1)) {
      return(NULL)
    }
    dn <- dirname(input$fn$datapath)
    fn2 <<- file.path(dn, input$fn$name)
    if (file.exists(fn2))
       file.remove(fn2)
    file.rename(input$fn$datapath, fn2)
    fn <<- fn1
    mydata <<- NULL
    sheets.nms <<- tryCatch(get_Sheets(fn2), error=function(e) return (e))
    if (inherits(sheets.nms, "error")) {
      updateSelectInput(session, "sheet", choices="")
      shinyjs::showElement(id="errorBox")
      output$errorText <- renderText(paste0("Cannot open Excel file.\nReason: ", sheets.nms$message), quoted=TRUE)
      fn1 <- input$fn$name
      sheets.nms <<- ""
      shinyjs::disable("downloadResults")
      return()
    } else {
      shinyjs::enable("downloadResults")
      updateSelectInput(session, "sheet", choices=sheets.nms, selected=sheets.nms[1])
      if (nchar(input$sheet > 0)) {
         output$myPlot <- renderPlot({
             plotIt(fn2, input$sheet, input, session)
         })
      }
    }
  }, ignoreNULL = TRUE)

  output$downloadResults <- downloadHandler(
    filename = function() {
      tmp <- basename(fn2)
      tmp <- strsplit(tmp, "\\.")[[1]][1]
      ext <- ".png"
      if (input$saveType=='pdf') {
        ext <- ".pdf"
      } else if (input$saveType=='svg') {
        ext <- ".svg"
      }
      outFile <- paste0("MyPlot_", tmp, "_", Sys.Date(), ext)
      outFile <- gsub(" ", "_", outFile)
      outFile
    },
    content <- function(file) {
      flagIt <- FALSE
      ratio <-  session$clientData$output_myPlot_width / session$clientData$output_myPlot_height
      if (input$saveType=='png') {
        png(file, 
          width = session$clientData$output_myPlot_width,
          height = session$clientData$output_myPlot_height
        )
      } else if (input$saveType=='pdf') {
        ratio <-  session$clientData$output_myPlot_width / session$clientData$output_myPlot_height
          width <- 15
          pdf(file, 
            width =  width,
            height = width / ratio
         )
      } else {
        ratio <-  session$clientData$output_myPlot_width / session$clientData$output_myPlot_height
        width <- 15
        flagIt <- TRUE
        CairoSVG(file, 
           width =  width,
           height = width / ratio
        )
      } 
      plotIt(currentFile, currentSheet, input, session, flagIt)
      dev.off()
    } 
  )
  
#  session$onSessionEnded(stopApp) # kill shiny on browser close
}

shinyApp(D_ui, D_server)
