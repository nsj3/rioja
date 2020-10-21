##
## Copyright (c) 2019, Steve Juggins
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

email <- tags$html( tags$body( a(href="mailto:Stephen.Juggins@ncl.ac.uk")))
riojaURL <- a("rioja", href="https://cran.r-project.org/web/packages/rioja/index.html")
SJURL <- a("Steve Juggins", href="mailto:Stephen.Juggins@ncl.ac.uk")
inputWidth <- "300px"
numericWidth <- "100px"
version <- "1.0"
colWidth <- "120px"

errorMsg <- ""
currentSheet <- ""
currentFile <- ""
mydata <- NULL
fn <- ""
fn2 <- ""

plotIt <- function(fName, sheet, input, session, fileFlag=FALSE) {
#  sheet <- input$sheet
  errorMsg <<- ""
  if (is.character(sheet) & nchar(sheet)==0)
    sheet = 1
  selTaxa <- input$selTaxa
  yvar <- input$yvar
  style <- list()
  style$scalePC <- FALSE
  if (3 %in% input$misc)
     TRUE  
  if (currentSheet != sheet || currentFile != fName) {
    d <- read_excel(fName, sheet=sheet)
    nMiss <- sum(is.na(d))
    if (nMiss > 0) {
      msg <- paste0("<p><text style='color:red'>Spreadsheet has ", nMiss, " missing values. ",
                     "Plotting silouette diagrams and zonation are disabled.<br>",  
                     "Check the help and example data for the correct data format.</text></p>")
      errorMsg <<- msg
      shinyjs::showElement(id="errorBox")
      shinyjs::disable(selector = "#showZones")
      shinyjs::disable(selector = "#nZones")
      shinyjs::disable(selector = "#doClust")
      shinyjs::disable(selector = "#style input[value='3']")
      updateCheckboxGroupInput(session, 'style', selected=as.numeric(1))
#        return("") 
    } else {
      shinyjs::hideElement(id="errorBox")
      errorMsg <<- ""
      shinyjs::enable(selector = "#showZones")
      shinyjs::enable(selector = "#nZones")
      shinyjs::enable(selector = "#doClust")
      shinyjs::enable(selector = "#style input[value='3']")
      updateCheckboxGroupInput(session, 'style', selected=as.numeric(3))
    }

    currentSheet <<- sheet
    currentFile <<- fName
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
    if (any(mx > 50 && mx < 150)) {
      sel <- apply(mydata$spec, 2, max, na.rm=TRUE) > 2
      selTaxa <- sel
      style$scalePC <- TRUE
      sel <- input$misc
      sel <- unique(c(sel, "3"))
      updateCheckboxGroupInput(session, 'misc', selected=as.numeric(sel))
    } 
    updatePickerInput(session, 'selTaxa', choices=colnames(mydata$spec), selected=colnames(mydata$spec)[selTaxa])
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
    } else {
       showZones = input$nZones
    }
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
  if (nchar(input$title) > 1)
    xLeft = xLeft + strheight(input$title, units='figure', cex=input$axisSize) + 0.02 / plotRatio
  
  if (names(dev.cur())[1]=="pdf" | fileFlag) {
    xLeft = xLeft + .02
  }
  poly.line <- NA 
  if (style$line)
    poly.line <- input$outlineCol

  x <- strat.plot(d, yvar = yvar, y.rev=style$yrev, scale.percent=style$scalePC, 
             plot.bar=style$bar, plot.line=style$line, plot.poly=style$poly, plot.symb=style$symbol, 
             col.poly=input$fillCol, col.bar=style$colBar, lwd.bar=style$lwdBar, col.symb=input$symbCol, 
             col.poly.line=poly.line, col.line=input$outlineCol, symb.cex=input$symbSize, exag=style$exag, 
             wa.order=style$autoOrder, bar.back=!input$barTop, clust=clust, cex.xlabel=input$nameSize, srt.xlabel=input$nameAngle,
             yTop=yTop, xRight=xRight, ylabel=input$title, cex.yaxis=input$axisSize, cex.axis=0.8*input$axisSize, ylabPos=ylabPos, 
             cex.ylabel=input$axisSize, tcl=-.4, mgp=c(3, input$axisSize/3, 0), xLeft=xLeft, scale.minmax=style$scaleMinMax)   
  if (showZones > 0) {
     addClustZone(x, clust, showZones, col="red")
  }

}

"riojaPlot @ Newcastle University"

header <- dashboardHeader(title="riojaPlot", tags$li(a(href='riojaPlot.pdf', icon('info-circle'), title='Help', width=8), class='dropdown'))

D_ui <- dashboardPage(header, 
        dashboardSidebar(disable = TRUE),                       
        dashboardBody(
        tags$head(tags$style(HTML('.skin-blue .main-header .logo {
                              background-color: #3c8dbc;
                              font-size: 2em;
                              text-align: left; }
                              .skin-blue .main-header .logo:hover {
                              background-color: #3c8dbc; }
                              .box {-webkit-box-shadow: none; 
                              -moz-box-shadow: none;
                              box-shadow: none; }
                              .dropdown {
                              padding: 0px;
                              margin: 0px;
                              font-size: 1.8em;
                              border: 0;
                              height: 50px;
                              }
                              .errorMsg {
                              padding-top: 0px;
                              margin-top: -15px;
                              margin-bottom: -5px;
                              height: 30px;
                              }
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
        shinyjs::hidden(wellPanel(id="errorBox", htmlOutput("errorText", class="errorMsg", height="30px"), width=NULL, solidHeader=TRUE)),
        box(plotOutput("myPlot", height="600px"), width=NULL, solidHeader=TRUE),

        fluidRow(
          tabBox(title='',id='tab1', width=8, height="250px",
            tabPanel('Select Variables',      
              box( 
               selectInput('yvar', 'Y axis', choices='', selected='', width=inputWidth, multiple=FALSE),
               textInput('title', 'Label for Y axis'),
              solidHeader=TRUE, width=4),
              box( 
                pickerInput('selTaxa', 'Select X vars', choices='', multiple=TRUE, options=pickerOptions(dropupAuto=FALSE)),
              solidHeader=TRUE, width=4),
            ),
            tabPanel(title='Settings', height="300px",
               box( 
                  checkboxGroupInput('style', 'Style', choices = c(Line=1, Symbols=2, Silouette=3), selected=3),
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
                                                   selected=c(1, 2, 3)),
               solidHeader=TRUE, width=3)
            ),
          tabPanel(title='Colours', 
            box(
              spectrumInput("outlineCol", "Line", "black", width=colWidth, choices = c(
                list('black' ),
                as.list(brewer.pal(9, "Blues")[-(1:3)]),
                as.list(brewer.pal(9, "Greens")[-(1:3)]),
                as.list(brewer.pal(9, "Reds")[-(1:3)]),
                as.list(brewer.pal(9, "Greys")))
              ),
              spectrumInput("fillCol", "Silouette", "darkgreen", width=colWidth, choices = c(
                list('black' ),
                as.list(brewer.pal(9, "Blues")[-(1:3)]),
                as.list(brewer.pal(9, "Greens")[-(1:3)]),
                as.list(brewer.pal(9, "Reds")[-(1:3)]),
                as.list(brewer.pal(9, "Greys")))
              ),
              solidHeader=TRUE, width=2),
            box(
              spectrumInput("hLineCol", "Bar", "#D9D9D9", width=colWidth, choices = c(
                 list('black' ),
                 as.list(brewer.pal(9, "Blues")[-(1:3)]),
                 as.list(brewer.pal(9, "Greens")[-(1:3)]),
                 as.list(brewer.pal(9, "Reds")[-(1:3)]),
                 as.list(brewer.pal(9, "Greys"))), 
             ), 
             spectrumInput("zoneCol", "Zone", "red", width=colWidth, choices = c(
               list('black' ),
               as.list(brewer.pal(9, "Blues")[-(1:3)]),
               as.list(brewer.pal(9, "Greens")[-(1:3)]),
               as.list(brewer.pal(9, "Reds")[-(1:3)]),
               as.list(brewer.pal(9, "Greys"))), 
             ),
           solidHeader=TRUE, width=2),
           box(
             spectrumInput("symbCol", "Symbol", "#000000", width=colWidth, choices = c(
               list('black' ),
               as.list(brewer.pal(9, "Blues")[-(1:3)]),
               as.list(brewer.pal(9, "Greens")[-(1:3)]),
               as.list(brewer.pal(9, "Reds")[-(1:3)]),
               as.list(brewer.pal(9, "Greys"))), 
             ), 
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

D_server <- function(input, output, session) {
   shinyjs::disable("downloadResults")
   shinyjs::hide(id="errorBox")
   errorMsg <<- ""
   observeEvent(input$exampleDF, {
    output$myPlot <- renderPlot({
      shinyjs::hideElement(id="errorBox")
      if (input$exampleDF) {
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
        }         
      }
    })
    
  }, ignoreNULL = TRUE)
  
  observeEvent(input$fn$name, {
    get_Sheets <- function(file) {
      sheets <- readxl::excel_sheets(file)
    }
    fn1 <- input$fn$name
    cat(paste("Input file = ", input$fn$name, "\n"))
    cat(paste("fn =", input$fn$type, "\n"))
    cat(paste("fn =", input$fn$datapath, "\n"))
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
      output$message1 <- renderText(paste0("Cannot open Excel file.\nReason: ", sheets.nms$message), quoted=TRUE)
      fn1 <- input$fn$name
      sheets.nms <<- ""
      shinyjs::disable("downloadResults")
      return()
    } else {
      shinyjs::enable("downloadResults")
    }
    updateSelectInput(session, "sheet", choices=sheets.nms)
    output$myPlot <- renderPlot({
      plotIt(fn2, input$sheet, input, session)
      if (nchar(errorMsg) > 1 ) {
#        shinyjs::showElement(id="errorBox")
        output$errorText <- renderText(errorMsg)
      } else {
#        shinyjs::hideElement(id="errorBox")
      }
    })
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
}

shinyApp(D_ui, D_server)
