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
suppressPackageStartupMessages(library(svglite))

header <- dashboardHeader(title = paste0("Rioja plotting"), titleWidth=750)
email <- tags$html( tags$body( a(href="mailto:Stephen.Juggins@ncl.ac.uk")))
riojaURL <- a("rioja", href="https://cran.r-project.org/web/packages/rioja/index.html")
SJURL <- a("Steve Juggins", href="mailto:Stephen.Juggins@ncl.ac.uk")

currentSheet <- ""
currentFile <- ""
colWidth <- "120px"
mydata <- NULL
currentEx <- FALSE

plotIt <- function(fName, sheet, input, session) {
#  sheet <- input$sheet
  print(input$sheet)
  print(fName)
  if (is.character(sheet) & nchar(sheet)==0)
    sheet = 1
  selTaxa <- input$selTaxa
  yvar <- input$yvar
  if (currentSheet != sheet || currentFile != fName) {
    d <- read_excel(fName, sheet=sheet)
    currentSheet <<- sheet
    currentFile <<- fName
    print("Read file")
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
    print(colnames(mydata$chron))
    updateSelectInput(session, 'yvar', choices=yvars, selected=yvars[1])
    yvar <- yvars[1]
    sel <- apply(mydata$spec, 2, max) > 2   
    selTaxa <- sel
    updatePickerInput(session, 'selTaxa', choices=colnames(mydata$spec), selected=colnames(mydata$spec)[sel])
  }
  
  style <- list()
  style$exag <- FALSE
  style$yrev <- FALSE
  style$scalePC <- FALSE
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
    style$autoOrder <- "bottomleft"

  doZone <- input$doClust
    
  style$lwdBar <- input$barSize;
  style$colBar <- input$hLineCol
  style$SymbSize <- input$symbSize
#  print(style$SymbSize)

  clust <- NULL
  showZones <- 0
  
  d <- mydata$spec
  
  d <- d[, selTaxa]
  yvar <- mydata$chron[, yvar, drop=TRUE]
    
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
  if (TRUE) {
     yMax <- max(sapply(nms, function(x) strwidth(x, units="figure", cex=input$nameSize))) 
     plotRatio<- session$clientData$output_myPlot_width / session$clientData$output_myPlot_height
     yTop <- 1.0 - (yMax * plotRatio * cos(pi/180 * (90-input$nameAngle)))
     yTop <- min(yTop, 0.95)
     if (is.null(clust))
        xRight <- 1 - (strwidth(nms[length(nms)], units='figure') * 0.9 * sin(pi/180 * (90-input$nameAngle)))
  }

  print ("yTop")
  print(yTop)
  print(xRight)
  
  x <- strat.plot(d, yvar = yvar, y.rev=style$yrev, scale.percent=style$scalePC, 
             plot.bar=style$bar, plot.line=style$line, plot.poly=style$poly, plot.symb=style$symbol, 
             col.poly=input$fillCol, col.bar=style$colBar, lwd.bar=style$lwdBar, col.symb=style$fillCol, 
             col.poly.line=input$outlineCol, col.line=input$outlineCol, symb.cex=style$symbSize, exag=style$exag, 
             wa.order=style$autoOrder, bar.back=!input$barTop, clust=clust, cex.xlabel=input$nameSize, srt.xlabel=input$nameAngle,
             yTop=yTop, xRight=xRight, ylabel=input$title)
  if (showZones > 0) {
     addClustZone(x, clust, showZones, col="red")
  }
}

D_ui <- dashboardPage(header, dashboardSidebar(disable = TRUE),
  dashboardBody(
    tags$head(tags$style(HTML('.skin-blue .main-header .logo {
                              background-color: #3c8dbc;
                              text-align: left; }
                              .skin-blue .main-header .logo:hover {
                              background-color: #3c8dbc; }
                              '))),
    # Boxes need to be put in a row (or column)
    fluidRow(shinyjs::useShinyjs(),
             
      column(width=3,
             box(fileInput("fn", "Select input file:",
              accept=c("application/vnd.ms-excel",
                       "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        ".xlsx", ".xls")),
              checkboxInput('exampleDF', 'Or use example data', value=FALSE), 
              width=NULL),
        box(textOutput("message1"),
            selectInput("sheet", "Select worksheet:", ""), width=NULL),
        box(downloadButton("downloadResults", "Save plot"), 
            radioButtons('saveType', "", choices=c('pdf', 'png', 'svg'), selected='png'), width=NULL),
        box(paste0("Powered by "), 
              riojaURL, 
            p(paste0(" Version ", utils::packageDescription("rioja", fields="Version"), 
                  " (", utils::packageDescription("rioja", fields="Date"), ")")),
            "Please email comments, bug reports etc to ", SJURL, ".", width=NULL)
      ),
      column(width=9,
        box(plotOutput("myPlot"), width=NULL, solidHeader=TRUE),
        wellPanel(
        fluidRow(
          column(3, 
            selectInput('yvar', 'Y axis', choices='', selected='', width="80%"),
            checkboxGroupInput('style', 'Style', choices = c(Line=1, Symbols=2, Shilouette=3), selected=3),
            radioButtons("hLine", "Show bar", choices=c("None", "Curve", "Full"), selected="Curve", inline=TRUE), 
            checkboxInput('barTop', 'Bars on top', value=TRUE),
            textInput('title', 'Y-axis title')
          ),
          column(3, 
            pickerInput('selTaxa', 'Select X vars', choices='', multiple=TRUE, options=pickerOptions(dropupAuto=FALSE), width="80%"),
            checkboxGroupInput('misc', 'Settings', c('Reverse Y axis'=1, 'Show 5x exag'=2, 'Scale for % data'=3, 'Auto order taxa'=4), 
                                                     selected=c(1, 2, 3)),
          ),
          column(2, 
            spectrumInput("outlineCol", "Line", "black", width=colWidth, choices = c(
                list('black' ),
                as.list(brewer.pal(9, "Blues")[-(1:3)]),
                as.list(brewer.pal(9, "Greens")[-(1:3)]),
                as.list(brewer.pal(9, "Reds")[-(1:3)]),
                as.list(brewer.pal(9, "Greys")))
             ),
             spectrumInput("fillCol", "Shilouette", "darkgreen", width=colWidth, choices = c(
                list('black' ),
                as.list(brewer.pal(9, "Blues")[-(1:3)]),
                as.list(brewer.pal(9, "Greens")[-(1:3)]),
                as.list(brewer.pal(9, "Reds")[-(1:3)]),
                as.list(brewer.pal(9, "Greys")))
             ),
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
             )
          ),
          column(2, 
             numericInput('symbSize', 'Symbol size', value=1, min=0.2, max=4, step=0.1, width="80%"),
             numericInput('barSize', 'Bar width', value=1, min=1, max=10, step=1, width="80%"), 
             numericInput('nameSize', 'Font size', value=1, min=.4, max=2, step=0.05, width="80%"), 
             numericInput('nameAngle', 'Rotate names', value=90, min=0, max=90, step=5, width="80%"), 
          ),
          column(2, 
             checkboxInput('doClust', 'Add zonation', value=FALSE),
             radioButtons("showZones", "Show zones", choices=c("No", "Auto", "Choose"), selected="Auto", inline=FALSE), 
             numericInput('nZones', 'Number of zones', value=2, min=2, max=10, step=1, width="80%"),
          )
        )
      )
    )
  )
)
)

fn <- ""
fn2 <- ""
mydata <- NULL
id <- NULL

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
  
  observeEvent(input$exampleDF, {
    output$myPlot <- renderPlot({
      if (input$exampleDF) {
        fn <- system.file("shiny_app/aber.xlsx", package="rioja")
        plotIt(fn, "Aber", input, session)
        shinyjs::enable("downloadResults")
      } else {
        if (nchar(fn2) > 5) {
           plotIt(fn2, input$sheet, input, session)
           shinyjs::enable("downloadResults")
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
    output$table2 <- renderText("")
    mydata <<- NULL
    sheets.nms <<- tryCatch(get_Sheets(fn2), error=function(e) return (e))
    if (inherits(sheets.nms, "error")) {
      updateSelectInput(session, "sheet", choices="")
      output$message1 <- renderText(paste0("Cannot open Excel file.\nReason: ", sheets.nms$message), quoted=TRUE)
      fn1 <- input$fn$name
      sheets.nms <<- ""
      output$table1 <- renderText("")
      shinyjs::disable("downloadResults")
      return()
    } else {
      shinyjs::enable("downloadResults")
    }
    updateSelectInput(session, "sheet", choices=sheets.nms)
    output$myPlot <- renderPlot({
      plotIt(fn2, input$sheet, input, session)
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
      print (session$clientData$output_myPlot_width)
      print(session$clientData$output_myPlot_height)
      ratio <-  session$clientData$output_myPlot_width / session$clientData$output_myPlot_height
      if (input$saveType=='png') {
        png(file, 
          width = session$clientData$output_myPlot_width,
          height = session$clientData$output_myPlot_height
        )
      } else if (input$saveType=='pdf') {
        ratio <-  session$clientData$output_myPlot_width / session$clientData$output_myPlot_height
          width <- 10
          pdf(file, 
            width =  width,
            height = width / ratio
         )
      } else {
        ratio <-  session$clientData$output_myPlot_width / session$clientData$output_myPlot_height
        width <- 10
        svglite(file, 
           width =  width,
           height = width / ratio
        )
      } 
      plotIt(currentFile, currentSheet, input, session)
      dev.off()
    } 
  )
}

shinyApp(D_ui, D_server)
