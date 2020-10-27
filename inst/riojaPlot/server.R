
D_server <- function(input, output, session) {
  shinyjs::disable("downloadResults")
  shinyjs::hide(id="errorBox")
  errorMsg <<- ""
  
  observeEvent(input$scalePC, {
    if (input$scalePC) {
      shinyjs::enable("autoSelect")
      shinyjs::enable("minCut")
      shinyjs::enable("maxCut")
    } else {
      shinyjs::disable("autoSelect")
      shinyjs::disable("minCut")
      shinyjs::disable("maxCut")
    }
  })
  
  observeEvent(input$autoSelect, {
    selectVars(session, input)
  })

  output$varList = renderRHandsontable({
    if (!is.null(input$selTaxa)) {
      sel <- which(colnames(mydata$spec) %in% input$selTaxa) - 1
      rht <- rhandsontable(mydata$group, selected=sel, readOnly=TRUE, height=160, width=320, contextMenu = FALSE, 
                           rowHeaders=FALSE) %>%              
        hot_validate_numeric("Group", min=1, max=6, choices=c("1", "2", "3", "4", "5", "6")) %>%
        hot_col("Group", readOnly=FALSE, format="0") %>%
        hot_col("Selected", width=0.1)  %>%
        hot_cols(colWidths=c(250, 50), renderer = " function (instance, td, row, col, prop, value, cellProperties) {
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
    isolate(updateSpectrumInput(session, 'groupCol1', selected=input$fillCol))
  }, ignoreInit=TRUE)
  
  observeEvent(input$groupCol1, {
    isolate(updateSpectrumInput(session, 'fillCol', selected=input$groupCol1))
  }, ignoreInit=TRUE)
  
  observeEvent(input$exampleDF, {
    output$myPlot <- renderPlot({
      if (input$exampleDF) {
        isolate(updateSelectInput(session, 'yvar', choices=integer(0)))
        isolate(updatePickerInput(session, 'selTaxa', choices=integer(0)))
        updateSelectInput(session, "sheet", choices=character(0))
        reset('fn')
        mydata$spec <<- NULL
        mydata$chron <<- NULL
        mydata$group <<- NULL
        fn2 <<- ""
        sheet <<- ""
        currentFile <<- ""
        currentSheet <<- ""
        removeMsg();
        fn <- system.file("riojaPlot/www/aber.xlsx", package="rioja")
        readIt(fn, "Aber", input, session)
        output$myPlot <- renderPlot({
          plotIt(mydata, input, session)
        })
        txt <- paste0("<a href='aber.xlsx'>Download</a> example data.")  
        output$exampleData <- renderUI({
          HTML(txt)})
        shinyjs::showElement("exampleData")
        shinyjs::enable("minCut")
        shinyjs::enable("maxCut")
        shinyjs::disable("sheet")
        shinyjs::enable("autoSelect")
        shinyjs::enable("downloadResults")
      } else {
        shinyjs::hide("exampleData")
        removeMsg();
        mydata$spec <<- NULL
        mydata$chron <<- NULL
        mydata$group <<- NULL
        fn2 <<- ""
        sheet <<- ""
        currentFile <<- ""
        currentSheet <<- ""
        isolate(updateSelectInput(session, 'yvar', choices=integer(0)))
        isolate(updatePickerInput(session, 'selTaxa', choices=integer(0)))
      }
    })
  }, ignoreNULL = TRUE)
  
  observeEvent(input$fn$name, {
    get_Sheets <- function(file) {
      sheets <- readxl::excel_sheets(file)
    }
    fn1 <- input$fn$name
#    cat(paste("Input file = ", input$fn$name, "\n"))
#    cat(paste("File type =", input$fn$type, "\n"))
#    cat(paste("Local path =", input$fn$datapath, "\n"))
    if (is.null(fn1)) {
      return(NULL)
    }
    dn <- dirname(input$fn$datapath)
    fn2 <<- file.path(dn, input$fn$name)
    if (file.exists(fn2))
      file.remove(fn2)
    file.rename(input$fn$datapath, fn2)
    fn <<- fn1
    isolate(mydata$spec <<- NULL)
    isolate(mydata$schron <<- NULL)
    isolate(mydata$group <<- NULL)
    sheets.nms <<- tryCatch(get_Sheets(fn2), error=function(e) return (e))
    if (inherits(sheets.nms, "error")) {
      updateSelectInput(session, "sheet", choices="")
      msg <- paste0("Cannot open Excel file.\nReason: ", sheets.nms$message)
      addNotification(msg, duration=30, type="error", closeButton=TRUE)
      fn1 <- input$fn$name
      sheets.nms <<- ""
      shinyjs::disable("downloadResults")
      return()
    } else {
      isolate(updateCheckboxInput(session, "exampleDF", value=FALSE))
      shinyjs::enable("downloadResults")
      shinyjs::enable("sheet")
      shinyjs::enable("autoSelect")
      shinyjs::enable("minCut")
      shinyjs::enable("maxCut")
      updateSelectInput(session, "sheet", choices=sheets.nms, selected=sheets.nms[1])
    }
  }, ignoreNULL = TRUE)
  
  observeEvent(input$sheet, {
    if (nchar(fn2) > 0 & nchar(input$sheet) > 0) {
      isolate(readIt(fn2, input$sheet, input, session))
      output$myPlot <- renderPlot({
        plotIt(mydata, input, session)
      })
    }
  }, ignoreNULL=TRUE)
  
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
      plotIt(mydata, input, session, flagIt)
      dev.off()
    } 
  )
  outputOptions(output, "varList", suspendWhenHidden=FALSE, priority=100)
  session$onSessionEnded(stopApp) # kill shiny on browser close
}
