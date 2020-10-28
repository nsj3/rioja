header <- dashboardHeader(title="riojaPlot", tags$li(a(href='riojaPlot.html', target='_blank', rel='noopener noreferrer', icon('question-circle-o'), title='Help', width=8, height=85), class='dropdown'))

D_ui <- dashboardPage(header, 
        dashboardSidebar(disable = TRUE),                       
        dashboardBody(
        tags$head(tags$style(HTML('
                           .skin-blue .main-header .logo { background-color: #3c8dbc; font-size: 2em; text-align: left; }
                              .skin-blue .main-header .logo:hover { background-color: #3c8dbc; }
                              .box {-webkit-box-shadow: none; -moz-box-shadow: none; box-shadow: none; }
                              .dropdown { padding-bottom: 0px: margin-bottom: 0px; font-size: 1.5em; border: 0; height: 50px;}
                              .main-header {max-height: 50px}
                              .main-header .logo {height: 50px;}
                              .sidebar-toggle {height: 50px; padding-top: 1px !important;}
                              .navbar {min-height:50px !important}
                              .error { color: red; padding: 0; margin-top: -20px;}
                              #symbSize, #barSize, #yMin, #yMax, #yInterval, #nameSize,
                              #nameAngle, #axisSize, #nZones, #yLabel, #exagMult, #minCut, #maxCut {height: 25px}
                              #autoSelect {margin-bottom: 10px} 
                              #shiny-notification-panel { position:fixed; top: calc(0%); left: calc(90%-150px); }
                              .shiny-notification {opacity: 1.0;}
                              .fa-test { background-image: question-circle.svg; height: 30px;}

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
             shinyjs::hidden(htmlOutput("exampleData")),
         width=NULL, solidHeader=TRUE),
         box(
            selectInput("sheet", "Select worksheet:", ""), 
            width=NULL, solidHeader=TRUE),
         box(downloadButton("downloadResults", "Save plot"), 
            radioButtons('saveType', "", choices=c('pdf', 'png', 'svg'), selected='png', inline=TRUE), 
         width=NULL, solidHeader=TRUE),
         box(p(tags$b(paste0("riojaPlot "), version)),
            paste0("Powered by "), 
              riojaURL, 
            p(paste0(" Version ", utils::packageDescription("rioja", fields="Version"), 
                  " (", utils::packageDescription("rioja", fields="Date"), ")")),
            "Comments & bug reports to ", SJURL, width=NULL, solidHeader=TRUE)
      ),
      column(width=10,
#        shinyjs::hidden(wellPanel(id="errorBox", htmlOutput("errorText", class="errorMsg"), width=NULL, solidHeader=TRUE)),
#        shinyjs::hidden(wellPanel(id="errorBox2", htmlOutput("errorText2", class="errorMsg"), width=NULL, solidHeader=TRUE)),
        box(plotOutput("myPlot", height="600px"), width=NULL, solidHeader=TRUE),

        fluidRow(
          tabBox(title='',id='tab1', width=10, height="250px",
            tabPanel('Variables', id='variables',      
              box( 
                 selectInput('yvar', 'Y axis', choices='', selected='', width=inputWidth, multiple=FALSE),
              solidHeader=TRUE, width=3),
              box( 
                 pickerInput('selTaxa', 'Select X variables', choices='', multiple=TRUE, 
                            options=pickerOptions(dropupAuto=FALSE, size=5, selectedTextFormat='count', 
                                                  countSelectedText="{0} of {1} variables selected")),
              solidHeader=TRUE, width=3),
              box(
                 checkboxInput("scalePC", "Scale for %", value=TRUE),
              solidHeader=TRUE, width=2),
              box(
                 actionButton('autoSelect', 'Auto Select X vars'),
                 fluidRow(
                   box(
                     splitLayout(
                       numericInput('minCut', 'Min cutoff', value=2, min=0, max=10, step=1, width=numericWidth2),
                       numericInput('maxCut', 'Max cutoff', value=150, min=50, max=200, step=10, width=numericWidth2)
                     ),
                   solidHeader=TRUE, width=12)
                 ), 
              solidHeader=TRUE, width=4),
            ),
            tabPanel(title='Settings', height="250px",
               box( 
                  checkboxGroupInput('style', 'Style', choices = c(Line=1, Symbols=2, Silhouette=3), selected=3),
               solidHeader=TRUE, width=2),
               box( 
                    radioButtons("hLine", "Show bar", choices=c("None", "Curve", "Full"), selected="Curve", inline=FALSE), 
                    checkboxInput('barTop', 'Bars on top', value=TRUE),
               solidHeader=TRUE, width=2),
               box(
                 numericInput('symbSize', 'Symbol size', value=0.5, min=0.2, max=4, step=0.1, width=numericWidth),
                 numericInput('barSize', 'Bar width', value=1, min=1, max=10, step=1, width=numericWidth), 
               solidHeader=TRUE, width=2),
               box( 
                  checkboxGroupInput('misc', 'Settings', c('Reverse Y axis'=1, 
                                                           'Show min/max'=4, 'Auto sort vars'=5), 
                                                   selected=c(1, 2)),
               solidHeader=TRUE, width=3),
               box( 
                  checkboxGroupInput('exag', 'Exaggeration', c('Show'=1, 'Auto Col'=2), selected=1), 
                  numericInput('exagMult', 'Exag mult.', value=2, min=1.2, max=10, step=0.2, width=numericWidth), 
               solidHeader=TRUE, width=2),
            ),
            tabPanel(title='Y axis', height="250px",
               box( 
                 numericInput('yMin', 'Y axis min', value=NA, width=numericWidth),
                 numericInput('yMax', 'Y axis max', value=NA, width=numericWidth),
              solidHeader=TRUE, width=2),
              box( 
                numericInput('yInterval', 'Y axis interval', value=NA, width=numericWidth),
                textInput('yLabel', 'Y axis Label'),
                solidHeader=TRUE, width=3)
            ),
          tabPanel(title='Colours', 
            box(
              spectrumInput("outlineCol", "Line", "black", width=colWidth, choices = colList, 
                            update_on='dragstop'),
              spectrumInput("fillCol", "Silhouette", "darkgreen", width=colWidth, choices = colList, 
                            update_on='dragstop'),
              solidHeader=TRUE, width=2),
            box(
              spectrumInput("hLineCol", "Bar", "#D9D9D9", width=colWidth, choices = colList, 
                            update_on='dragstop'), 
             spectrumInput("zoneCol", "Zone", "red", width=colWidth, choices = colList, 
                           update_on='dragstop'),
           solidHeader=TRUE, width=2),
           box(
             spectrumInput("symbCol", "Symbol", "#000000", width=colWidth, choices=colList, 
                           update_on='dragstop'),
             spectrumInput("exagCol", "Exaggeration", "#E0E0E0", width=colWidth, choices=colList, 
                           update_on='dragstop'),
             solidHeader=TRUE, width=2),
          ),
          tabPanel(title='Sizes',  
            box(
              numericInput('nameSize', 'X-var font size', value=1, min=.4, max=2, step=0.05, width=numericWidth), 
              numericInput('nameAngle', 'Rotate names', value=90, min=0, max=90, step=5, width=numericWidth), 
            solidHeader=TRUE, width=3),
            box(
              numericInput('axisSize', 'Axis font size', value=1, min=.4, max=2, step=0.05, width=numericWidth), 
              checkboxInput('italicise', 'Italicise names', value=FALSE),
              solidHeader=TRUE, width=3)
          ),
          tabPanel(title='Zonation', 
            box(
             checkboxInput('doClust', 'Add zonation', value=FALSE),
             radioButtons("showZones", "Show zones", choices=c("No", "Auto", "Choose"), selected="Auto", inline=FALSE), 
            solidHeader=TRUE, width=3),
            box(
             numericInput('nZones', 'Number of zones', value=2, min=2, max=10, step=1, width="120px"),
             radioButtons("dataTransformation", "Transform data?", choices=c("No", "Sqrt", "Scale"), selected="Sqrt", 
                          inline=FALSE), 
             solidHeader=TRUE, width=4),
          ),
          tabPanel(title='Groups', 
                   box(
                     spectrumInput("groupCol1", "Group 1", "darkgreen", width=colWidth, choices = colList, 
                                   update_on='dragstop'),
                     spectrumInput("groupCol2", "Group 2", "darkkhaki", width=colWidth, choices = colList, 
                                   update_on='dragstop'),
                     solidHeader=TRUE, width=2),
                   box(
                     spectrumInput("groupCol3", "Group 3", "darkorange", width=colWidth, choices = colList, 
                                   update_on='dragstop'),
                     spectrumInput("groupCol4", "Group 4", "darkred", width=colWidth, choices = colList, 
                                   update_on='dragstop'),
                     solidHeader=TRUE, width=2),
                   box(
                     spectrumInput("groupCol5", "Group 5", "deepskyblue", width=colWidth, choices = colList, 
                                   update_on='dragstop'),
                     spectrumInput("groupCol6", "Group 6", "firebrick", width=colWidth, choices = colList, 
                                   update_on='dragstop'),
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