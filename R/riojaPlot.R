utils::globalVariables(c("groupData", "cumulLine", "cumulLineCol",
                         "groupColours", "groupCex", "groupNames", ""))

riojaPlot <- function(x, y, selVars=NULL, groups=NULL, ...) {
   plotdata <- list()
   plotdata$spec <- x
   plotdata$chron <- y
   plotdata$selVars <- selVars
   plotdata$groups <- groups

   args <- list(...)
   argNames <- names(args)
   style <- makeStyle()
   validStyles <- names(style)
   for (i in argNames) {
      if (!(i %in% validStyles)) 
         stop(paste("Style ", i, "is not a valid riojaPlot style"))
      style[i] <- args[i]
   }
   splot1(plotdata, style)  
} 

listStyles <- function() {
  styles <- unlist(makeStyle())
  x <- data.frame(Style=names(styles), Value=styles)
  rownames(x) <- NULL
  invisible(x)
}

makeStyle <- function() {
   style <- list()
   style$yvarName <- ""
   style$secYvarName <- ""
   style$yLabel <- ""
   style$secyLabel <- ""
   style$showSecAxis <- FALSE
   style$scalePC <- TRUE
   style$scaleMinMax <- FALSE
   style$yrev <- TRUE
   style$yMin <- NA
   style$yMax <- NA
   style$yInterval <- NA
   style$secYMin <- NA
   style$secYMax <- NA
   style$secYInterval <- NA
   style$autoOrder <- "none"
   style$showBars <- TRUE
   style$showLines <- TRUE
   style$showPoly <- TRUE
   style$showSymbol <- FALSE
   style$showClust <- FALSE
   style$showZones <- "auto"
   style$ClustUseSelected <- FALSE
   style$lwdBar <- 1
   style$colBar <- "grey"
   style$barTop <- TRUE
   style$symbCol <- "black"
   style$fillCol <- "darkgreen"
   style$outlineCol <- "black"
   style$zoneCol <- "red"
   style$symbSize <- 0.5
   style$xAxisFontSize <- 0.7
   style$yAxisFontSize <- 0.8
   style$yLabelFontSize <- 1
   style$nameFontSize <- 1
   style$nameAngle <- 90
   style$tickLen <- -0.3
   style$cumulFontSize <- 0.7
   style$dataTransformation <- "sqrt"
   style$showExag <- FALSE
   style$exagCol <- "auto"
   style$exagMult <- 2
   style$nameStyleBreakLong <- TRUE
   style$nameStylenBreak <- 20
   style$nameStyleItalicise <- TRUE
   style$showGroups <- FALSE
   style$showCumul <- FALSE
   style$groupCol1 <- "darkgreen"
   style$groupCol2 <- "darkkhaki"
   style$groupCol3 <- "darkorange"
   style$groupCol4 <- "darkred"
   style$groupCol5 <- "deepskyblue"
   style$groupCol6 <- "darkgrey"
   style
}

splot1 <- function(mydata, style) 
{
   if (is.null(mydata$spec) | is.null(mydata$chron) )
      return();
  
   if (!is.null(mydata$selVars) & length(mydata$selVars) > 2) {
      if (!(any(mydata$selVars %in% colnames(mydata$spec))))
         stop("Some variables listed in selVars are not found in the data.")
      d <- mydata$spec[, mydata$selVars]
   } else {
     d <- mydata$spec
   }
  
   yvarName <- style$yvarName
   if (nchar(stringr::str_trim(yvarName)) > 1) {
     if (!(yvarName %in% colnames(mydata$chron))) {
        stop(paste(yvarName, " is not a column in the chronology data."))
     }
     yvar <- mydata$chron[, yvarName, drop=FALSE]
   } else {
     yvar <- mydata$chron[, 1, drop=FALSE]
     yvarName <- colnames(mydata$chron)[1]
   }


   style$poly.line <- NA 
   if (style$showLine)
      poly.line <- style$outlineCol
     
   clust <- NULL
   if (style$showClust) {
     if (style$ClustUseSelected) {
        d2 <- d
     } else {
        d2 <- mydata$spec
     }
     if (style$dataTransformation == "sqrt") {
        d2 <- sqrt(d2)
     } 
     if (style$dataTransformation == "scale") {
        d2 <- scale(d2, TRUE, TRUE)
     }
     diss <- dist(d2)
     clust <- chclust(diss)
     if (style$showZones == "auto") {
        x <- bstick(clust, plot=FALSE)
        x2 <- x$dispersion <= x$bstick
        style$showZones <- which(x2)[1]
        if (style$showZones < 2) {
          print("There are no significant zones in these data.")
        }
     } 
   } 

   if (style$exagMult < 1) {
      style$exagMult <- 1.0
      style$showExag <- FALSE
   }
  
   style$xNames <- colnames(d)
   if (style$nameStyleBreakLong) {
     style$xNames <- sjmisc::word_wrap(style$xNames, style$nameStylenBreak)
   }
   if (style$nameStyleItalicise) {
     style$xNames <- as.expression(sapply(style$xNames, function(x) bquote(italic(.(x))) ))
   }
   
   style$groupColours <- rep(style$fillCol, ncol(d))
   groupID <- rep(1, ncol(d))

   if ((style$showGroups | style$showCumul) & !is.null(mydata$groups)) {
      tmp <- data.frame(Names=colnames(d))
      colnames(mydata$groups)[1] <- "Names"
      tmp <- left_join(tmp, mydata$groups, by="Names")
      if (any(is.na(tmp[, 2, drop=TRUE])))
        stop("Some variable names not found in the grouping.")
      if (!is.factor(tmp[, 2])) {
            stop("Grouping variable must be a factor.")
      }
      groupID <- as.integer(tmp[, 2, drop=TRUE])
      groupColours <<- c(style$groupCol1, style$groupCol2, style$groupCol3, 
                        style$groupCol4, style$groupCol5)
      groupColours <- c(style$groupCol1, style$groupCol2, style$groupCol3, 
                        style$groupCol4, style$groupCol5)
      groupNames <<- levels(tmp[, 2])
      groupNames <- levels(tmp[, 2])
      if (style$showGroups)
         style$groupColours <- groupColours[groupID]
   }

   if (!is.null(dim(yvar)))
     yvar <- yvar[, 1, drop=FALSE]
   style$yLabels <- NULL
   style$ytks <- NULL
   style$ytks1 <- NULL
   style$ytks2 <- NULL
   ylim <- NULL
   if (is.character(yvar[, 1, drop=TRUE])) {
      style$yLabels <- yvar[, 1, drop=TRUE]
      yvar <- data.frame(SampleNo=1:length(yvar[, 1, drop=TRUE]))
      style$ytks1 <- yvar[, 1, drop=TRUE]
   } else {
      ylim <- range(yvar[, 1], na.rm=TRUE)
      if (!is.na(style$yMin) & is.na(style$yMax)) {
         ylim[1] <- style$yMin
         ylim[2] <- max(yvar[, 1], na.rm=TRUE)
      } else if (!is.na(style$yMax) & is.na(style$yMin)) {
         ylim[2] <- style$yMax
         ylim[1] <- min(yvar[, 1], na.rm=TRUE)
      } else if (!is.na(style$yMax) & !is.na(style$yMin)) {
         ylim[1] <- style$yMin
         ylim[2] <- style$yMax
      }
      if (!is.null(ylim) & !is.na(style$yInterval)) {
         style$ytks1 <- seq(ylim[1], ylim[2], by=style$yInterval)
      } 
   }

   secYvarName <- style$secYvarName
   doSecYvar <- FALSE
   secYvar <- NULL
   yLab <- yvarName
   if (nchar(style$yLabel)>0) {
      yLab <- style$yLabel
   }
   if (style$showSecAxis & nchar(stringr::str_trim(secYvarName)) > 0) {
      if (yvarName != secYvarName) {
         if (!(secYvarName %in% colnames(mydata$chron))) {
           stop(paste(secYvarName, " is not a column in the chronology data."))
         }
         secYvar <- mydata$chron[, secYvarName, drop=FALSE]
         if (!is.numeric(secYvar[, 1, drop=TRUE])) {
            print("Secondary Y axis variable must be numeric, not character.")
         } else {
            doSecYvar <- TRUE
            if (nchar(style$secyLabel)>0) {
               secYvarName <- style$secyLabel
            }
            yvar <- as.data.frame(cbind(yvar, secYvar))
            yLab <- c(yLab, secYvarName)
         }
      }
   }
   ylim2 <- NULL
   style$yks <- style$yks1
   if (doSecYvar) {
      ylim2 <- range(yvar[, 2], na.rm=TRUE)
      if (!is.na(style$secYMin) & is.na(style$secYMax)) {
         ylim2[1] <- style$secYMin
         ylim2[2] <- max(yvar[, 2], na.rm=TRUE)
      } else if (!is.na(style$secYMax) & is.na(style$secYMin)) {
         ylim2[2] <- style$secYMax
         ylim2[1] <- min(yvar[, 2], na.rm=TRUE)
      } else if (!is.na(style$secYMax) & !is.na(style$secYMin)) {
         ylim2[1] <- style$secYMin
         ylim2[2] <- style$secYMax
      }
      if (!is.null(ylim2)) {
         if (is.na(style$secYInterval)) {
            style$ytks2 <- pretty(ylim2, n=10)
         } else {
            style$ytks2 <- seq(ylim2[1], ylim2[2], by=style$secYInterval)
         }
         style$ytks <- list(style$`ytks1`, style$`ytks2`)
      } 
   }
#   size <<- dev.size()
   nms <- colnames(d)

# Groups   
   funlist <- lapply(1:ncol(d), function(x) NULL)
   
   if (style$showCumul) {
      groupData <<- t(apply(d, 1, 
                            function(x) cumsum(tapply(unlist(x), 
                            groupID, sum, na.rm=TRUE))))
      groupData <- t(apply(d, 1, 
                            function(x) cumsum(tapply(unlist(x), 
                            groupID, sum, na.rm=TRUE))))
      tt <- table(groupID)
      if (length(tt) == 1) {
         groupData <<- t(groupData)
         groupData <- t(groupData)
         colnames(groupData) <<- names(tt)
         colnames(groupData) <- names(tt)
      }
      d <- data.frame(d, Cumulative=c(100, rep(0, nrow(d)-1)))
      style$xNames <- c(style$xNames, "Cumulative")
      nCol <- ncol(d)
      funlist <- lapply(1:(nCol), function(x) NULL)
      funlist[[nCol]] <- plotCumul
      style$groupColours <- c(style$groupColours, NA)
      cumulLine <<- style$showLine
      cumulLineCol <<- style$outlineCol
      cumulLine <- style$showLine
      cumulLineCol <- style$outlineCol
   }
   
   fin <- par("fin")
   xSpace <- 0.1 / fin[1]
   groupCex <<- style$cumulFontSize
   groupCex <- style$cumulFontSize
   
   style$lineColour <- style$outlineCol
   if (style$showGroups & !style$showPoly & style$showLine) {
      style$lineColour <- style$groupColours
   }

   oldfig <- par("fig")
   oldmar <- par("mar")
   oldusr <- par("usr")

   on.exit({ par(mar=oldmar); par(fig=oldfig); par(usr=oldusr) }) 
   
#   if (style$nameStyleBreakLong) {
#      yLab <- sjmisc::word_wrap(yLab, style$nameStylenBreak)
#   }

   x <- splot2(d, yvar = yvar, y.rev=style$yrev, scale.percent=style$scalePC, 
                plot.bar=style$showBars, plot.line=style$showLine, plot.poly=style$showPoly, plot.symb=style$showSymbol, 
                col.poly=style$groupColours, col.bar=style$colBar, lwd.bar=style$lwdBar, 
                col.symb=style$symbCol, col.poly.line=style$poly.line, col.line=style$lineColour, 
                symb.cex=style$symbSize, exag=style$showExag, wa.order=style$autoOrder, bar.back=!style$barTop, 
                clust=clust, cex.xlabel=style$nameFontSize, srt.xlabel=style$nameAngle, 
                ylabel=yLab, cex.yaxis=style$yAxisFontSize, cex.axis=style$xAxisFontSize, 
                cex.ylabel=style$yLabelFontSize, scale.minmax=style$scaleMinMax, ylim=ylim, y.tks=style$ytks, 
                y.tks.labels=style$yLabels, col.bg=NULL, col.exag=style$exagCol, exag.mult=style$exagMult, 
                exag.alpha=0.15, x.names=style$xNames, fun2=funlist, xSpace=xSpace, tcl=style$tickLen)

      
#   retVal <- tryCatch(x <- splot2(d, yvar = yvar, y.rev=style$yrev, scale.percent=style$scalePC, 
#                plot.bar=style$bar, plot.line=style$line, plot.poly=style$poly, plot.symb=style$symbol, 
#                col.poly=style$groupColours, col.bar=style$colBar, lwd.bar=style$lwdBar, 
#                col.symb=style$symbCol, col.poly.line=style$poly.line, col.line=style$lineColour, 
#                symb.cex=style$symbSize, exag=style$exag, wa.order=style$autoOrder, bar.back=!style$barTop, 
#                clust=clust, cex.xlabel=style$nameSize, srt.xlabel=style$nameAngle, 
#                ylabel=yLab, cex.yaxis=style$yAxisSize, cex.axis=style$xAxisSize, 
#                cex.ylabel=style$yLabelSize, scale.minmax=style$scaleMinMax, ylim=ylim, y.tks=style$ytks, 
#                y.tks.labels=style$yLabels, col.bg=NULL, col.exag=style$exagCol, exag.mult=style$exagMult, 
#                exag.alpha=0.15, x.names=style$xNames, fun2=funlist, xSpace=xSpace),
#      error=function(e) return(e))

   
#   if (inherits(retVal, "error")) {
#      print("Error")
#      return()
#   }
   
   if (style$showClust & style$showZones > 1) {
      addClustZone2(x, clust, style$showZones, col=style$zoneCol, yaxs="i")
   }
   invisible(x)
}


splot2 <- function(d, yvar = NULL, scale.percent = FALSE, graph.widths=1, minmax=NULL, 
                  scale.minmax=TRUE, xLeft=NULL, xRight=NULL, yBottom=NULL, yTop=NULL, 
                  title="", cex.title=1.8, y.axis=TRUE, x.axis=TRUE, min.width=5, 
                  ylim=NULL, y.rev=FALSE, y.tks=NULL, y.tks.labels=NULL, ylabel=NULL,
                  cex.ylabel=1, cex.yaxis=0.8, xSpace=0.01, x.pc.inc=10, x.pc.lab=TRUE, 
                  x.pc.omit0=TRUE, wa.order="none", plot.line=TRUE, col.line="black", 
                  lwd.line=1, col.symb="black", plot.bar=TRUE, lwd.bar=1, col.bar="grey",
                  sep.bar=FALSE, bar.back=FALSE, plot.poly=FALSE, col.poly="grey", 
                  col.poly.line=NA, lwd.poly=1, plot.symb=FALSE, symb.pch=19, symb.cex=1,
                  x.names=NULL, cex.xlabel=1.0, srt.xlabel=90, mgp=c(3, cex.axis/3, 0.2),
                  ylabPos=NULL, cex.axis=0.8, clust=NULL, clust.width=0.1, orig.fig=NULL, 
                  exag=FALSE, exag.mult=5, col.exag="grey90", exag.alpha=0.2, 
                  col.bg=NULL, fun1=NULL, fun2=NULL, add=FALSE, omitMissing=TRUE, ...)
{
   errorMsg <- function(msg) {
     if (shiny_running()) {
       stop(msg)
     } else {
       stop(msg)
     }
   }
   d <- as.data.frame(d)
   fcall <- match.call(expand.dots=TRUE)
   if (!is.null(clust)) {
     if (!is(clust, "chclust"))
        errorMsg("clust must be a chclust object")
   }
   if (!is.null(clust)) {
      if (is.null(xRight))
         xRight <- 1
      xRight = xRight - clust.width
   }
   doSecYvar <- FALSE

   if (is.null(yvar)) {
      yvar <- data.frame(SampleNo=1:nrow(d))
      if (is.null(ylim)) {
         ylim=c(1, nrow(d))
      }
   } else {
      if (is.null(dim(yvar))) {
         nm <- substitute(yvar)
         yvar <- data.frame(tmp=yvar)
         colnames(yvar) <- as.character(nm)
      }
      else {
         yvar <- as.data.frame(yvar)
         if (ncol(yvar)>1)
            doSecYvar <- TRUE
      }
   }
  
   yNames <- c("", "")
   if (!is.null(ylabel)) {
      if (length(ylabel)==1)
         yNames <- c(ylabel, "")
      else
         yNames <- ylabel[1:2]
   } else {
      yNames <- colnames(yvar)     
   }
  
   if (is.null(x.names))
      x.names=colnames(d)   
   if (is.null(ylim)) {
      ylim = range(yvar[, 1], na.rm=TRUE)
   } else {
      if (is.na(yvar[1, 1]))
         ylim[1] <- min(yvar[, 1], na.rm=TRUE)
      if (is.na(ylim[2]))
         ylim[2] <- max(yvar[, 1], na.rm=TRUE)
   }
   oldfig = par("fig")
   oldmai <- par("mai")
   if (is.null(orig.fig)) {
      orig.fig = par("fig")
   }
   if (exag.mult < 1.0)
      exag <- FALSE
   nsp <- ncol(d)
   nsam <- nrow(d)
   if (scale.percent==TRUE & length(x.pc.inc) > 1) {
      if (length(x.pc.inc) != nsp) 
         errorMsg("length of x.pc.inc should equal number of curves")
   } else {
      x.pc.inc <- rep(x.pc.inc[1], nsp)
   }
   if (!is.null(minmax)) {
     if (ncol(minmax) != 2) 
        errorMsg("minmax should have 2 columns")
     if (nrow(minmax) != nsp) 
        errorMsg("number of rows of minmax should equal number of curves")
   }
   par(mai = c(0, 0, 0, 0))
   if (length(graph.widths) == 1)
      graph.widths <- rep(1, nsp)
   if (length(graph.widths) != nsp) 
      errorMsg("Length of graph.widths should equal number of curves")
   if (length(exag) == 1)
      exag <- rep(exag[1], nsp)
   if (length(exag) != nsp) 
      errorMsg("Length of exag should equal number of curves")
   if (length(exag.mult) == 1)
      exag.mult <- rep(exag.mult[1], nsp)
   if (length(exag.mult) != nsp) 
      errorMsg("Length of exag.mult should equal number of curves")
   if (length(col.exag) == 1)
      col.exag <- rep(col.exag[1], nsp)
   if (length(col.exag) != nsp) 
      errorMsg("Length of col.exag should equal number of curves")
   if (!is.null(fun1)) {
      if (length(fun1) == 1)
         fun1 <- lapply(1:nsp, function(x) fun1)
      if (length(fun1) != nsp)
         errorMsg("Length of fun1 should equal number of curves")
   }
   if (!is.null(fun2)) {
      if (length(fun2) == 1)
         fun2 <- lapply(1:nsp, function(x) fun2)
      if (length(fun2) != nsp)
         errorMsg("Length of fun2 should equal number of curves")
   }
   if (length(x.axis) == 1)
      x.axis <- rep(x.axis[1], nsp)
   if (length(x.axis) != nsp)
      errorMsg("Length of x.axis should equal number of curves")
   cc.line <- rep(col.line, length.out=nsp)
   if (sep.bar)
      cc.bar <- rep(col.bar, length.out=nsam)
   else      
      cc.bar <- rep(col.bar, length.out=nsp)
   cc.poly <- rep(col.poly, length.out=nsp)
   cc.poly.line <- rep(col.poly.line, length.out=nsp)
#  if(plot.poly)
#    plot.line <- FALSE
   make.col <- function(x, alpha) {
      apply(col2rgb(x)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha))
   }
   if (col.exag[1] == "auto")
      col.exag <- make.col(cc.poly, exag.alpha)
   inc <- 0.002
   if (wa.order == "topleft" || wa.order == "bottomleft") {
      V1 <- 1:nrow(d)
      colsum <- base::colSums(d, na.rm=TRUE)
#    opt <- (t(d) %*% yvar)/colsum
     opt <- (t(d) %*% V1)/colsum
     if ((wa.order == "topleft" & !y.rev) | (wa.order == "bottomleft" & y.rev))
        opt.order <- rev(order(opt))
     else opt.order <- order(opt)
     d <- d[, opt.order]
     if (!is.null(minmax)) 
        minmax <- minmax[opt.order, ]
     if (!is.null(x.names))
        x.names <- x.names[opt.order]
     graph.widths <- graph.widths[opt.order]
     exag <- exag[opt.order]
     exag.mult <- exag.mult[opt.order]
     if (!is.null(fun1))
        fun1 <- fun1[opt.order]
     if (!is.null(fun2))
        fun2 <- fun2[opt.order]
     x.axis <- x.axis[opt.order]
     cc.poly <- cc.poly[opt.order]
     cc.poly.line <- cc.poly.line[opt.order]
     cc.line <- cc.line[opt.order]
     if (!sep.bar)
        cc.bar <- cc.bar[opt.order]
   }
   if (scale.percent) {
      colM <- apply(d, 2, max, na.rm=TRUE)
      colM <- floor((colM + 4.9)/5) * 5
      colM[colM < min.width] <- min.width
      colM.sum <- sum(colM, na.rm=TRUE)
   } else {
      colM.sum <- sum(graph.widths, na.rm=TRUE)
      colM <- graph.widths
   }
  
# determine fig margins  
  
   maxlen <- max(sapply(x.names, function(x) strwidth(x, units="figure", cex=cex.xlabel))) 
   maxlen <- max(maxlen, strwidth(yNames[1], units="figure", cex=cex.xlabel))
   if (doSecYvar)
      maxlen <- max(maxlen, strwidth(yNames[2], units="figure", cex=cex.xlabel))
   fin <- par("fin")
   plotRatio <- fin[1] / fin[2]
   xLabSpace <- 0.1
   if (is.null(yTop)) {
      xlSpace <- xLabSpace / fin[2]
      yTop <- 1.0 - (maxlen * plotRatio * cos(pi/180 * (90-srt.xlabel))) - xlSpace - 0.01
      yTop <- min(yTop, 0.95)
   }
   
   if (is.null(xLeft)) {
      if (!is.null(y.tks.labels))
        ylabs <- y.tks.labels
      else 
#        ylabs <- as.character(yvar[, 1])
        ylabs <- pretty(yvar[, 1], n=10)
      incX <- strheight("M", units="figure", cex=1) / plotRatio # distance to axis values
      mx1 <- max(sapply(ylabs, function(x) strwidth(x, units='figure', cex=cex.yaxis))) # width of axis labels
      xLeft <- incX + mx1 + 0.02 / plotRatio

# without label
      incX <- strwidth("0", units='figure', cex=1)
      if (doSecYvar)
         xLeft <- mx1 + incX * 4
      else 
         xLeft <- mx1 + incX * 3
      
# now label
#      if (nchar(stringr::str_trim(yNames[1])) > 0 & !doSecYvar) {
      if (nchar(yNames[1]) > 0 & !doSecYvar) {
         line2fig <- strheight(yNames[1], units='figure', cex=1) / plotRatio
         if (is.null(ylabPos)) {
            ylabPos <- 1 + mx1 / line2fig
         }
         xLeft <- xLeft + (line2fig + line2fig * cex.ylabel) 
      }
   }
   ylab2 <- NULL
   yAxis2Pos <- 0
   tcll <- -.3
   spc <- 0
   
#   if ("tcl" %in% names(fcall))
#      tcll <- eval(fcall$tcl)
#   spc <- 0
#   if ("las" %in% names(fcall)) {
#      if ((eval(fcall$las)) == 2)
#        spc = 0.3
#  }
   args <- list(...)
   if ("tcl" %in% names(args)) {
       tcll <- args[["tcl"]]
   }
   if ("las" %in% names(args)) {
       if(args[["las"]]==2) {
         spc <- 0.3
       }
   }

   if (y.axis & doSecYvar) {
      if (is(y.tks, "list") & !is.null(y.tks[[2]])) {
         y.tks2 <- y.tks[[2]]
         xout <- y.tks2
      } else {
         xout <- pretty(yvar[, 2], 10)
      }
      if (as.integer(R.Version()$major) > 3)
         ylab2 <- stats::approx(yvar[, 2, drop=TRUE], yvar[, 1, drop=TRUE], xout=xout, na.rm=TRUE)
      else 
         ylab2 <- stats::approx(yvar[, 2, drop=TRUE], yvar[, 1, drop=TRUE], xout=xout)
      mx1 <- max(sapply(as.character(ylab2$x), function(x) strwidth(x, units='figure', cex=cex.yaxis)))
      yAxis2Pos <- mx1 + incX * 3 
      incX <- strwidth("0", units='figure', cex=1)
      xLeft <- xLeft + yAxis2Pos
   }

   if (is.null(clust) & is.null(xRight)) {
        xRight <- 1.0
        xLen <- xRight - xLeft
        xInc <- xLen - ((nsp + 1) * xSpace)
        n <- length(colM)
        inc <- xInc * colM[n]/colM.sum
        wid <- strwidth(x.names[length(x.names)], units='figure', 
                               cex=cex.xlabel) * 0.9 * sin(pi/180 * (90-srt.xlabel))
        if (wid > inc) {
          xRight <- 1 - (wid-inc)
        }
   } 
   if (is.null(yBottom)) {
      yBottom <- 0.06
   }

   xLen <- xRight - xLeft
   xInc <- xLen - ((nsp + 1) * xSpace)
   inc <- xInc/colM.sum
   if (inc < 0.0)
     errorMsg("Too many variables, curves will be too small.")
   x1 <- xLeft
    #    par(fig = c(x1, x1+0.4, yStart, yTop))
   if (y.rev) {
     tmp <- ylim[1]
     ylim[1] <- ylim[2]
     ylim[2] <- tmp
   }
   if (y.axis) {

     if (doSecYvar) {
       par(fig = figCnvt(orig.fig, c(yAxis2Pos, yAxis2Pos+0.2, yBottom, yTop)), new=add)
       plot(0, cex = 0.5, xlim = c(0, 1), axes = FALSE, type = "n", xaxs="i", yaxs = "i", ylim = ylim, tcl=tcll, ...)
       axis(side = 2, las = 1, at = ylab2$y, labels = as.character(ylab2$x), cex.axis=cex.yaxis, xpd=FALSE, 
            tcl=tcll, mgp=c(3, 0.6, 0))
       addName(yNames[2], xLabSpace, srt.xlabel, cex.xlabel, y.rev, offset=-2)     
       add <- TRUE
     }

     par(fig = figCnvt(orig.fig, c(x1, x1+0.2, yBottom, yTop)), new=add)
     plot(NA, cex = 0.5, xlim = c(0, 1), axes = FALSE, type = "n", xaxs="i", yaxs = "i", ylim = ylim, tcl=tcll, ...)
     usr1 <- par("usr")
     if (is.null(y.tks))
       y.tks <- axTicks(2)
     else if (mode(y.tks)=="list") {
       y.tks <- y.tks[[1]]
     }
     if (is.null(y.tks.labels))
       y.tks.labels <- as.character(y.tks)
     else
       y.tks.labels <- y.tks.labels
     ax <- axis(side = 2, las = 1, at = y.tks, labels = as.character(y.tks.labels), cex.axis=cex.yaxis, xpd=NA, 
                tcl=tcll, mgp=c(3, 0.6, 0))
     x1 <- x1 + xSpace
#     mtext(title, adj = 0, line = 5, cex = cex.title)
#     if (nchar(stringr::str_trim(yNames[1])) > 0) {
     if (nchar(yNames[1]) > 0) {
        if (!doSecYvar) {
           mtext(yNames[1], side=2, line=ylabPos, cex=cex.ylabel)
        } else {
           addName(yNames[1], xLabSpace, srt.xlabel, cex.xlabel, y.rev, offset=-2)     
        }
     }
   }
   ty <- ifelse(plot.line, "l", "n")
   
 
   figs <- vector("list", length=nsp)
   usrs <- vector("list", length=nsp)
  
   for (i in 1:nsp) {
# omit missing values  
     y_var <- yvar[, 1, drop=TRUE]
     x_var <- d[, i, drop=TRUE]
  
     cumulPlot <- FALSE
     if (toupper(x.names[i])=="CUMULATIVE") {
        cumulPlot <- TRUE
     }
     
     nsam2 <- nsam
     if (omitMissing) {
        miss <- is.na(y_var) | is.na(x_var)
        nsam2 <- sum(!miss)
        if (nsam2 < nsam) {
           y_var <- y_var[!miss]
           x_var <- x_var[!miss]
           if (sep.bar) {
              cc.bar <- cc.bar[!miss]
           }
        }
     }
     par(new = TRUE)
     par(lend = "butt")
     if (scale.percent) {
        inc2 <- inc * colM[i]
        par(fig = figCnvt(orig.fig, c(x1, x1 + inc2, yBottom, yTop)))
        plot(0, 0, cex = 0.5, xlim = c(0, colM[i]), axes = FALSE, 
           xaxs = "i", type = "n", yaxs = "i", ylim = ylim, xlab="", ylab="", ...)
        if (!is.null(col.bg))
           rect(par("usr")[1],ylim[1],par("usr")[2],ylim[2], col=col.bg, border=NA)
        if (!is.null(fun1[[i]])) {
           fun1[[i]](x=x_var, y=y_var, i=i, nm=x.names[i])
        }
        if (plot.poly & exag[i] & !cumulPlot) {
           y <- c(y_var[1], y_var, y_var[nsam2])
           x2 <- c(0, x_var*exag.mult[i], 0)
           polygon(x2, y, col = col.exag[i], border = NA, xpd=FALSE)
        }        
        if (bar.back  & !cumulPlot) {
           if (is.logical(plot.bar)) {
              if (plot.bar) {
                if (sep.bar) {
                   segments(rep(0, nsam2), y_var, x_var, y_var, lwd = lwd.bar, col = cc.bar)
                } else {
                   segments(rep(0, nsam2), y_var, x_var, y_var, lwd = lwd.bar, col = cc.bar[i])
                }
              }
           } else {
              if (plot.bar=="full") {
                 abline(h=y_var, col=cc.bar, lwd=lwd.bar)
              }
           }
        }
        if (plot.poly) {
           y <- c(y_var[1], y_var, y_var[nsam2])
           x <- c(0, x_var, 0)
           polygon(x, y, col = cc.poly[i], border = NA, lwd=lwd.poly)
        }
        if (!bar.back & !cumulPlot) {
           if (is.logical(plot.bar)) {
              if (plot.bar) {
                 if (sep.bar) {
                    segments(rep(0, nsam2), y_var, x_var, y_var, lwd = lwd.bar, col = cc.bar)
                 } else {
                    segments(rep(0, nsam2), y_var, x_var, y_var, lwd = lwd.bar, col = cc.bar[i])
                 }
              }
           } else {
              if (plot.bar=="full") {
                 abline(h=y_var, col=cc.bar, lwd=lwd.bar)
              }
           }
       }
       lines(c(0, 0), c(min(y_var, na.rm=TRUE), max(y_var, na.rm=TRUE)), ...)
       if (ty == "l") 
          lines(x_var, y_var, col = cc.line[i], lwd = lwd.line)
       if (plot.symb & !cumulPlot) {
          points(x_var, y_var, pch=symb.pch, cex=symb.cex, col=col.symb, xpd=FALSE)
       }
       if (!is.null(fun2[[i]])) {
          fun2[[i]](x=x_var, y=y_var, i=i, nm=x.names[i])
       }
       xlabb <- seq(0, colM[i], by = x.pc.inc[i])
       if (x.axis[i]) {
          if (x.pc.lab) {
             xlabbt <- as.character(xlabb)
             if (x.pc.omit0)
                xlabbt[1] <- ""
             mgpX <- if (is.null(mgp)) { c(3,max(0.0, spc-tcll), 0.3 ) } else { mgp }
             axis(side = 1, at = xlabb, labels = xlabbt, mgp=mgpX, cex.axis=cex.axis, tcl=tcll, ...)
         } else {
             axis(side = 1, at = xlabb, labels = FALSE, mgp=mgpX, ...)
         }
       }
       x1 <- x1 + inc2 + xSpace
     } else {
       inc2 <- inc * colM[i]
       par(fig = figCnvt(orig.fig, c(x1, min(1, x1 + inc2, na.rm=TRUE), yBottom, yTop)))
       if (!is.null(minmax)) {
          plot(x_var, y_var, cex = 0.5, axes = FALSE, xaxs = "i", 
               type = "n", yaxs = "i", ylim = ylim, xlim=c(minmax[i, 1], minmax[i,2]), tcl=tcll, ...)
       } else {
          plot(x_var, y_var, cex = 0.5, axes = FALSE, xaxs = "i", 
             type = "n", yaxs = "i", ylim = ylim, tcl=tcll, ...)
       }
       if (!is.null(col.bg))
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col=col.bg)
       tks <- axTicks(1)
       us <- par("usr")
       if (!is.null(fun1[[i]])) {
          fun1[[i]](x=x_var, y=y_var, i=i, nm=x.names[i])
       }
       if (plot.poly & exag[i] & !cumulPlot) {
          y <- c(y_var[1], y_var, y_var[nsam2])
          x2 <- c(us[1], x_var*exag.mult[i], us[1])
          polygon(x2, y, col = col.exag[i], border = NA)
       }
       if (bar.back & !cumulPlot) {
          if (is.logical(plot.bar)) {
            if (plot.bar) {
              if (sep.bar) {
                 segments(rep(us[1], nsam2), y_var, x_var, y_var, lwd = lwd.bar, col = cc.bar)
              } else {
                 segments(rep(us[1], nsam2), y_var, x_var, y_var, lwd = lwd.bar, col = cc.bar[i])                            }
            }
          } else {
             if (plot.bar=="full") {
                abline(h=y_var, col=cc.bar, lwd=lwd.bar)
             }
          }
       }
       if (plot.poly) {
          y <- c(y_var[1], y_var, y_var[nsam2])
          x <- c(us[1], x_var, us[1])
          if (exag[i]) {
             x2 <- c(us[1], x_var*exag.mult[i], us[1])
             polygon(x2, y, col = col.exag[i], border = NA)
          }
          polygon(x, y, col = cc.poly[i], border = NA, lwd=lwd.poly)
       }
       if (!bar.back & !cumulPlot) {
          if (is.logical(plot.bar)) {
             if (plot.bar) {
               if (sep.bar) {
                  segments(rep(us[1], nsam2), y_var, x_var, y_var, lwd = lwd.bar, col = cc.bar)
               } else {
                  segments(rep(us[1], nsam2), y_var, x_var, y_var, lwd = lwd.bar, col = cc.bar[i])                              }
             }
          } else {
             if (plot.bar=="full") {
                abline(h=y_var, col=cc.bar, lwd=lwd.bar)
             }
          }
       }
       lines(c(us[1], us[1]), c(min(y_var, na.rm=TRUE), max(y_var, na.rm=TRUE)), ...)
       if (ty == "l") 
          lines(x_var, y_var, col = cc.line[i], lwd = lwd.line)
       if (plot.symb & !cumulPlot) {
          points(x_var, y_var, pch=symb.pch, cex=symb.cex, col=col.symb, xpd=FALSE)
       }
       if (!is.null(fun2[[i]])) {
          fun2[[i]](x=x_var, y=y_var, i=i, nm=x.names[i])
       }
       mgpX <- if (is.null(mgp)) { c(3, max(0.0, spc-tcll), 0.3) } else { mgp }
       if (x.axis[i]) {
          if (scale.minmax) {
             nn <- length(axTicks(1))
             tk <- c(axTicks(1)[1], axTicks(1)[nn])
             axis(side = 1, at = tk, labels = as.character(tk), cex.axis=cex.axis, mgp=mgpX, tcl=tcll, ...)
          } else {
             axis(side = 1, cex.axis=cex.axis, mgp=mgpX, ...)
          }
       }
       x1 <- x1 + inc2 + xSpace
     }
     usr2 <- par("usr")
     tks1 <- usr2[1]
     fin <- par("fin")
     rD <- abs((usr2[4] - usr2[3]))
     rF <- fin[2]
     r <- rD/rF*xLabSpace
     pos <- usr2[4] + r
     if (y.rev)
        pos <- usr2[4]-r
     if (!cumulPlot) {
        par("lheight" = 0.7)
        if (srt.xlabel < 90)
           text(tks1[1], pos, labels=x.names[i], adj=c(0, 0), srt=srt.xlabel, cex = cex.xlabel, xpd=NA)
        else
           text(tks1[1], pos, labels=x.names[i], adj=c(0, 1), srt=srt.xlabel, cex = cex.xlabel, xpd=NA)
        par("lheight" = 1)
     }
     usrs[[i]] <- usr2   
     figs[[i]] <- par("fig")
   }
   if (!is.null(clust)) {
      par(fig = figCnvt(orig.fig, c(x1, xRight+clust.width, yBottom, yTop)))
      par(mar=c(0,0,0,0))
      par(new = TRUE)
      if (y.rev)
         xl <- rev(ylim)
      else
         xl <- ylim
     plot(clust, xvar=yvar[, 1, drop=TRUE], horiz=TRUE, x.rev=y.rev, labels=rep("", length(yvar[, 1, drop=TRUE])), 
          hang=-1, mgp=mgpX, cex.axis=cex.axis, xlim=xl, yaxs="i", xpd=FALSE, tcl=tcll, ...)
   }
   par(mai = oldmai)
   oldfig[oldfig < 0] <- 0
   par(fig = oldfig)
   ll <- list(call=fcall, box=c(xLeft=xLeft, xRight=xRight, yBottom=yBottom, yTop=yTop), usr = usr1, 
              yvar=yvar[, 1, drop=TRUE], ylim=ylim, y.rev=y.rev, figs=figs, usrs=usrs)
   invisible(ll)
}

addClustZone2 <- function(x, clust, nZone, ...) {
  oldpar <- par(c("fig", "mar", "usr"))
  par(fig=x$box)
  par(mar=c(0,0,0,0))
  par(usr=c(0, 1, x$usr[3], x$usr[4]))
  cc <- cutree(clust, k=nZone)
  zn <- which(diff(cc)>0)
  #   if (x$yaxt.rev)
  #      x$yvar <- rev(x$yvar)
  zone <- (x$yvar[zn] + x$yvar[zn+1]) / 2
  r <- range(c(x$usr[3], x$usr[4]))
  sel <- which (zone >= r[1] & zone <= r[2])
  if (length(sel) > 0) {
    zone <- zone[sel]
    segments(0, zone, 1, zone, xpd=NA, ...)
  }
  par(oldpar)
}

addName <- function(x, xLabSpace, srt.xlabel, cex.xlabel, y.rev, offset=0)
{
    usr2 <- par("usr")
    tks1 <- usr2[1]
    fig <- par("fin")
    rD <- abs((usr2[4] - usr2[3]))
    rF <- fig[2]
    r <- rD/rF * xLabSpace
    yPos <- usr2[4] + r 
    if (y.rev)
      yPos <- usr2[4]-r 
    rD <- abs((usr2[1] - usr2[2]))
    rf <- fig[1]
    r <- rD/rF * .4 # offset
    xPos <- tks1[1] - r

    if (srt.xlabel < 90)
      text(xPos, yPos, labels=x, adj = c(0, 0), srt=srt.xlabel, cex = cex.xlabel, xpd=NA)
    else
      text(xPos, yPos, labels=x, adj = c(0, 1), srt=srt.xlabel, cex = cex.xlabel, xpd=NA)
}

shiny_running = function () {
  # Look for runApp call somewhere in the call stack.
  # from https://stackoverflow.com/questions/32806974/detecting-whether-shiny-runs-the-r-code
  
  frames = sys.frames()
  calls = lapply(sys.calls(), `[[`, 1)
  call_name = function (call)
    if (is.function(call)) '<closure>' else deparse(call)
  call_names = vapply(calls, call_name, character(1))
  
  #  target_call = grep('^runApp$', call_names)
  target_call = grep('runApp$', call_names)
  
  if (length(target_call) == 0)
    return(FALSE)
  
  # Found a function called runApp, verify that it's Shiny's.
  target_frame = frames[[target_call]]
  namespace_frame = parent.env(target_frame)
  isNamespace(namespace_frame) && environmentName(namespace_frame) == 'shiny'
}

plotCumul <- function(x, y, i, nm) 
{
  nG <- ncol(groupData)
  groupN <- as.integer(colnames(groupData))
  N <- length(x)
  usr <- par("usr")
  segments(usr[2], usr[3], usr[2], usr[4], col="grey")
  for (j in nG:1) {
    y2 <- c(y[1], y, y[N])
    x2 <- c(usr[1], groupData[, j, drop=TRUE], usr[1])
    bord <- NA
    if (exists("cumulLine") & exists("cumulLineCol") & cumulLine)
       bord <- cumulLineCol
    polygon(x2, y2, col=groupColours[groupN[j]], border = bord, xpd=FALSE)
  } 
  fig <- par("fig")
  oldmar <- par("mar")
  par(mar=c(0, 0, 0, 0))
  oldusr <- par("usr")
  oldfig <- fig
  fig[3] <- fig[4]
  fig[4] <- 1
  par(fig=fig, new=TRUE)
  plot(0, xlim=c(0,1), ylim=c(0, 1), axes=FALSE, type="n", xlab="", ylab="", xaxs="i", yaxs="i")
  fin <- par("fin")
  scale <- 1.0 / fin[2] 
  lineHeight_in <- strheight("M", units="figure", cex=groupCex) * fin[2]
  inc <- min(0.25, lineHeight_in) * scale
  for (i in nG:1) {
     y <- (nG-i)*inc*1.3 + (0.1* scale)
     rect(0.8, y, 1.0, y+inc, col=groupColours[groupN[i]])
     if (!is.null(groupNames)) {
        text(0.75, y+inc/2, groupNames[i], adj=c(1, 0.5), cex=groupCex)
     }
  }
  par(mar=oldmar)
  par(fig=oldfig)
  par(usr=oldusr)
}

addZone2 <- function(rp, upper, lower=NULL, ...) {
  oldpar <- par(c("fig", "mar", "usr"))
  par(fig=rp$box)
  par(mar=c(0,0,0,0))
  par(usr=c(0, 1, rp$usr[3], rp$usr[4]))
  if (is.null(lower))
    segments(0, upper, 1, upper, xpd=NA, ...)
  else
    rect(0, lower, 1, upper, ...)
  par(oldpar)
}
