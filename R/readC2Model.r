read.C2Model <- function(fName) {
  if (!file.exists(fName)) 
      stop(paste("Cannot open file", fName))
  if(.Platform$OS.type == "windows" & .Machine$sizeof.pointer < 5) {
#     if (require(RODBC, quietly=TRUE)==FALSE) {
#        stop("This function requires package RODBC")
#     }
     on.exit(odbcCloseAll())
#     channel <- odbcConnectExcel(fName)
#     fp <- RODBC:::full.path(fName)
#     con <- paste("Driver={Microsoft Excel Driver (*.xls)};DriverId=790;Dbq=", fp, ";DefaultDir=", dirname(fp), ";", sep = "")
#     channel <- odbcDriverConnect(con, tabQuote = c("[", "]"))
     channel <- odbcConnectExcel(fName, tabQuote = c("[", "]"))
     tabs <- sqlTables(channel, errors=TRUE)
     tabs <- tabs[tabs$TABLE_TYPE == "TABLE", ]
     ntabs <- nrow(tabs)
     X <- vector("list", length=ntabs)
     names(X) <- tabs$TABLE_NAME
     for (i in 1:ntabs) {
        X1 <- sqlFetch(channel, tabs$TABLE_NAME[i], colnames=FALSE, rownames=FALSE)
        if (tabs$TABLE_NAME[i] == "Summary") {
           X[[i]] <- as.character(X1[!is.na(X1), 1])
        } else {
           rownames(X1) <- X1[, 2]
           X[[i]] <- X1[, -c(1:3), drop=FALSE]
        }
     }
     odbcCloseAll()
  } else {
#     if (require(gdata, quietly=TRUE)==FALSE) {
#        stop("This function requires package gdata")
#     }
     X <- list()
     nSheet <- 1
     while (TRUE) {
        cat(paste("Reading sheet", nSheet, "\n"))
        ret <- try(read.xls(fName, nSheet, verbose=FALSE), silent=TRUE)
        if (class(ret) == "try-error") {
           break
        } else {
          X[[nSheet]] <- ret
          names(X)[nSheet] <- paste("Sheet", nSheet, sep="")
          nSheet <- nSheet + 1
        }
        flush.console()
     }
  }
  class(X) <- "C2"
  X
}

print.C2 <- function(x, ...) {
   if ("Summary" %in% names(x)) {
     n <- grep(":", x$Summary)
     cat(as.character(x$Summary[1:max(n)]), sep="\n")  
   } else {
     n <- grep(":", x[[1]][, 1])
     cat(as.character(x[[1]][1:max(n), 1]), sep="\n")  
   }
   cat("\nModel contains the following components:\n\n")
   cat(names(x), sep="\n")
}

summary.C2 <- function(object, ...) {
   if ("Summary" %in% names(object))
     cat(as.character(object$Summary), sep="\n")  
   else 
     cat(as.character(object[[1]][, 1]), sep="\n")  
}


