#############################################
#
#  Initialize PT conditions
#
############################################
cat("Initializing PT conditions...\n")
# fix-tag: error handling
input_pt <- rep(list(rep(list(c(Pressure = NULL, Temperature = NULL)), x_n)), y_n)
#############################################
#
#  PT definition from configuration file
#
############################################
if (pt_def == "input") {
  cat("Calculating PT conditions from inputs...\n")
  pnts <- unlist(strsplit(names(pt_definitions), "_"))
  for (h in 1:length(pt_definitions)) {
    a <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", pnts[h * 2 - 1])), split = ";"))
    b <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", pnts[h * 2])), split = ";"))
    for (x_i in a[1]:b[1]) {
      for (y_i in a[2]:b[2]) {
        press <- eval(parse(text = paste0(pt_definitions[[h]][1])))
        temp <- eval(parse(text = paste0(pt_definitions[[h]][2])))
        input_pt[[y_i]][[x_i]] <- matrix(c(press, temp), nrow = 1, dimnames = list("PT", c("Pressure", "Temperature")))
      }
    }
  }
}
#############################################
#
#  PT definition from input file
#
############################################
if (pt_def == "file") {
  # mod-tag: cleanup below
  cat("Reading PT conditions from file...\n")
  PT_file <- paste(work_dir, "/Projects/", workingfile, "/Inputs/PT.txt", sep = "")
  PT0 <- read.table(PT_file, sep = "\t")
  if (any(is.na(PT0[1, ]))) {
    warning("Error at init_pt.r\nError in PT file, NA found in column names, probably a tab left at the end of the first row, will try to rectify and continue")
    # mod-tag: Surely there's a better way to make a matrix numeric?
    make_numeric_mat <- function(x) {
      y <- matrix(as.numeric(as.vector(as.matrix(x))), nrow(x), )
      rownames(y) <- rownames(x)
      colnames(y) <- colnames(x)
      return(y)
    }
    err <- length(PT0[1, ])
    mat <- PT0[-1, -1]
    colnames(mat) <- as.character(unlist(PT0[1, 1:(err - 1)]))
    rownames(mat) <- 1:(nrow(PT0) - 1)
    PT0 <- make_numeric_mat(mat)
  }
  for (i in 1:nrow(PT0)) {
    PT[[PT0[i, 1]]][[PT0[i, 2]]]$press <- PT0[i, 3]
    PT[[PT0[i, 1]]][[PT0[i, 2]]]$temp <- PT0[i, 4]
  }
}
#############################################
#
#  Error validation
#
############################################
for (y_i in 1:y_n) {
  for (x_i in 1:x_n) {
    if (is.null(input_pt[[y_i]][[x_i]])) {
      cat("Error: No pt defined for x_i =", x_i, " y_i =", y_i, " \n")
      stop()
    }
  }
}
cat("Done with PT conditions\n")
cat("..............................................\n")
