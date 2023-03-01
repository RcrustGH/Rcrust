#############################################
#
#  Initialize bulk composition data (only major elements for now)
#
############################################
cat("Initializing bulk composition...\n")
if (exists("calculate_traces")) {
  if (calculate_traces) {
    if (!require(GCDmodel)) {
      # try to install GCDmodel
      install.packages(paste0(gsub("/code", "/plugins", getwd()), "/GCDmodel_0.5.02.zip"))
      if (!require(GCDmodel)) {
        for (i in 1:length(.libPaths())) {
          # if install was unsuccessful, try to manually place package in library/libraries
          try(
            unzip(paste0(gsub("/code", "/plugins", getwd()), "/GCDmodel.zip"),
              exdir = paste0(.libPaths()[[i]], ),
              overwrite = TRUE
            ),
            silent = TRUE
          )
          if (dir.exists(paste0(.libPaths()[[i]], "/GCDmodel"))) {
            cat("\nGCDmodel manually unpacked at ", paste0(.libPaths()[[i]], "/GCDmodel \n"))
            break
          }
        }
      }
    }
    library(GCDmodel)
    if (apply_trace_correction == "Apatite saturation") {
      if (!kd_file == "") {
        kd.in <- read.table(paste0(gsub("/code", "/data/", getwd()), kd_file), sep = "\t")
      }
    } else {
      kd.in <- read.table(paste0(gsub("/code", "/data/", getwd()), kd_file), sep = "\t")
    }
  }
}
# Error handling
if (exists("bulk_def")) {
  if (bulk_def == "input") {
    cat("Bulk composition defined from inputs...\n")
  } else {
    if (bulk_def == "file") {
      cat("Bulk composition defined from file...\n")
    } else {
      cat("Error ! No valid option specified for major composition...\n")
      stop()
    }
  }
} else {
  cat("Error ! bulk_def not specified...\n")
  stop()
}
# Create data structure
cat("P-T-X space under investigation with x =", x_n, "and y =", y_n, "\n")
input_bulk <- rep(list(rep(list(NULL), x_n)), y_n)
#############################################
#
#  Bulk definition from input
#
############################################
## Majors and traces
# temporary setting to FALSE
if (class(try(set_oxygen_fugacity, silent = TRUE)) == "try-error") {
  set_oxygen_fugacity <- FALSE
}
if (set_oxygen_fugacity) {
  # modtag log10(fugacity) is labelled here as O2 but input as log10(fugacity)
  all_elements <- c(major_elements, trace_elements, "O2")
} else {
  all_elements <- c(major_elements, trace_elements)
}
all_elements <- setdiff(all_elements, "")

if (bulk_def == "input") {
  cat("Creating bulk compositions from definitions in configuration file\n")
  # fix-tag should set Mass and volume somewhere but can only normalise once expressions have been evaluated before parsing to Rcrust()
  # c0_normalised<-c(bulk_definitions[[h]]/sum(bulk_definitions[[h]])*100,100,100)
  # input_bulk[[y_i]][[x_i]]<-matrix(c0_normalised,nrow=1,dimnames = list("Bulk",c(all_elements,"mass","volume")))
  # all_elements<-toupper(all_elements)
  pnts <- unlist(strsplit(names(bulk_definitions), "_"))
  for (h in 1:(length(bulk_definitions))) {
    a <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", pnts[h * 2 - 1])), split = ";"))
    b <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", pnts[h * 2])), split = ";"))
    for (x_i in a[1]:b[1]) {
      for (y_i in a[2]:b[2]) {
        if (length(bulk_definitions[[h]]) == 1) {
          input_bulk[[y_i]][[x_i]] <- rep(bulk_definitions[[h]], length(all_elements) + 1)
        } else {
          input_bulk[[y_i]][[x_i]] <- c(bulk_definitions[[h]])
        }
        if (!length(c(all_elements, "mass")) == length(input_bulk[[y_i]][[x_i]])) {
          cat(paste0("Error number of entries in input_bulk (", length(input_bulk[[y_i]][[x_i]]), ") is not equal to all elements + 1 (", length(c(all_elements, "mass")), ") for x_i = ", x_i, "; y_i = ", y_i, "\n"))
          stop()
        }
        names(input_bulk[[y_i]][[x_i]]) <- c(all_elements, "mass")
        # Apply Scaling on x_i - Note if assigning dependence do not leave gaps while scaling
        # mod-tag: work this into GUI
        scale_bulk_x_i <- FALSE
        if (scale_bulk_x_i) {
          if (x_i > 1) {
            leftmost <- 0
            for (i in 2:(x_i - 1)) {
              if (is.null(input_bulk[[y_i]][[i]])) {
                leftmost <- i
                break
              }
            }
            incr <- 1
            for (i in leftmost:(x_i - 1)) {
              # scale proportionate value to distance
              b_prop <- incr / (length(leftmost:(x_i - 1)) + 1)
              input_bulk[[y_i]][[i]] <- as.numeric(input_bulk[[y_i]][[x_i]]) * b_prop + as.numeric(input_bulk[[y_i]][[leftmost - 1]]) * (1 - b_prop)
              names(input_bulk[[y_i]][[i]]) <- c(all_elements, "mass")
              class(input_bulk[[y_i]][[i]]) <- "character"
              incr <- incr + 1
            }
          }
        }
      }
    }
  }
}
#############################################
#
#  Bulk definition from file
#
############################################
## Majors
if (bulk_def == "file") {
  cat("Bulk composition for majors defined from bulk file:", bulk_file, "\n")
  bulk_file_loc <- paste0(projects_directory, "/", working_file, "/Inputs/", bulk_file)
  if (file.exists(bulk_file_loc)) {
    table_in <- as.matrix(read.table(bulk_file_loc, sep = "\t"))
    table_out <- table_in[-1, , drop = FALSE]
    colnames(table_out) <- table_in[1, ]
    all_elements <- setdiff(table_in[1, ], c("from", "to", "mass"))
    for (h in 1:nrow(table_out)) {
      a <- unlist(strsplit(table_out[h, 1], split = ";"))
      b <- unlist(strsplit(table_out[h, 2], split = ";"))
      for (x_i in a[1]:b[1]) {
        for (y_i in a[2]:b[2]) {
          input_bulk[[y_i]][[x_i]] <- table_out[h, c(-1, -2)]
          names(input_bulk[[y_i]][[x_i]]) <- c(all_elements, "mass")
        }
      }
    }
  } else {
    cat("Error during bulk initialisation\nBulk file: ", bulk_file_loc, " not found")
    stop()
  }
}
#############################################
#
#  Error validation
#
############################################
for (y_i in 1:y_n) {
  for (x_i in 1:x_n) {
    if (is.null(input_bulk[[y_i]][[x_i]])) {
      cat("Error: No bulk defined for x_i =", x_i, " y_i =", y_i, " \n")
      stop()
    }
  }
}
cat("Done with bulk composition preparation\n")
cat("..............................................\n")
