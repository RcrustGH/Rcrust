###############################
## Rcrust (server.r)
###############################
##
# Functions
##
# fix-tag remove working directory and projects directory and replace with gsub argument to auto locate relative to opened file
# function-def: exists_and_numeric(x)
exists_and_numeric <- function(x) {
  chk <- try(x, silent = TRUE)
  if (class(chk) == "try-error") {
    return(FALSE)
  } else {
    if (is.na(as.numeric(chk))) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}
# function-def: get_val(y_i,x_i,phase,variable,crust=crust,oxy_num=24,site_ocup="biotite",select=1)
get_val <- function(y_i, x_i, phase, variable, crust = crust, oxy_num = 24, site_ocup = "biotite", select = 1) {
  if (variable == "Temperature" | variable == "Pressure") {
    chk <- try(input_pt[[y_i]][[x_i]][which(colnames(input_pt[[1]][[1]]) == variable)], silent = TRUE)
    if (class(chk) == "try-error") {
      val <- 0
    } else {
      val <- input_pt[[y_i]][[x_i]][which(colnames(input_pt[[1]][[1]]) == variable)]
    }
    return(val)
  }
  if (variable == "Min_formula") {
    return(Min_formula(y_i, x_i, phase, oxy_num, site_ocup, crust = crust)[as.numeric(select)])
  }
  chk <- try(crust[[y_i]][[x_i]][phase, variable], silent = TRUE)
  if (class(chk) == "try-error") {
    val <- 0
  } else {
    if (class(chk) == "NULL") {
      val <- NA
    } else {
      val <- crust[[y_i]][[x_i]][phase, variable]
    }
  }
  return(val)
}
# function-def break_on_comma(x)
# Break on comma unless enclosed in quotes
break_on_comma <- function(x) {
  # find commas
  commas <- which(strsplit(x, "")[[1]] == ",")
  # name commas
  names(commas) <- rep("break", length.out = length(commas))
  # find quotes
  quotes <- which(strsplit(x, "")[[1]] == "\"")
  # number quotes
  names(quotes) <- rep(c("odd", "even"), length.out = length(quotes))
  # if comma is between an odd and an even quotes set to not break
  if (length(quotes > 0)) {
    new_names <- NULL
    for (comma in commas) {
      new_name <- "break"
      if (comma > min(quotes) & comma < max(quotes)) {
        below <- which(quotes < comma)[length(which(quotes < comma))]
        above <- which(quotes > comma)[1]
        if (names(below) == "odd" & names(above) == "even") {
          new_name <- "dont"
        }
      }
      new_names <- c(new_names, new_name)
    }
    names(commas) <- new_names
  }
  # break by commas tagged as break
  if (any(names(commas) == "break")) {
    breaks <- c(-1, commas[which(names(commas) == "break")] - 1, nchar(x))
    new_x <- NULL
    for (i in 1:(length(breaks) - 1)) {
      new_x <- c(new_x, substr(x, breaks[i] + 2, breaks[i + 1]))
    }
    x <- new_x
  }
  return(x)
}
# function-def:sub_brackets(x)
# substitute brackets for obrac and cbrac
sub_brackets <- function(x) {
  return(gsub("\\(", "_obrac_", gsub("\\)", "_cbrac_", x)))
}
# function-def:ret_brackets(x)
# return brackets from place holders obrac and cbrac
ret_brackets <- function(x) {
  return(gsub("_obrac_", "\\(", gsub("_cbrac_", "\\)", x)))
}
# function-def:.wtd.add(thelines,prop="mass",avname="Averaged")
.wtd.add <<- function(thelines, prop = "mass", avname = "Averaged") {
  if (nrow(thelines) == 1) {
    foo <- thelines
  } else {
    # fix-tag: is N(g) extensive or intensive
    extensive.cn <- c("wt%", "vol%", "mol%", "mol", "mass")
    extensive.cn <- intersect(extensive.cn, colnames(thelines))
    intensive.cn <- setdiff(colnames(thelines), c(prop, extensive.cn, "Phase"))
    foo <- matrix(rep(0, length(colnames(thelines))), nrow = 1)
    colnames(foo) <- colnames(thelines)
    suppressWarnings(class(thelines) <- "numeric")
    # Intensive
    foo[, intensive.cn] <- thelines[, prop] %*% thelines[, intensive.cn, drop = F] / sum(thelines[, prop])
    # Extensive
    for (i in extensive.cn) {
      foo[1, i] <- sum(thelines[, i])
    }
  }
  rownames(foo) <- avname
  return(foo)
}
# function-def:check_tuple(tuple)
check_tuple <- function(tuple = "") {
  if (tuple == "") {
    return(list("Valid tuple", ""))
  }
  split_tuple <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", tuple)), split = ",|;"))
  if (!length(split_tuple) == 2) {
    return(list("Error: Tuple is invalid", NULL))
  }
  for (k in 1:length(split_tuple)) {
    if (is.na(suppressWarnings(as.numeric(split_tuple[k])))) {
      return(list("Error: Tuple argument is non-numeric", NULL))
    }
    if (as.numeric(split_tuple[k]) < 1) {
      return(list("Error: Tuple argument is less than 1", NULL))
    }
    if (!as.numeric(split_tuple[k]) %% 1 == 0) {
      return(list("Error: Tuple argument is not a whole number", NULL))
    }
  }
  return(list("Valid tuple", paste0("{", paste(split_tuple, collapse = ";"), "}")))
}
# function-def:flip_y(mat)
# matrix reflection in y direction
flip_y <- function(mat) {
  mat_out <- mat[nrow(mat):1, , drop = FALSE]
  return(mat_out)
}
# function-def:flip_x(mat)
# matrix reflection in x direction
flip_x <- function(mat) {
  mat_out <- mat[, ncol(mat):1, drop = FALSE]
  return(mat_out)
}
# function-def:rotate(x)
# Matrix rotation
rotate <- function(x) t(apply(x, 2, rev))
substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}
# function-def:phase_abundance(crust,axis,path=1,p_a=1,p_b=p_a,path_label="Point",input_pt=NULL)
phase_abundance <- function(crust, axis, path = 1, p_a = 1, p_b = p_a, path_label = "Point", input_pt = NULL, proportion = "mass") {
  # Header
  outname <- switch(axis,
    x = paste0("Phase adundance vs ", path_label, " for {", p_a, ";", path, "} to {", p_b, ";", path, "}"),
    y = paste0("Phase adundance vs ", path_label, " for {", path, ";", p_a, "} to {", path, ";", p_b, "}")
  )
  # Abundance rows
  switch(axis,
    x = {
      pnt_y <- "path"
      pnt_x <- "p_i"
    },
    y = {
      pnt_y <- "p_i"
      pnt_x <- "path"
    }
  )
  all_phases <- NULL
  for (p_i in p_a:p_b) {
    all_phases <- union(all_phases, rownames(crust[[eval(parse(text = pnt_y))]][[eval(parse(text = pnt_x))]]))
  }
  # reorder into rs,es,es_cumul,as
  rs_phases <- all_phases[grepl(".*_rs", all_phases)]
  rs_bulk <- grep("Bulk_rs", rs_phases)
  es_phases <- all_phases[grepl(".*_es", all_phases)]
  es_cumul <- grep(".*_es_cumul", es_phases)
  as_phases <- all_phases[grepl(".*_as", all_phases)]
  all_phases <- c(rs_phases[-rs_bulk], rs_phases[rs_bulk], es_phases[-es_cumul], es_phases[es_cumul], as_phases)
  abundance_rows <- NULL
  for (ph in all_phases) {
    abundace_phase <- NULL
    for (p_i in p_a:p_b) {
      chk_phase <- try(crust[[eval(parse(text = pnt_y))]][[eval(parse(text = pnt_x))]][ph, proportion], silent = TRUE)
      if (class(chk_phase) == "try-error") {
        abundace_phase <- c(abundace_phase, 0)
      } else {
        abundace_phase <- c(abundace_phase, crust[[eval(parse(text = pnt_y))]][[eval(parse(text = pnt_x))]][ph, proportion])
      }
    }
    abundance_rows <- rbind(abundance_rows, matrix(abundace_phase, 1))
  }
  rownames(abundance_rows) <- all_phases
  col_nms <- p_a:p_b
  if (!path_label == "Point") {
    if (path_label == "Pressure(kbar)" | path_label == "Temperature(C)") {
      switch(path_label,
        "Pressure(kbar)" = slct <- 1,
        "Temperature(C)" = slct <- 2
      )
      if (!is.null(input_pt)) {
        col_nms <- NULL
        for (p_i in p_a:p_b) {
          col_nms <- c(col_nms, input_pt[[eval(parse(text = pnt_y))]][[eval(parse(text = pnt_x))]][slct])
        }
      }
    }
  }
  colnames(abundance_rows) <- col_nms
  return(list(outname, abundance_rows))
}
# function-def:write_phase_abundance(data,working_file=working_file,projects_directory=projects_directory,file_type)
write_phase_abundance <- function(data, working_file = working_file, projects_directory = projects_directory, file_type) {
  outfile_path <- paste0(projects_directory, "/", working_file, "/Outputs/", working_file, " ", data[[1]], file_type)
  if (file_type == ".txt") {
    write.table(data[[2]], outfile_path, sep = "\t", quote = F, row.names = TRUE)
  }
  if (file_type == ".csv") {
    write_test <- try(write.csv(data[[2]], outfile_path, row.names = TRUE), silent = TRUE)
    if (class(write_test) == "try-error") {
      cat("Error cannot write to ", outfile_path, ", please close all programs that may be accessing the file then try again\n")
      return(paste0("Error could not save phase abundance file: ", outfile_path, ", file may be open in another program, please
	  close all programs that may be accessing the file then try again\n"))
    }
  }
  cat("File written to ", data[[1]], "\n")
  return(paste0("Phase abundance file saved to ", projects_directory, "/", working_file, "/Outputs/\n"))
}
# function-def:get_PAM_names(crust,PAM_system)
get_PAM_names <- function(crust, PAM_system) {
  y_n <- length(crust)
  x_n <- length(crust[[1]])
  # Group phases
  PAM <- matrix(0, y_n, x_n)
  for (x_i in 1:x_n) {
    for (y_i in 1:y_n) {
      # Grab phases from system under consideration
      if (PAM_system == "Reactive Subsystem") {
        grab <- rownames(crust[[y_i]][[x_i]])[which(substrRight(rownames(crust[[y_i]][[x_i]]), 3) == "_rs")]
      }
      if (PAM_system == "Extract Subsystem") {
        grab <- rownames(crust[[y_i]][[x_i]])[which(substrRight(rownames(crust[[y_i]][[x_i]]), 3) == "_es")]
      }
      if (PAM_system == "Full System") {
        grab <- rownames(crust[[y_i]][[x_i]])
      }
      if (length(grab) > 0) {
        # Remove bulk arguments
        grab <- grab[-which(grab == "Bulk_rs" | grab == "Bulk_es" | grab == "Bulk_es_cumul")]
        if (!PAM_system == "Full System") {
          # remove identifiers
          for (i in 1:length(grab)) {
            grab[i] <- substr(grab[i], 1, nchar(grab[i]) - 3)
          }
        }
        # Reorder alphabetically
        grab <- sort(grab)
      }
      collapsed <- paste(grab, collapse = "+")
      if (collapsed == "") {
        collapsed <- "No Phases"
      }
      PAM[y_i, x_i] <- collapsed
    }
  }
  PAM_names <- setdiff(intersect(PAM, PAM), "0")
  return(list(PAM_names, PAM))
}
# function-def:PAM_calc(crust,PAM_system,compile_PAM=FALSE,PAM_compilation=NULL)
PAM_calc <- function(crust, PAM_system, compile_PAM = FALSE, PAM_compilation = NULL) {
  PAM_data <- get_PAM_names(crust, PAM_system)
  PAM <- PAM_data[[2]]
  if (compile_PAM) {
    validate(need(
      file.exists(
        paste0(sub("/code", "/Projects", getwd()), "/Compile/", PAM_compilation, " compilation legend.txt")
      ),
      paste0(PAM_compilation, " compilation legend.txt not found in ", sub("/code", "/Projects", getwd()), "/Compile/", "\nPlease compile legend first")
    ))
    compilation_names <- read.table(paste0(sub("/code", "/Projects", getwd()), "/Compile/", PAM_compilation, " compilation legend.txt"))
    class(compilation_names) <- "vector"
    PAM_names <- as.character(compilation_names[[2]])
  } else {
    PAM_names <- PAM_data[[1]]
  }
  # Shade by number of phases (variance)
  no_phs <- NULL
  for (i in 1:length(PAM_names)) {
    no_phs <- c(no_phs, length(strsplit(PAM_names, split = "+", fixed = TRUE)[[i]]))
  }
  PAM_names <- PAM_names[rev(order(no_phs, PAM_names))]
  # Populate grid
  for (y_i in 1:length(crust)) {
    for (x_i in 1:length(crust[[1]])) {
      if (!PAM[y_i, x_i] == 0) {
        PAM[y_i, x_i] <- which(PAM_names == PAM[y_i, x_i])
      }
    }
  }
  mode(PAM) <- "numeric"
  # Pull_common phases from legend
  nam <- strsplit(PAM_names, split = "+", fixed = TRUE)
  inter <- nam[[1]]
  if (length(PAM_names) > 1) {
    for (i in 1:length(PAM_names)) {
      inter <- intersect(inter, nam[[i]])
    }
  }
  all_pres <- paste(inter, collapse = "+")
  internam <- list()
  for (i in 1:length(PAM_names)) {
    internam[[i]] <- setdiff(nam[[i]], inter)
  }
  new_nam <- list()
  for (i in 1:length(PAM_names)) {
    if (length(internam[[i]]) > 0) {
      new_nam[i] <- paste(internam[[i]], collapse = "+")
    } else {
      new_nam[i] <- "-"
    }
  }
  PAM_names <- unlist(new_nam)
  PAM_legend <- 1:length(PAM_names)
  names(PAM_legend) <- PAM_names
  # if compiling save individual legend (phases relevant to this specific section)
  compile_PAM_legend <- NULL
  if (compile_PAM) {
    compile_PAM_legend <- PAM_legend
    PAM_legend <- PAM_legend[sort((setdiff(intersect(PAM, PAM), "0")))]
  }
  # Flip y-axis for matrix drawing
  PAM <- flip_y(PAM)
  library(raster)
  library(rgeos)
  x <- raster::raster(PAM)
  pol <- raster::rasterToPolygons(x, dissolve = TRUE)
  # create label ids
  pol_ex <- extract(x, pol)
  pol_id <- NULL
  for (i in 1:length(pol_ex)) {
    pol_id <- c(pol_id, pol_ex[[i]][1])
  }
  return(list(PAM, PAM_legend, all_pres, pol, pol_id, compile_PAM_legend))
}
# function-def:data_file(crust,x_n=length(crust[[1]]),y_n=length(crust),choose_columns=NULL,choose_rows=NULL,choose_points="All")
data_file <- function(crust, x_n = length(crust[[1]]), y_n = length(crust), choose_columns = NULL, choose_rows = NULL, choose_points = "All", assign_label) {
  #######################################################################
  # Outputs select_data list
  # Settings for outputing data_file

  # choose_columns    options =   from crust ("wt%",comps,"mass","V(J/bar)","H(J)","Gruneisen_T","Ks(bar)","Mu(bar)","V0(km/s)","Vp(km/s)","Vs(km/s)","Vp/Vs","Rho(kg/m3)","Cp(J/K)","alpha(1/K)","beta(1/bar)","S(J/K)","N(g)","Cp/Cv")
  #                               from other ("Phase","y_i","x_i","Pressure(kbar)","Temperature(C)")
  #                               Brief
  #                               All
  #                   default = All
  # choose_columns<-c("All")
  # choose_columns<-c("Brief")

  # choose_rows       default = All
  # choose_rows       options =  Reactive subsystem
  #                             Extract subsystem
  #                             Addition subsystem
  # choose_rows<-c("All")
  # choose_rows<-c("Bio(TCC)_rs","Bulk_rs")

  # choose_points           default = All

  #######################################################################
  # fix-tag: number duplicate names (feldspar, mica)
  validate(need(!choose_points == "", "To select points enter arguments seperated by commas of the form {x_a;y_a} for single points or {x_a;y_a}:{x_b;y_b} for ranges where a<=b<=n"))
  if (choose_points == "All") {
    choose_points <- paste("{1;1}:{", x_n, ";", y_n, "}", sep = "")
  }
  choose_points <- unlist(strsplit(choose_points, split = ","))
  choose_points <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", choose_points)), split = ","))
  choose_points <- strsplit(choose_points, split = ":|;")
  data_out <- NULL
  for (i in 1:length(choose_points)) {
    x_a <- as.numeric(choose_points[[i]][1])
    y_a <- as.numeric(choose_points[[i]][2])
    if (length(choose_points[[i]]) > 2) {
      x_b <- as.numeric(choose_points[[i]][3])
      y_b <- as.numeric(choose_points[[i]][4])
    } else {
      x_b <- as.numeric(choose_points[[i]][1])
      y_b <- as.numeric(choose_points[[i]][2])
    }
    # error validation on selection (range must be possible i.e. b>=a,b<=n)
    validate(need(all(x_a <= x_b, y_a <= y_b, x_b <= x_n, y_b <= y_n), "To select points enter arguments seperated by commas of the form {x_a;y_a} for single points or {x_a;y_a}:{x_b;y_b} for ranges where a<=b<=n"))
    for (y_i in y_a:y_b) {
      for (x_i in x_a:x_b) {
        # fix-tag: create validation here for erronous pixels
        # error validation for point existence
        if (is.null(crust[[y_i]][[x_i]])) {
          # fix-tag: only works if point 1;1 is populated - dont know how many columns otherwise
          new_pnt <- matrix(NA, 1, (ncol(crust[[1]][[1]]) + 6))
          new_pnt[1, 1] <- paste("Blank", y_i, x_i, sep = "_")
        } else {
          ID <- matrix(paste(rownames(crust[[y_i]][[x_i]]), y_i, x_i, sep = "_"), ncol = 1)
          colnames(ID) <- "ID"
          phase <- matrix(rownames(crust[[y_i]][[x_i]]), ncol = 1)
          colnames(phase) <- "Phase"
          pnt <- matrix(c(y_i, x_i), nrow = nrow(crust[[y_i]][[x_i]]), ncol = 2, byrow = TRUE)
          colnames(pnt) <- c("y_i", "x_i")
          if (exists("input_pt")) {
            p_t <- matrix(input_pt[[y_i]][[x_i]], nrow = nrow(crust[[y_i]][[x_i]]), ncol = 2, byrow = TRUE)
          } else {
            p_t <- matrix(0, nrow = nrow(crust[[y_i]][[x_i]]), ncol = 2, byrow = TRUE)
          }
          colnames(p_t) <- c("Pressure(kbar)", "Temperature(C)")
          # new_pnt<-cbind(ID,phase,pnt,p_t,signif(crust[[y_i]][[x_i]],4))
          new_pnt <- cbind(ID, phase, pnt, p_t, crust[[y_i]][[x_i]])
          rownames(new_pnt) <- ID
        }
        data_out <- rbind(data_out, new_pnt)
      }
    }
  }
  if (is.null(choose_columns)) {
    choose_columns <- colnames(data_out)
  }
  if (any(choose_columns == "Brief")) {
    choose_columns <- union(
      choose_columns[-which(choose_columns == "Brief")],
      c("ID", "Phase", "y_i", "x_i", "Pressure(kbar)", "Temperature(C)", "wt%", comps, "mass")
    )
  }
  if (is.null(choose_rows)) {
    select_rows <- 1:nrow(data_out)
  } else {
    if (choose_rows == "All") {
      select_rows <- 1:nrow(data_out)
    } else {
      select_rows <- NULL
      if (!is.na(match("Reactive subsystem", choose_rows))) {
        rs_rows <- grep("_rs", rownames(data_out))
        choose_rows <- choose_rows[-match("Reactive subsystem", choose_rows)]
        select_rows <- union(select_rows, rs_rows)
      }
      if (!is.na(match("Extract subsystem", choose_rows))) {
        rs_rows <- grep("_es", rownames(data_out))
        choose_rows <- choose_rows[-match("Extract subsystem", choose_rows)]
        select_rows <- union(select_rows, rs_rows)
      }
      for (row_arg in choose_rows) {
        chk_names <- NULL
        for (i in 1:length(rownames(data_out))) {
          chk_names <- c(chk_names, paste(strsplit(rownames(data_out), "_")[[i]][c(-3, -4)], collapse = "_"))
        }
        select_rows <- union(select_rows, which(chk_names == row_arg))
      }
    }
  }
  data_out <- data_out[sort(select_rows), choose_columns, drop = FALSE]
  rownames(data_out) <- NULL
  return(data_out)
}
# function-def:write_data_file(data_out,working_file=working_file,projects_directory=projects_directory,file_type=".csv")
write_data_file <- function(data_out, working_file = working_file, projects_directory = projects_directory, file_type = ".csv") {
  # file_type           options = ".csv", ".txt"
  #                    default = ".csv"
  outfile_path <- paste0(projects_directory, "/", working_file, "/Outputs/", working_file, "_data_file", file_type)
  if (length(data_out) > 0) {
    if (file_type == ".txt") {
      write.table(data_out, outfile_path, sep = "\t", quote = F, row.names = FALSE)
      cat("File written to ", outfile_path, "\n")
    }
    if (file_type == ".csv") {
      write_test <- try(write.csv(data_out, outfile_path, row.names = FALSE), silent = TRUE)
      if (class(write_test) == "try-error") {
        cat("Error cannot write to ", outfile_path, ", please close all programs that may be accessing the file then try again\n")
        return(paste0("Error could not save Data File: file may be open in another program, please close all programs that may be accessing the file then try again\n"))
      }
      cat("File written to ", outfile_path, "\n")
    }
    if (file_type == ".ps") {
      cat("Error cannot write Data File to .ps format\n")
      return(paste0("Error cannot write Data File to .ps format\n"))
    }
  }
  return(paste0("Data File saved to ", paste0(projects_directory, "/", working_file, "/Outputs/"), "\n"))
}
# function-def:grid_data(Grid_variable,Grid_variable_phase="Bulk_rs",crust_in=crust_out(),input_pt_in=input_pt,oxy_num=24,site_ocup="biotite",select=1)
grid_data <- function(Grid_variable, Grid_variable_phase = "Bulk_rs", crust_in = crust_out(), input_pt_in = input_pt, oxy_num = 24, site_ocup = "biotite", select = 1) {
  x_n <- length(crust_in[[1]])
  y_n <- length(crust_in)
  if (!is.null(Grid_variable)) {
    grid_out_mat <- matrix(0, y_n, x_n)
    for (x_i in 1:x_n) {
      for (y_i in 1:y_n) {
        if (Grid_variable == "A/CNK") {
          if (any(rownames(crust_in[[y_i]][[x_i]]) == Grid_variable_phase)) {
            grid_out_mat[y_i, x_i] <- calcACNK(crust_in[[y_i]][[x_i]][Grid_variable_phase, ])
          } else {
            grid_out_mat[y_i, x_i] <- 0
          }
        } else {
          grid_out_mat[y_i, x_i] <- get_val(y_i, x_i, Grid_variable_phase, Grid_variable, crust_in, oxy_num, site_ocup, select)
        }
      }
    }
    # Flip matrix so origin is bottom left (psuedosection convention)
    grid_out_mat <- flip_y(grid_out_mat)
    if (Grid_variable == "y_i" | Grid_variable == "x_i" | Grid_variable == "Temperature" | Grid_variable == "Pressure") {
      grid_out_title <- Grid_variable
    } else {
      grid_out_title <- paste(Grid_variable_phase, Grid_variable)
    }
    suppressWarnings(mode(grid_out_mat) <- "numeric")
    return(list(grid_out_title, grid_out_mat))
  }
}
# function-def:.First()
# Function for shortcut
.First <- function() {
  cat("Loading Dependencies\n")
  load_dependencies()
  library(shiny)
  # Launch with GUI function
  Rcrust <<- function() {
    # If working directory is x\Projects\y then set to x\code
    if (length(grep("Rcrust/Projects/", getwd())) == 1) {
      setwd(paste0(strsplit(getwd(), split = "Projects")[[1]][1], "code"))
    }
    runApp()
  }
  # function-def:manual_load(working_file,projects_directory=paste0(substring(getwd(),1,nchar(getwd())-4),"Projects"))
  # Launch without GUI function
  manual_load <<- function(working_file, projects_directory = paste0(substring(getwd(), 1, nchar(getwd()) - 4), "Projects")) {
    source(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    source("main.r")
  }
  # launch the GUI
  # If working directory is x\Projects\y then set to x\code
  if (length(grep("Rcrust/Projects/", getwd())) == 1) {
    setwd(paste0(strsplit(getwd(), split = "Projects")[[1]][1], "code"))
  }
  runApp()
}
# function-def:Rcrust()
# Launch with GUI function
Rcrust <<- function() {
  # If working directory is x\Projects\y then set to x\code
  if (length(grep("Rcrust/Projects/", getwd())) == 1) {
    setwd(paste0(strsplit(getwd(), split = "Projects")[[1]][1], "code"))
  }
  runApp()
}
# Launch without GUI function
manual_load <<- function(working_file, projects_directory = paste0(substring(getwd(), 1, nchar(getwd()) - 4), "Projects")) {
  source(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
  source("main.r")
}
# function-def:error_handling(working_file,projects_directory)
# error handling on button press
error_handling <- function(working_file, projects_directory) {
  if (projects_directory == "") {
    return("Error: no projects directory specified")
  }
  if (!dir.exists(projects_directory)) {
    return("Error: projects directory does not exist")
  }
  if (working_file == "") {
    return("Error: no working file specified")
  }
  return("error handling passed")
}
# Validate inputs
# must stay in server.R
check_inputs <- function(input) {
  # GUI level validation check that have at least minimum of P,T,and X definitions in GUI
  var_missing <- NULL
  chk_variables <- c("input$x_n", "input$y_n", "input$n_pt_def", "input$n_bulk_def", "input$major_elements")
  for (i in rev(chk_variables)) {
    if (any(eval(parse(text = i)) == "")) {
      var_missing <- i
    }
  }
  if (!is.null(var_missing)) {
    # reactive_message$data <- paste0("Error: Cannot run calculation. Missing ", var_missing)
    return(paste0("Error: Cannot run calculation. Missing ", var_missing))
    # Check that P2O5 selected for certain saturation corrections.
  }
  if (!input$x_n == "") {
    if (is.na(suppressWarnings(as.numeric(input$x_n)))) {
      return("Error: X must be numeric")
    }
    if (as.numeric(input$x_n) < 1) {
      return("Error: X must be greater than 1")
    }
    if (!as.numeric(input$x_n) %% 1 == 0) {
      return("Error: X must be a whole number")
    }
  }
  if (!input$y_n == "") {
    if (is.na(suppressWarnings(as.numeric(input$y_n)))) {
      return("Error: Y must be numeric")
    }
    if (as.numeric(input$y_n) < 1) {
      return("Error: Y must be greater than 1")
    }
    if (!as.numeric(input$y_n) %% 1 == 0) {
      return("Error: Y must be a whole number")
    }
  }
}
check_pt_def <- function(input) {
  if (!input$n_pt_def == "") {
    # Error validation
    if (is.na(suppressWarnings(as.numeric(input$n_pt_def)))) {
      return("Error: Number of PT definitions must be numeric")
    }
    if (as.numeric(input$n_pt_def) < 1) {
      return("Error: Number of PT definitions must be greater than 0")
    }
    if (!as.numeric(input$n_pt_def) %% 1 == 0) {
      return("Error: Number of PT definitions must be a whole number")
    }
    list_pt <- NULL
    for (i in 1:as.numeric(input$n_pt_def)) {
      # check tuples
      from <- check_tuple(eval(parse(text = paste0("input$pt_from_", i))))
      if (!from[[1]] == "Valid tuple") {
        return(paste0("Error in PT Definition:     ", from[[1]]))
      }
      to <- check_tuple(eval(parse(text = paste0("input$pt_to_", i))))
      if (!to[[1]] == "Valid tuple") {
        return(paste0("Error in PT Definition:     ", to[[1]]))
      }
      list_pt <- c(
        list_pt,
        paste0(
          "\"", from[[2]], "_", to[[2]], "\"=c(",
          pasteq(eval(parse(text = paste0("input$pressure_", i)))),
          ",",
          pasteq(eval(parse(text = paste0("input$temperature_", i)))),
          ")"
        )
      )
    }
    pt_definitions <- paste0("list(", paste0(list_pt, collapse = ","), ")")
  } else {
    pt_definitions <- ""
  }
  return(pt_definitions)
}
check_comp_trans <- function(input) {
  if (!input$n_comp_trans == "") {
    # Error validation
    if (is.na(suppressWarnings(as.numeric(input$n_comp_trans)))) {
      return("Error: Number of Component transformations must be numeric")
    }
    if (as.numeric(input$n_comp_trans) < 1) {
      return("Error: Number of Component transformations must be greater than 0")
    }
    if (!as.numeric(input$n_comp_trans) %% 1 == 0) {
      return("Error: Number of Component transformations must be a whole number")
    }
    list_trans <- NULL
    for (i in 1:as.numeric(input$n_comp_trans)) {
      list_trans <- c(list_trans, paste0("\"", eval(parse(text = paste0("input$old_comp_", i))), "_", eval(parse(text = paste0("input$new_comp_", i))), "\"=c(", pasteq(eval(parse(text = paste0("input$comp_", i)))), ")"))
    }
    comp_transformations <- paste0("list(", paste0(list_trans, collapse = ","), ")")
  } else {
    comp_transformations <- ""
  }
  return(comp_transformations)
}
check_bulk_def <- function(input) {
  if (!input$n_bulk_def == "") {
    # Error validation
    if (is.na(suppressWarnings(as.numeric(input$n_bulk_def)))) {
      return("Error: Number of Bulk definitions must be numeric")
    }
    if (as.numeric(input$n_bulk_def) < 1) {
      return("Error: Number of Bulk definitions must be greater than 0")
    }
    if (!as.numeric(input$n_bulk_def) %% 1 == 0) {
      return("Error: Number of Bulk definitions must be a whole number")
    }
    list_bulk <- NULL
    for (i in 1:as.numeric(input$n_bulk_def)) {
      # check tuples
      from <- check_tuple(eval(parse(text = paste0("input$bulk_from_", i))))
      if (!from[[1]] == "Valid tuple") {
        return(paste0("Error in Bulk Definition:     ", from[[1]]))
      }
      to <- check_tuple(eval(parse(text = paste0("input$bulk_to_", i))))
      if (!to[[1]] == "Valid tuple") {
        return(paste0("Error in Bulk Definition:     ", to[[1]]))
      }
      bulk_eval <- strsplit(eval(parse(text = paste0("input$bulk_", i))), split = ",")[[1]]
      # place in quotes
      for (j in 1:length(bulk_eval)) {
        bulk_eval[j] <- pasteq(bulk_eval[j])
      }
      list_bulk <- c(list_bulk, paste0("\"", from[[2]], "_", to[[2]], "\"=c(", paste(bulk_eval, collapse = ","), ")"))
    }
    bulk_definitions <- paste0("list(", paste0(list_bulk, collapse = ","), ")")
  } else {
    bulk_definitions <- ""
  }
  return(bulk_definitions)
}
check_traces <- function(input) {
  if (input$apply_trace_correction == "Apatite saturation" || input$apply_trace_correction == "Apatite & Monazite Saturation") {
    if (!any(input$trace_elements == "P2O5")) {
      # reactive_message$data <- paste0("P2O5 not selected under trace elements for saturation corrections")
      return(paste0("P2O5 not selected under trace elements for saturation corrections"))
    }
    if (input$apply_trace_correction == "Apatite & Monazite Saturation") {
      if (!exists_and_numeric(input$D_ApMelt_LREE)) {
        # reactive_message$data <- paste0("Please enter a numeric value under D Ap/Melt LREE")
        return(paste0("Please enter a numeric value under D Ap/Melt LREE"))
      }
      if (!exists_and_numeric(input$Xmz)) {
        # reactive_message$data <- paste0("Please enter a numeric value under Xmz")
        return(paste0("Please enter a numeric value under Xmz"))
      }
    }
  }
  if (!input$apply_trace_correction == "Apatite saturation") {
    if (!file.exists(paste0(gsub("/code", "/data", getwd()), "/", input$kd_file))) {
      return(paste0("Error: Kd file not found at ", paste0(gsub("/code", "/data", getwd()), "/", input$kd_file)))
    }
  }
  if (any(is.na(input$trace_elements)) || is.null(input$trace_elements) || input$trace_elements == "") {
    # reactive_message$data <- paste0("You have selected Partition traces above solidus, but have not selected any trace elements.")
    return(paste0("You have selected Partition traces above solidus, but have not selected any trace elements."))
  }
  saturation_input <- ""
  if (input$apply_trace_correction == "Apatite saturation") {
    saturation_input <- paste0(
      saturation_input,
      "apatite_saturation_Ap<-", pasteq(input$apatite_saturation_Ap), "\n"
    )
  }
  if (input$apply_trace_correction == "Apatite & Monazite Saturation") {
    if (!exists_and_numeric(input$D_ApMelt_LREE)) {
      return(paste0("Error: Please enter a numeric value under D Ap/Melt LREE"))
    }
    if (!exists_and_numeric(input$Xmz)) {
      return(paste0("Error: Please enter a numeric value under Xmz"))
    }
    saturation_input <- paste0(
      saturation_input,
      "apatite_saturation_ApMnz<-", pasteq(input$apatite_saturation_ApMnz), "\n",
      "Xmz<-", paste0(input$Xmz), "\n",
      "D_ApMelt_LREE<-", paste0(input$D_ApMelt_LREE), "\n"
    )
  }
  return(saturation_input)
}
check_ph_add_def <- function(input) {
  if (!input$n_ph_add_def == "") {
    # Error validation
    if (is.na(suppressWarnings(as.numeric(input$n_ph_add_def)))) {
      return("Error: Number of Phase Addition definitions must be numeric")
    }
    if (as.numeric(input$n_ph_add_def) < 1) {
      return("Error: Number of Phase Addition definitions must be greater than 0")
    }
    if (!as.numeric(input$n_ph_add_def) %% 1 == 0) {
      return("Error: Number of Phase Addition definitions must be a whole number")
    }
    list_ph_add <- NULL
    for (i in 1:as.numeric(input$n_ph_add_def)) {
      ph_add_defs <- NULL
      # check tuples
      from <- check_tuple(eval(parse(text = paste0("input$ph_add_from_", i))))
      if (!from[[1]] == "Valid tuple") {
        return(paste0("Error in Phase Addition Definition:     ", from[[1]]))
      }
      to <- check_tuple(eval(parse(text = paste0("input$ph_add_to_", i))))
      if (!to[[1]] == "Valid tuple") {
        return(paste0("Error in Phase Addition Definition:     ", to[[1]]))
      }
      ph_add_con <- paste0("condition=", pasteq(eval(parse(text = paste0("input$ph_add_con_", i)))))
      phases <- unlist(strsplit(eval(parse(text = paste0("input$ph_add_phs_", i))), split = ","))
      for (j in 1:length(phases)) {
        ph_add_defs <- c(ph_add_defs, paste0(pasteq(phases[j]), "=", pasteq(eval(parse(text = paste0("input$ph_add_phs_", i, "_", phases[j]))))))
      }
      list_ph_add <- c(list_ph_add, paste0("\"", from[[2]], "_", to[[2]], "\"=c(", paste0(c(ph_add_con, ph_add_defs), collapse = ","), ")"))
    }
    ph_add_definitions <- paste0("list(", paste0(list_ph_add, collapse = ","), ")")
  } else {
    ph_add_definitions <- ""
  }
  return(ph_add_definitions)
}
check_ph_extr_def <- function(input) {
  if (!input$n_ph_extr_def == "") {
    # Error validation
    if (is.na(suppressWarnings(as.numeric(input$n_ph_extr_def)))) {
      return("Error: Number of Phase Extraction definitions must be numeric")
    }
    if (as.numeric(input$n_ph_extr_def) < 1) {
      return("Error: Number of Phase Extraction definitions must be greater than 0")
    }
    if (!as.numeric(input$n_ph_extr_def) %% 1 == 0) {
      return("Error: Number of Phase Extraction definitions must be a whole number")
    }
    list_ph_extr <- NULL
    for (i in 1:as.numeric(input$n_ph_extr_def)) {
      ph_extr_defs <- NULL
      # check tuples
      from <- check_tuple(eval(parse(text = paste0("input$ph_extr_from_", i))))
      if (!from[[1]] == "Valid tuple") {
        return(paste0("Error in Phase Extraction Definition:     ", from[[1]]))
      }
      to <- check_tuple(eval(parse(text = paste0("input$ph_extr_to_", i))))
      if (!to[[1]] == "Valid tuple") {
        return(paste0("Error in Phase Extraction Definition:     ", to[[1]]))
      }
      ph_extr_con <- paste0("condition=", pasteq(eval(parse(text = paste0("input$ph_extr_con_", i)))))
      # seperate on quotes then on comma
      phases <- gsub('"', "", break_on_comma(eval(parse(text = paste0("input$ph_extr_phs_", i)))))
      for (j in 1:length(phases)) {
        ph_extr_defs <- c(
          ph_extr_defs,
          paste0(
            pasteq(phases[j]),
            "=",
            pasteq(eval(parse(text = paste0("input$", '\"', "ph_extr_phs_", i, "_", sub_brackets(phases[j]), '\"'))))
          )
        )
      }
      list_ph_extr <- c(
        list_ph_extr,
        paste0(
          "\"", from[[2]], "_", to[[2]],
          "\"=c(",
          paste0(c(ph_extr_con, ph_extr_defs), collapse = ","),
          ")"
        )
      )
    }
    ph_extr_definitions <- paste0("list(", paste0(list_ph_extr, collapse = ","), ")")
  } else {
    ph_extr_definitions <- ""
  }
  return(ph_extr_definitions)
}
# loads component packet inputs
load_component_packet <- function(cp_components_r, cp_phases_r) {
  if (!all(cp_components_r == "")) {
    list_cp_isolate <- NULL
    cp_phases_combine <- NULL
    for (i in 1:length(cp_components_r)) {
      cp_phases_defs <- NULL
      # stores components and amounts to isolate from each
      list_cp_isolate <- c(list_cp_isolate, paste0("\"", names(cp_components_r[i]), "\"=\"", cp_components_r[i], "\""))
      # phases to partition component into
      for (j in 1:length(cp_phases_r[[i]])) {
        # phases and corresponding values
        cp_phases_defs <- c(cp_phases_defs, paste0(pasteq(names(cp_phases_r[[i]][j])), "=", pasteq(cp_phases_r[[i]][j])))
      }
      # stores phases and amounts with respective components as cp_phases_[component]
      if (i > 1) {
        cp_phases_combine <- c(
          cp_phases_combine, "\n",
          paste0(
            "cp_phases_",
            names(cp_components_r[i]),
            "<-c(",
            paste0(cp_phases_defs, collapse = ","),
            ")"
          )
        )
      } else {
        cp_phases_combine <- c(
          cp_phases_combine,
          paste0("cp_phases_", names(cp_components_r[i]), "<-c(", paste0(cp_phases_defs, collapse = ","), ")")
        )
      }
    }
    # combined into a string to store in text file
    cp_packet_definitions <- paste0(cp_phases_combine, collapse = "")
    list_cp_isolate <- paste0(list_cp_isolate, collapse = ",")
  } else {
    cp_packet_definitions <- "cp_phases_<-c(\"\")"
    list_cp_isolate <- ""
  }
  return(list(cp_packet_definitions, list_cp_isolate))
}
get_component_packet_inputs <- function(input, available_components_r) {
  list_cp_isolate <- NULL
  cp_phases_combine <- NULL
  packet_components <- c(break_on_comma(input$cp_components))
  for (i in 1:as.numeric(length(packet_components))) {
    cp_phases_defs <- NULL
    # stores components and amounts to isolate from each
    list_cp_isolate <- c(
      list_cp_isolate,
      paste0("\"", packet_components[i], "\"=\"", eval(parse(text = paste0("input$cp_isolate_", i))), "\"")
    )
    # phases to partition component into
    phases <- gsub('"', "", break_on_comma(eval(parse(text = paste0("input$cp_phases_", i)))))
    for (j in 1:length(phases)) {
      # phases and corresponding values
      cp_phases_defs <- c(
        cp_phases_defs,
        paste0(
          pasteq(phases[j]),
          "=",
          pasteq(eval(parse(text = paste0("input$", '\"', "cp_phases_", i, "_", sub_brackets(phases[j]), '\"'))))
        )
      )
    }
    # stores phases and amounts with respective components as cp_phases_[component]
    if (i > 1) {
      cp_phases_combine <- c(
        cp_phases_combine,
        "\n",
        paste0(
          "cp_phases_",
          packet_components[i],
          "<-c(",
          paste0(cp_phases_defs, collapse = ","),
          ")"
        )
      )
    } else {
      cp_phases_combine <- c(
        cp_phases_combine,
        paste0(
          "cp_phases_",
          packet_components[i],
          "<-c(",
          paste0(cp_phases_defs, collapse = ","),
          ")"
        )
      )
    }
  }
  # combined into a string to store in text file
  cp_packet_definitions <- paste0(cp_phases_combine, collapse = "")
  list_cp_isolate <- paste0(list_cp_isolate, collapse = ",")
  # checking for new components added by component packet, not important for ui, important for running
  # Mod-tag: tried to make it only show up in text file when component created, but problem with components remaining loaded in workspace
  new_components <- "\nnew_components<-c(\""
  count_new_cp <- 0
  for (i in 1:length(packet_components)) {
    if (is.na(match(packet_components[i], available_components_r))) {
      count_new_cp <- count_new_cp + 1
      if (count_new_cp > 1) {
        new_components <- paste0(new_components, "\",\"", packet_components[i])
      } else {
        new_components <- paste0(new_components, packet_components[i])
      }
    }
  }
  cp_packet_definitions <- paste0(cp_packet_definitions, new_components, "\")")
  return(list(cp_packet_definitions, list_cp_isolate))
}
check_phases_to_rename <- function(input) {
  if (!input$n_phases_to_rename == "") {
    # Error validation
    if (is.na(suppressWarnings(as.numeric(input$n_phases_to_rename)))) {
      return("Error: Number of phases to rename must be numeric")
    }
    if (as.numeric(input$n_phases_to_rename) < 1) {
      return("Error: Number of phases to rename must be greater than 0")
    }
    if (!as.numeric(input$n_phases_to_rename) %% 1 == 0) {
      return("Error: Number of phases to rename must be a whole number")
    }
    list_phases_to_rename <- NULL
    for (i in 1:as.numeric(input$n_phases_to_rename)) {
      list_phases_to_rename <- c(list_phases_to_rename, paste0("\"", eval(parse(text = paste0("input$old_name_", i))), "~", eval(parse(text = paste0("input$new_name_", i))), "\"=c(", pasteq(eval(parse(text = paste0("input$condition_", i)))), ")"))
    }
    phases_to_rename <- paste0("list(", paste0(list_phases_to_rename, collapse = ","), ")")
  } else {
    phases_to_rename <- ""
  }
  return(phases_to_rename)
}
# function-def:bl(x)
# blank function : returns a value or quotes "" if blank
bl <- function(x) {
  if (is.null(x)) {
    return("\"\"")
  } else {
    if (x == "") {
      return("\"\"")
    } else {
      return(x)
    }
  }
}
# function-def:bl(str_to_phases)
# string to phases function : takes a string e.g. "melt(HP),Pl,q" and returns a vector of the phases e.g. "melt(HP)","Pl","q"
str_to_phases <- function(x) {
  x_list <- strsplit(x, ",")
  if (length(x_list) == 0) {
    x_out <- NULL
  } else {
    x_out <- x_list[[1]][1]
    if (length(x_list[[1]]) > 1) {
      for (mi in 2:length(x_list[[1]])) {
        x_out <- c(x_out, x_list[[1]][mi])
      }
    }
  }
  return(x_out)
}
# function-def:pasteq(x)
# pasteq fubction: collapse vector into comma seperated terms with quotes around them
pasteq <- function(x) {
  if (length(x) == 0) {
    paste0("\"", x, "\"", collapse = ",")
  } else {
    a <- ""
    y <- NULL
    for (i in 1:length(x)) {
      if (is.null(x[i])) {
        y <- c(y, "")
      } else {
        if (x[i] == a) {
          y <- c(y, "")
        } else {
          # if want to remove all quotes (dont know if I want this) then gsub('"',"",x)
          # remove enclosing quotes if present
          x_split <- strsplit(as.character(x[i]), split = "")[[1]]
          if (x_split[1] == "\"") {
            x_split[1] <- ""
          }
          if (x_split[length(x_split)] == "\"") {
            x_split[length(x_split)] <- ""
          }
          x[i] <- paste(x_split, collapse = "")

          y <- c(y, x[i])
        }
      }
    }
    paste0("\"", y, "\"", collapse = ",")
  }
}
##
# Shiny Server
##
shinyServer(function(input, output, session) {
  # Stop app when browser closed
  session$onSessionEnded(stopApp)
  # Global Variables (all variables that must be kept between load events needs to be in the vector keep_on_load)
  keep_on_load <<- c("from_clear_button", "from_copy", "first_load")
  from_clear_button <<- FALSE
  from_copy <<- FALSE
  # Reactive variables (globally accessed)
  # Reactive stores
  pt_definitions_r <- reactiveValues(data = "")
  bulk_definitions_r <- reactiveValues(data = NULL)
  ph_add_definitions_r <- reactiveValues(data = NULL)
  ph_extr_definitions_r <- reactiveValues(data = NULL)
  cp_components_r <- reactiveValues(data = NULL)
  cp_phases_r <- reactiveValues(data = NULL)
  solution_models_file_r <- reactiveValues(data = NULL)
  available_components_r <- reactiveValues()
  comp_transformations_r <- reactiveValues(data = NULL)
  current_components_r <- reactiveValues(data = NULL)
  all_elements_r <- reactiveValues(data = NULL)
  major_elements_r <- reactiveValues(data = NULL)
  trace_elements_r <- reactiveValues(data = NULL)
  phases_to_rename_r <- reactiveValues(data = NULL)
  # mod-tag - this is a quick fix
  use_sol_models_r <- reactiveValues(data = NULL)
  # store_r contains output data from runs
  store_r <- reactiveValues(crust_r = NULL, input_pt_r = NULL, input_bulk_r = NULL, all_elements_r = NULL, major_elements_r = NULL, trace_elements_r = NULL, phases_to_rename_r = NULL)
  # initialise a passing message for error checking and reporting
  reactive_message <- reactiveValues(data = NULL)
  traces_message <- reactiveValues(data = NULL)
  load_pt_r <- reactiveValues(data = NULL)
  # mod-tag load elements for backwards compatability
  if (!exists("major_elements")) {
    major_elements <- ""
  }
  # if(!is.null(major_elements)){if(major_elements[1]==""){major_elements<-NULL}}
  if (!exists("trace_elements")) {
    trace_elements <- ""
  }
  if (!exists("all_elements")) {
    all_elements <- c(major_elements, trace_elements)
    all_elements <- setdiff(all_elements, "")
  }
  # save data reactive
  save_data_file <- reactive({
    # Grab working file and projects directory
    working_file <- input$working_file
    projects_directory <- input$projects_directory
    # error handling
    reactive_message$data <- error_handling(working_file, projects_directory)
    if (reactive_message$data == "error handling passed") {
      # if error handling is passed
      # Return success message
      return(paste0("File saved to ", projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    } else {
      # Return error message
      return(reactive_message$data)
    }
  })
  # Save function
  on_save <- reactive({
    # Grab working file and projects directory
    working_file <- input$working_file
    projects_directory <- input$projects_directory
    # Replace projects directory if from copy
    if (!exists("from_copy")) {
      from_copy <<- FALSE
    }
    if (from_copy) {
      projects_directory <- copy_directory
      from_copy <<- FALSE
    }
    # error handling
    reactive_message$data <- error_handling(working_file, projects_directory)
    if (reactive_message$data == "error handling passed") {
      # if error handling is passed
      # Check basic GUI inputs
      reactive_message$data <- check_inputs(input)
      if (length(grep("Error", reactive_message$data)) != 0) {
        return(reactive_message$data)
      }
      # Start building the save file
      w_file <- paste0(
        "###############\n#\n#   Rcrust input file\n#\n###############\n# Location of project files\n",
        "working_file<-", pasteq(working_file), "\n",
        "projects_directory<-", pasteq(projects_directory), "\n"
      )
      w_size <- paste0(
        "###############\n#\n#   Size data\n#\n###############\n# number of points in x and y directions\n",
        "x_n<-", bl(input$x_n), "\n",
        "y_n<-", bl(input$y_n), "\n"
      )
      # Create pt_definitions
      if (input$n_pt_def == "load") {
        if (!all(pt_definitions_r$data == "")) {
          list_pt <- NULL
          for (i in 1:length(pt_definitions_r$data)) {
            list_pt <- c(list_pt, paste0("\"", names(pt_definitions_r$data)[i], "\"=c(", pasteq(pt_definitions_r$data[[i]]), ")"))
          }
          pt_definitions <- paste0("list(", paste0(list_pt, collapse = ","), ")")
        } else {
          pt_definitions <- ""
        }
      } else {
        pt_definitions <- paste0(check_pt_def(input))
        if (length(grep("Error", pt_definitions)) != 0) {
          return(pt_definitions)
        }
      }
      w_pt <- paste0(
        "###############\n#\n#   PT data\n#\n###############\n",
        "pt_def<-\"input\"                         #input,file\n",
        "pt_definitions<-", bl(pt_definitions), "\n"
      )
      # component transformations
      if (input$n_comp_trans == "load") {
        if (!all(comp_transformations_r$data == "")) {
          list_trans <- NULL
          for (i in 1:length(comp_transformations_r$data)) {
            list_trans <- c(list_trans, paste0("\"", names(comp_transformations_r$data)[i], "\"=c(", pasteq(comp_transformations_r$data[[i]]), ")"))
          }
          comp_transformations <- paste0("list(", paste0(list_trans, collapse = ","), ")")
        } else {
          comp_transformations <- ""
        }
      } else {
        comp_transformations <- check_comp_trans(input)
        if (length(grep("Error", comp_transformations)) != 0) {
          return(comp_transformations)
        }
      }
      # bulk definitions
      if (input$n_bulk_def == "load") {
        if (!all(bulk_definitions_r$data == "")) {
          list_bulk <- NULL
          for (i in 1:length(bulk_definitions_r$data)) {
            list_bulk <- c(list_bulk, paste0("\"", names(bulk_definitions_r$data)[i], "\"=c(", pasteq(bulk_definitions_r$data[[i]]), ")"))
          }
          bulk_definitions <- paste0("list(", paste0(list_bulk, collapse = ","), ")")
        } else {
          bulk_definitions <- ""
        }
      } else {
        bulk_definitions <- check_bulk_def(input)
        if (length(grep("Error", bulk_definitions)) != 0) {
          return(bulk_definitions)
        }
      }
      # major elements
      if (!is.null(input$major_elements)) {
        if (any(input$major_elements == "load")) {
          major_elements <- major_elements_r$data
        } else {
          major_elements <- unlist(strsplit(input$major_elements, split = ","))
        }
      } else {
        major_elements <- ""
      }
      # trace elements
      if (!is.null(input$trace_elements)) {
        if (any(input$trace_elements == "load")) {
          trace_elements <- trace_elements_r$data
        } else {
          trace_elements <- unlist(strsplit(input$trace_elements, split = ","))
        }
      } else {
        trace_elements <- ""
      }
      # Accessory phase saturation
      if (input$Xmz == "" || is.null(input$Xmz)) {
        Xmz <- 0.83
      } else {
        Xmz <- input$Xmz
      }
      saturation_input <- ""
      if (input$calculate_traces) {
        saturation_input <- check_traces(input)
        if (length(grep("Error", saturation_input)) != 0) {
          return(saturation_input)
        }
      }
      all_elements <- c(major_elements, trace_elements)
      if (input$bulk_def_file == FALSE) {
        bulk_def <- "input"
      } else {
        bulk_def <- "file"
      }
      w_bulk_composition <- paste0(
        "###############\n#\n#   Bulk composition data \n#\n###############\n",
        "comp_transformations<-c(", bl(comp_transformations), ")\n",
        "bulk_def<-", pasteq(bulk_def), "                         #input,file\n",
        "major_elements<-c(", pasteq(major_elements), ")\n",
        "set_oxygen_fugacity<-", input$set_oxygen_fugacity, "\n",
        "calculate_traces<-", input$calculate_traces, "\n",
        "apply_trace_correction<-", pasteq(input$apply_trace_correction), "\n",
        saturation_input,
        "kd_file<-", pasteq(input$kd_file), "\n",
        "trace_elements<-c(", pasteq(trace_elements), ")\n",
        "bulk_definitions<-c(", bl(bulk_definitions), ")\n",
        "bulk_file<-", pasteq(input$bulk_file), "\n"
      )
      # phase additions
      if (input$ph_add) {
        if (input$n_ph_add_def == "load") {
          if (!all(ph_add_definitions_r$data == "")) {
            list_ph_add <- NULL
            for (i in 1:length(ph_add_definitions_r$data)) {
              qt_names <- unlist(lapply(names(ph_add_definitions_r$data[[i]]), pasteq))
              qt_values <- unlist(lapply(ph_add_definitions_r$data[[i]], pasteq))
              list_ph_add <- c(list_ph_add, paste0("\"", names(ph_add_definitions_r$data)[i], "\"=c(", paste(qt_names, qt_values, collapse = ",", sep = "="), ")"))
            }
            ph_add_definitions <- paste0("list(", paste0(list_ph_add, collapse = ","), ")")
          } else {
            ph_add_definitions <- ""
          }
        } else {
          ph_add_definitions <- check_ph_add_def(input)
          if (length(grep("Error", ph_add_definitions)) != 0) {
            return(ph_add_definitions)
          }
        }
      } else {
        ph_add_definitions <- ""
      }
      w_phase_addition <- paste0(
        "###############\n#\n#   Phase addition\n#\n###############\n",
        "ph_add<-", input$ph_add, "\n",
        "ph_add_definitions<-c(", bl(ph_add_definitions), ")\n"
      )
      # phase extractions
      if (input$ph_extr) {
        if (input$n_ph_extr_def == "load") {
          if (!all(ph_extr_definitions_r$data == "")) {
            list_ph_extr <- NULL
            for (i in 1:length(ph_extr_definitions_r$data)) {
              qt_names <- unlist(lapply(names(ph_extr_definitions_r$data[[i]]), pasteq))
              qt_values <- unlist(lapply(ph_extr_definitions_r$data[[i]], pasteq))
              list_ph_extr <- c(list_ph_extr, paste0("\"", names(ph_extr_definitions_r$data)[i], "\"=c(", paste(qt_names, qt_values, collapse = ",", sep = "="), ")"))
            }
            ph_extr_definitions <- paste0("list(", paste0(list_ph_extr, collapse = ","), ")")
          } else {
            ph_extr_definitions <- ""
          }
        } else {
          ph_extr_definitions <- check_ph_extr_def(input)
          if (length(grep("Error", ph_extr_definitions)) != 0) {
            return(ph_extr_definitions)
          }
        }
      } else {
        ph_extr_definitions <- ""
      }
      w_phase_extraction <- paste0(
        "###############\n#\n#   Phase extraction\n#\n###############\n",
        "ph_extr<-", input$ph_extr, "\n",
        "reequilibrate_steps<-", input$reequilibrate_steps, "\n",
        "ph_extr_definitions<-c(", bl(ph_extr_definitions), ")\n"
      )
      # component packet
      if (input$component_packet) {
        if (input$cp_components == "load") {
          # load inputs and store
          cp_loaded <- load_component_packet(cp_components_r$data, cp_phases_r$data)
          cp_packet_definitions <- cp_loaded[[1]]
          list_cp_isolate <- cp_loaded[[2]]
        } else {
          if (!isTRUE(input$cp_components == "")) {
            # save inputs
            cp_loaded <- get_component_packet_inputs(input, available_components_r$current)
            cp_packet_definitions <- cp_loaded[[1]]
            list_cp_isolate <- cp_loaded[[2]]
          } else {
            cp_packet_definitions <- "cp_phases_<-c(\"\")"
            list_cp_isolate <- ""
          }
        }
      } else {
        cp_packet_definitions <- "cp_phases_<-c(\"\")"
        list_cp_isolate <- ""
        new_components <- ""
      }
      w_component_packet <- paste0(
        "###############\n#\n#   Component packet\n#\n###############\n",
        "component_packet<-", input$component_packet, "\n",
        "cp_components<-c(", bl(list_cp_isolate), ")\n",
        bl(cp_packet_definitions), "\n"
      )
      # phase renaming
      if (input$n_phases_to_rename == "load") {
        if (!all(phases_to_rename_r$data == "")) {
          list_phases_to_rename <- NULL
          for (i in 1:length(phases_to_rename_r$data)) {
            list_phases_to_rename <- c(list_phases_to_rename, paste0("\"", names(phases_to_rename_r$data)[i], "\"=c(", pasteq(phases_to_rename_r$data[[i]]), ")"))
          }
          phases_to_rename <- paste0("list(", paste0(list_phases_to_rename, collapse = ","), ")")
        } else {
          phases_to_rename <- ""
        }
      } else {
        phases_to_rename <- check_phases_to_rename(input)
        if (length(grep("Error", phases_to_rename)) != 0) {
          return(phases_to_rename)
        }
      }
      if (input$solution_models_file == "load") {
        if (!exists("solution_models_file")) {
          solution_models_file <- ""
        }
      } else {
        solution_models_file <- input$solution_models_file
        use_sol_models_r$data <- input$use_sol_models
      }
      # error validation for meemum path
      if (input$meemum_path == "") {
        return("Error: No meemum path defined")
      }
      if (!file.exists(paste0(gsub("/code", "/data", getwd()), "/", input$meemum_path))) {
        return(paste0("Error: Meemum not found at ", paste0(gsub("/code", "/data", getwd()), "/", input$meemum_path)))
      }
      w_modelling_options <- paste0(
        "###############\n#\n#   Modelling Options\n#\n###############\n",
        "thermodynamic_data_file<-", pasteq(input$thermodynamic_data_file), "\n",
        "solution_models_file<-", pasteq(solution_models_file), "\n",
        "meemum_path<-", pasteq(input$meemum_path), "\n",
        "perplex_option_file<-", pasteq(input$perplex_option_file), "\n",
        "use_sol_models<-c(", pasteq(use_sol_models_r$data), ")\n",
        "saturated_components<-", pasteq(input$saturated_components), "\n",
        "saturated_phase_components<-", pasteq(input$saturated_phase_components), "\n",
        "independent_potential_fugacity_activity<-", pasteq(input$independent_potential_fugacity_activity), "\n",
        "exclude_phases<-c(", pasteq(break_on_comma(input$exclude_phases)), ")\n",
        "calculate_activities<-", input$calculate_activities, "\n",
        "phases_to_rename<-c(", bl(phases_to_rename), ")\n",
        "G_pure_phases<-", pasteq(input$G_pure_phases), "\n",
        "print_meem<-", input$print_meem, "\n",
        "export_meemum_output<-", input$export_meemum_output, "\n",
        "end_of_calc<-", pasteq(input$end_of_calc), "\n"
      )
      w_output_options <- paste0(
        "###############\n#\n#   Output Options\n#\n###############\n",
        "phase_aliases<-", pasteq(input$phase_aliases), "\n",
        "PAM_compilation<-", pasteq(input$PAM_compilation), "\n",
        "compile_PAM<-", pasteq(input$compile_PAM), "\n"
      )
      # Compile all tabs into a page
      thepage <- c(
        w_file,
        w_size,
        w_pt,
        w_bulk_composition,
        w_phase_addition,
        w_phase_extraction,
        w_component_packet,
        w_modelling_options,
        w_output_options
      )
      # If directory doesnt exist, create it
      if (!dir.exists(paste0(projects_directory, "/", working_file))) {
        dir.create(paste0(projects_directory, "/", working_file))
      }
      if (!dir.exists(paste0(projects_directory, "/", working_file, "/Inputs/"))) {
        dir.create(paste0(projects_directory, "/", working_file, "/Inputs/"))
      }
      if (!dir.exists(paste0(projects_directory, "/", working_file, "/Outputs/"))) {
        dir.create(paste0(projects_directory, "/", working_file, "/Outputs/"))
      }
      # Grab additional parameters if file already exists
      add_text <- NULL
      if (file.exists(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))) {
        scanned <- scan(
          file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"),
          what = "character",
          sep = "\n",
          quiet = TRUE
        )
        break_line <- which(scanned == "#   Additional Parameters")
        if (length(break_line) == 1) {
          add_text <- scanned[break_line:length(scanned)]
        }
      }
      # Save .txt file
      write(
        c(thepage, add_text),
        file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt")
      )
      # Save workspace
      working_file <<- input$working_file
      save.image(file = paste0(projects_directory, "/", working_file, "/", working_file, ".RData"))
      # Return success message
      return(paste0("File saved to ", projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    } else { # error handing not passed
      # Return error message
      return(reactive_message$data)
    }
  })
  # Load function
  on_load <- reactive({
    projects_directory <- input$projects_directory
    working_file <- input$working_file
    # error handling
    reactive_message$data <- error_handling(working_file, projects_directory)
    if (reactive_message$data == "error handling passed") {
      # Remove current workspace but keeping marker: from clear button
      # mod-tag: should we really be clearing workspace? can't we just overwrite values for load?
      rm(list = ls())
      # fix-tag - have to remove major elements here to allow load, must fix this
      rm("major_elements", inherits = TRUE)
      rm(list = setdiff(ls(envir = .GlobalEnv), c("keep_on_load", keep_on_load)), envir = .GlobalEnv)
      # clear reactive stores
      pt_definitions_r$data <- NULL
      bulk_definitions_r$data <- NULL
      ph_add_definitions_r$data <- NULL
      ph_extr_definitions_r$data <- NULL
      cp_components_r$data <- NULL
      cp_phases_r$data <- list()
      # Load .txt unless from clear button
      if (from_clear_button) {
        projects_directory <- input$projects_directory
        working_file <- input$working_file
      } else {
        source(paste0(input$projects_directory, "/", input$working_file, "/Inputs/", input$working_file, ".txt"))
      }
      # Values to load
      load_variables <- c(
        "x_n" = "inp",
        "y_n" = "inp",
        "n_pt_def" = "reactive",
        "n_comp_trans" = "reactive",
        "n_phases_to_rename" = "reactive",
        "bulk_def_file" = "checkbox",
        "set_oxygen_fugacity" = "checkbox",
        "calculate_traces" = "checkbox",
        "apply_trace_correction" = "select",
        "Xmz" = "inp",
        "D_ApMelt_LREE" = "inp",
        "apatite_saturation_Ap" = "select",
        "apatite_saturation_ApMnz" = "select",
        "major_elements" = "reactive_select",
        "trace_elements" = "reactive_select",
        "kd_file" = "inp", "n_bulk_def" = "reactive",
        "bulk_file" = "inp", "ph_add" = "checkbox",
        "n_ph_add_def" = "reactive",
        "ph_extr" = "checkbox",
        "reequilibrate_steps" = "checkbox",
        "n_ph_extr_def" = "reactive",
        "thermodynamic_data_file" = "inp",
        "saturated_components" = "inp",
        "saturated_phase_components" = "inp",
        "independent_potential_fugacity_activity" = "inp",
        "calculate_activities" = "checkbox",
        "G_pure_phases" = "inp",
        "exclude_phases" = "inp",
        "component_packet" = "checkbox",
        "cp_components" = "reactive",
        "print_meem" = "checkbox",
        "export_meemum_output" = "checkbox",
        "end_of_calc" = "inp",
        "solution_models_file" = "reactive",
        "perplex_option_file" = "inp",
        "meemum_path" = "inp",
        "phase_aliases" = "inp",
        "PAM_compilation" = "inp",
        "compile_PAM" = "checkbox"
      )
      # load reactive stores
      if (exists("pt_definitions")) {
        pt_definitions_r$data <- pt_definitions
      }
      if (exists("bulk_definitions")) {
        bulk_definitions_r$data <- bulk_definitions
      }
      if (exists("Xmz")) {
        if (Xmz == "") {
          Xmz <- 0.83
        }
      }
      if (exists("comp_transformations")) {
        comp_transformations_r$data <- comp_transformations
      }
      if (exists("ph_add_definitions")) {
        ph_add_definitions_r$data <- ph_add_definitions
      }
      if (exists("ph_extr_definitions")) {
        ph_extr_definitions_r$data <- ph_extr_definitions
      }
      if (exists("cp_components")) {
        cp_components_r$data <- cp_components
        cp_phases_r$data <- list()
        for (k in 1:length(cp_components)) {
          cp_phases_r$data[[k]] <- eval(parse(text = paste0("cp_phases_", names(cp_components[k]))))
          # stores in a list without naming elements in the list
        }
      }
      if (exists("solution_models_file")) {
        solution_models_file_r$data <- solution_models_file
      }
      if (exists("use_sol_models")) {
        use_sol_models_r$data <- use_sol_models
      }
      if (exists("major_elements")) {
        major_elements_r$data <- major_elements
      }
      if (exists("trace_elements")) {
        trace_elements_r$data <- trace_elements
      }
      if (exists("all_elements")) {
        all_elements_r$data <- all_elements
      }
      if (exists("phases_to_rename")) {
        phases_to_rename_r$data <- phases_to_rename
      }
      # custom loads
      # bulk_def
      if (exists("bulk_def")) {
        if (bulk_def == "input") {
          bulk_def_file <- FALSE
        } else {
          bulk_def_file <- TRUE
        }
      } else {
        bulk_def_file <- FALSE
      }
      # get addition phases
      if (is.null(ph_add_definitions_r$data[[1]][1])) {
        add_phases <- ""
      } else {
        add_phases <- paste(setdiff(names(ph_add_definitions_r$data[[1]]), "condition"), collapse = ",")
      }
      # get extract phases
      if (is.null(ph_extr_definitions_r$data[[1]][1])) {
        extr_phases <- ""
      } else {
        extr_phases <- paste(setdiff(names(ph_extr_definitions_r$data[[1]]), "condition"), collapse = ",")
      }
      # load values
      for (i in 1:length(load_variables)) {
        if (load_variables[i] == "reactive") {
          updateTextInput(session, names(load_variables)[i], value = "load")
        }
        if (load_variables[i] == "reactive_select") {
          updateSelectizeInput(session, names(load_variables)[i], selected = "load")
        }
        if (load_variables[i] == "inp") {
          if (exists(names(load_variables)[i])) {
            if (names(load_variables)[i] == "exclude_phases") {
              updateTextInput(session, names(load_variables)[i], value = pasteq(eval(parse(text = names(load_variables)[i]))))
            } else {
              updateTextInput(session, names(load_variables)[i], value = eval(parse(text = names(load_variables)[i])))
            }
          } else {
            updateTextInput(session, names(load_variables)[i], value = "")
          }
        }
        # modtag - remove this
        if (load_variables[i] == "def_count") {
          def_name <- sub("n_", "", sub("_def", "_definitions", names(load_variables)[i]))
          if (exists(def_name)) {
            if (is.null(eval(parse(text = def_name))[[1]][1])) {
              val_out <- ""
            } else {
              if (eval(parse(text = def_name))[[1]][1] == "") {
                val_out <- ""
              } else {
                val_out <- length(eval(parse(text = def_name)))
              }
            }
          } else {
            val_out <- ""
          }
          updateTextInput(session, names(load_variables)[i], value = val_out)
        }
        if (load_variables[i] == "checkbox") {
          if (exists(names(load_variables)[i])) {
            updateCheckboxInput(session, names(load_variables)[i], value = eval(parse(text = names(load_variables)[i])))
          } else {
            # False is default for all load checkboxes
            updateCheckboxInput(session, names(load_variables)[i], value = FALSE)
          }
        }
        if (load_variables[i] == "select") {
          if (exists(names(load_variables)[i])) {
            updateSelectInput(session, names(load_variables)[i], selected = eval(parse(text = names(load_variables)[i])))
          } else {
            updateSelectInput(session, names(load_variables)[i], selected = "")
          }
        }
      }
      # Return success message
      if (from_clear_button) {
        from_clear_button <<- FALSE
        return(paste0("Inputs and workspace cleared"))
      } else {
        return(paste0("Loaded ", projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      }
    } else {
      # Return error message
      return(reactive_message$data)
    }
  })
  # Toolbar buttons
  # Load button
  observeEvent(input$load, {
    # load the file from input$working_file
    working_file <- input$working_file
    projects_directory <- input$projects_directory
    # error handling
    reactive_message$data <- error_handling(working_file, projects_directory)
    if (reactive_message$data == "error handling passed") {
      # load
      if (file.exists(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))) {
        reactive_message$data <- paste0(on_load())
      } else {
        reactive_message$data <- paste0("No input file found at ", paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      }
      # Load workspace if it exists (previous calculation results)
      if (file.exists(paste0(projects_directory, "/", working_file, "/", working_file, ".RData"))) {
        load(paste0(projects_directory, "/", working_file, "/", working_file, ".RData"), envir = .GlobalEnv)
      }
      # Refresh reactive outputs if they exist
      if (exists("crust")) {
        store_r$crust_r <- crust
      }
      if (exists("input_pt")) {
        store_r$input_pt_r <- input_pt
      }
      if (exists("input_bulk")) {
        store_r$input_bulk_r <- input_bulk
      }
      if (exists("major_elements")) {
        store_r$major_elements_r <- major_elements
      }
      if (exists("trace_elements")) {
        store_r$trace_elements_r <- trace_elements
      }
      if (exists("all_elements")) {
        store_r$all_elements_r <- all_elements
      }
    } else {
      reactive_message$data
    }
  })
  # Clear button
  observeEvent(input$clear, {
    # delete selection inputs
    from_clear_button <<- TRUE
    reactive_message$data <- paste0(on_load())
  })
  # Console button
  observeEvent(input$console, {
    cat("To regain access to the Rcrust GUI type \'c\' then press enter\n")
    browser()
  })
  check_root <- function(proj_directory) {
    if (!length(unlist(gregexpr("[/]", proj_directory))) == 2) {
      showModal(modalDialog(
        tags$p(paste0(
          "Warning Rcrust is not located in a root directory. Your project is located in ",
          proj_directory, "."
        )),
        tags$p(paste0("This should be e.g. C:/Rcrust/ or D:/Rcrust. Please relocate Rcrust to a root directory.")),
        footer = tagList(
          modalButton("OK")
        )
      ))
    }
  }
  # Run Button
  observeEvent(input$run, {
    # Grab working file and projects directory
    working_file <- input$working_file
    projects_directory <- input$projects_directory
    # check for project root directory. (Displays slightly annoying popup)
    check_root(projects_directory)
    # error handling
    reactive_message$data <- error_handling(working_file, projects_directory)
    if (reactive_message$data == "error handling passed") {
      # if error handling is passed
      # Save
      reactive_message$data <- paste0(on_save())
      # source the saved variables into the workspace
      source(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      # Run
      source("main.r")
      reactive_message$data <- paste0(
        "Calculation complete, Results saved to ",
        projects_directory, "/",
        working_file, "/",
        working_file, ".RData\n Select outputs throught the 'Outputs' tab"
      )
      # Save copy to directory
      if (FALSE) {
        copy_directory <<- "H:/Rcrust/Projects"
        from_copy <<- TRUE
        reactive_message$data <- paste0(on_save())
      }
      # End of Calculation
      if (!exists("end_of_calc")) {
        end_of_calc <- "Return to Interface"
      }
      # Email report
      if (FALSE) {
        library(gmailr)
        mime() %>%
          to("mjmayne@outlook.com") %>%
          from("matt.mayne.1992@gmail.com") %>%
          html_body("Body text") -> html_msg
        html_msg %>%
          subject(paste0(working_file, " completed successfully after ", run_time, "with finish option ", end_of_calc)) %>%
          # attach_file(paste0(projects_directory,"/",working_file,"/",working_file,".RData")) -> file_attachment
          # attach_file(paste0(projects_directory,"/Yak_PT fh any_all/Yak_PT fh any_all.RData")) -> file_attachment
          send_message()
      }
      if (end_of_calc == "Logout") {
        system("shutdown -l")
      }
      if (end_of_calc == "Shutdown") {
        system("shutdown -s")
      }
      # Refresh reactive outputs if they exist
      if (exists("crust")) {
        store_r$crust_r <- crust
      }
      if (exists("input_pt")) {
        store_r$input_pt_r <- input_pt
      }
      if (exists("input_bulk")) {
        store_r$input_bulk_r <- input_bulk
      }
      if (exists("major_elements")) {
        store_r$major_elements_r <- major_elements
      }
      if (exists("trace_elements")) {
        store_r$trace_elements_r <- trace_elements
      }
    }
  })
  # Save button
  observeEvent(input$save, {
    reactive_message$data <- paste0(on_save())
  })
  # Send error/success messages to GUI
  output$print_message <- renderText({
    if (is.null(reactive_message$data)) {
      return()
    }
    paste0(reactive_message$data)
  })
  # Dyanmically use number of PT definitions to create the correct number of From,To,P,T inputs
  output$pt <- renderUI({
    if (input$n_pt_def == "load") {
      if (all(pt_definitions_r$data == "") | is.null(pt_definitions_r$data)) {
        def_num <- ""
      } else {
        def_num <- length(pt_definitions_r$data)
      }
      updateTextInput(session, "n_pt_def", value = def_num)
    } else {
      def_num <- input$n_pt_def
    }
    if (!(is.null(def_num) | def_num == "" | def_num == 0)) {
      validate(if (is.na(suppressWarnings(as.numeric(def_num)))) {
        "Error: Number of PT definitions must be numeric"
      } else {
        NULL
      })
      validate(if (as.numeric(def_num) < 1) {
        "Error: Number of PT definitions must be greater than 0"
      } else {
        NULL
      })
      validate(if (!as.numeric(def_num) %% 1 == 0) {
        "Error: Number of PT definitions must be a whole number"
      } else {
        NULL
      })
      fixedRow(
        lapply(1:(as.numeric(def_num) * 4), function(i) {
          a <- 1
          ii <- i
          while (ii > 4) {
            ii <- ii - 4
            a <- a + 1
          }
          if (ii == 1) {
            chk_From <- try(unlist(strsplit(names(pt_definitions_r$data)[a], split = "_"))[1], silent = TRUE)
            if (class(chk_From) == "try-error") {
              chk_From <- NULL
            }
            column(2, textInput(paste0("pt_from_", a), "From", value = chk_From))
          } else {
            if (ii == 2) {
              chk_To <- try(unlist(strsplit(names(pt_definitions_r$data)[a], split = "_"))[2], silent = TRUE)
              if (class(chk_To) == "try-error") {
                chk_To <- NULL
              }
              column(2, textInput(paste0("pt_to_", a), "To", value = chk_To))
            } else {
              if (ii == 3) {
                chk_pressure <- try(pt_definitions_r$data[[a]][1], silent = TRUE)
                if (class(chk_pressure) == "try-error") {
                  chk_pressure <- NULL
                }
                column(4, textInput(paste0("pressure_", a), "Pressure (kbar)", value = chk_pressure))
              } else {
                if (ii == 4) {
                  chk_temperature <- try(pt_definitions_r$data[[a]][2], silent = TRUE)
                  if (class(chk_temperature) == "try-error") {
                    chk_temperature <- NULL
                  }
                  column(4, textInput(paste0("temperature_", a), "Temperature (C)", value = chk_temperature))
                }
              }
            }
          }
        })
      )
    }
  })
  # Import P-T definitions button
  observeEvent(input$import_pt, {
    # bugtag - try rather updating text to have only a local load #updateTextInput(session,input$n_ph_add_def,value="load")
    reactive_message$data <- paste0(on_save())
    # Import definitions from file (csv of definitions or txt as Rcrust input file)
    # read Rcrust Input file
    if (input$file_pt["type"] == "text/plain") {
      data_in <- readLines(paste(input$file_pt[4]))
      input_file <- readLines(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      if (!length(grep("pt_definitions<-", data_in)) == 0) {
        if (length(grep("pt_definitions<-", input_file)) == 0) {
          input_file <- c(input_file, grep("pt_definitions<-", data_in))
        } else {
          input_file[grep("pt_definitions<-", input_file)] <- data_in[grep("pt_definitions<-", data_in)]
        }
      }
      write(input_file, file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    } else {
      # read csv of definitions
      data_in <- as.matrix(read.csv(paste(input$file_pt[4])))
      lines_all <- NULL
      for (line_no in 1:nrow(data_in)) {
        line_i <- paste(
          "\"{",
          data_in[line_no, 1], ";",
          data_in[line_no, 2], "}_{",
          data_in[line_no, 3], ";",
          data_in[line_no, 4], "}\"=c(\"",
          data_in[line_no, 5], "\",\"",
          data_in[line_no, 6], "\")",
          sep = ""
        )
        if (line_no == 1) {
          lines_all <- paste("pt_definitions<-list(", line_i, sep = "")
        } else {
          lines_all <- paste(lines_all, line_i, sep = ",")
        }
        if (line_no == nrow(data_in)) {
          lines_all <- paste(lines_all, ")", sep = "")
        }
      }
      input_file <- readLines(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      if (length(grep("pt_definitions<-", input_file)) == 0) {
        input_file <- c(input_file, lines_all)
      } else {
        input_file[grep("pt_definitions<-", input_file)] <- lines_all
      }
      write(input_file, file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    }
    reactive_message$data <- paste0(on_load())
  })
  # Dynamically use number of component transformations to create the correct number of inputs
  output$trans <- renderUI({
    if (input$n_comp_trans == "load") {
      if (all(comp_transformations_r$data == "") | is.null(comp_transformations_r$data)) {
        def_num <- ""
      } else {
        def_num <- length(comp_transformations_r$data)
      }
      updateTextInput(session, "n_comp_trans", value = def_num)
    } else {
      def_num <- input$n_comp_trans
    }
    if (!(is.null(def_num) | def_num == "" | def_num == 0)) {
      validate(if (is.na(suppressWarnings(as.numeric(def_num)))) {
        "Error: Number of Component transformations must be numeric"
      } else {
        NULL
      })
      validate(if (as.numeric(def_num) < 1) {
        "Error: Number of Component transformations must be greater than 0"
      } else {
        NULL
      })
      validate(if (!as.numeric(def_num) %% 1 == 0) {
        "Error: Number of Components transformations must be a whole number"
      } else {
        NULL
      })
      fixedRow(
        lapply(1:(as.numeric(def_num) * 3), function(i) {
          a <- 1
          ii <- i
          while (ii > 3) {
            ii <- ii - 3
            a <- a + 1
          }
          if (ii == 1) {
            chk_old <- try(unlist(strsplit(names(comp_transformations_r$data)[a], split = "_"))[1], silent = TRUE)
            if (class(chk_old) == "try-error") {
              chk_old <- NULL
            }
            column(2, textInput(paste0("old_comp_", a), "Replace component", value = chk_old))
          } else {
            # mod-tag: better to have a select input but for some reason it doesnt work well with column()
            # column(2,selectizeInput(paste0('old_comp_',a),'Replace',c("",available_components_r$data),selected=chk_old))}else{
            if (ii == 2) {
              chk_new <- try(unlist(strsplit(names(comp_transformations_r$data)[a], split = "_"))[2], silent = TRUE)
              if (class(chk_new) == "try-error") {
                chk_new <- NULL
              }
              column(2, textInput(paste0("new_comp_", a), "New component", value = chk_new))
            } else {
              if (ii == 3) {
                comp_label <- paste(available_components_r$current, collapse = ",")
                chk_comp <- try(comp_transformations_r$data[[a]], silent = TRUE)
                if (class(chk_comp) == "try-error") {
                  chk_comp <- NULL
                }
                column(8, textInput(paste0("comp_", a), comp_label, value = paste(chk_comp, collapse = ",")))
              } else {
              }
            }
          }
        })
      )
    }
  })
  # Component transformation
  observe({
    available_components <- NULL
    suppressWarnings(
      qq_try <- try(
        scan(
          file = gsub(
            "Rcrust/code",
            paste0("Rcrust/data/", input$thermodynamic_data_file),
            getwd()
          ),
          what = "character",
          sep = "\n",
          quiet = TRUE
        ),
        silent = TRUE
      )
    )
    if (class(qq_try) != "try-error") {
      qq <- qq_try
      start_foo <- grep("begin_components", qq)
      end_foo <- grep("end_components", qq)
      qq <- qq[(start_foo + 1):(end_foo - 1)]
      qq <- strsplit(qq, " ")
      for (i in 1:length(qq)) {
        available_components <- c(available_components, qq[[i]][1])
      }
    }
    available_components_r$current <- available_components
  })
  # Component transformation
  observe({
    # use current component transformation inputs to update current components_r and selection in major elements
    current_components <- available_components_r$current
    chk_trans <- FALSE
    if (!is.na(suppressWarnings(as.numeric(input$n_comp_trans)))) {
      if (!as.numeric(input$n_comp_trans) < 1) {
        if (as.numeric(input$n_comp_trans) %% 1 == 0) {
          chk_trans <- TRUE
        }
      }
    }
    if (chk_trans) {
      for (i in 1:input$n_comp_trans) {
        pos <- which(current_components == eval(parse(text = paste0("input$old_comp_", i))))
        current_components[pos] <- eval(parse(text = paste0("input$new_comp_", i)))
      }
    }
    current_components_r$data <- current_components
  })
  # Observe component packet
  # if component is not available_components_r$current then adds to current_components_r$data
  observe({
    current_components <- available_components_r$current
    chk_packet_component <- FALSE
    if (!isTRUE(input$cp_components == "load")) {
      if (!isTRUE(input$cp_components == "")) {
        packet_components <- c(break_on_comma(input$cp_components))
        if (any(is.na(match(packet_components, available_components_r$current)))) {
          chk_packet_component <- TRUE
        }
      }
    }
    if (chk_packet_component) {
      for (i in 1:length(packet_components)) {
        if (is.na(match(packet_components[i], available_components_r$current))) {
          current_components <- c(current_components, packet_components[i])
        }
      }
    }
    current_components_r$data <- current_components
  }) # end of Observe component packet
  observe({
    # use current component transformation inputs to update current components_r and selection in major elements
    current_components <- available_components_r$current
    chk_trans <- FALSE
    if (!is.na(suppressWarnings(as.numeric(input$n_comp_trans)))) {
      if (!as.numeric(input$n_comp_trans) < 1) {
        if (as.numeric(input$n_comp_trans) %% 1 == 0) {
          chk_trans <- TRUE
        }
      }
    }
    if (chk_trans) {
      for (i in 1:input$n_comp_trans) {
        pos <- which(current_components == eval(parse(text = paste0("input$old_comp_", i))))
        current_components[pos] <- eval(parse(text = paste0("input$new_comp_", i)))
      }
    }
    current_components_r$data <- current_components
  })
  output$maj <- renderUI({
    # fix-tag: does not work if have multiple component transformations, scoping means we come here after first transformation is complete because we alter "current components"
    if (!input$n_comp_trans == "load") {
      if (any(input$major_elements == "load")) {
        selectizeInput(
          "major_elements",
          "Major elements",
          c(
            major_elements_r$data,
            setdiff(current_components_r$data, major_elements_r$data),
            "load"
          ),
          selected = major_elements_r$data,
          multiple = TRUE
        )
      } else {
        selectizeInput(
          "major_elements",
          "Major elements",
          c(
            input$major_elements,
            setdiff(current_components_r$data, input$major_elements),
            "load"
          ),
          selected = input$major_elements,
          multiple = TRUE
        )
      }
    }
  })
  output$traces <- renderUI({
    traces_file <- ""
    # If valid kd file is input, trace elements from columns of kd file become available.
    if (!input$kd_file == "") {
      if (!file.exists(paste0(gsub("/code", "/data", getwd()), "/", input$kd_file))) {
        traces_message$data <- paste0("Provide valid Kd file in order to select trace elements")
        show_traces <- TRUE
      } else {
        traces_message$data <- paste0("")
        traces_file <- colnames(read.table(paste0(gsub("/code", "/data", getwd()), "/", input$kd_file), sep = "\t"))
      }
    }
    # When Apatite saturation is selected, the trace elements input field has P2O5 as option.
    if (input$apply_trace_correction == "Apatite saturation" || input$apply_trace_correction == "Apatite & Monazite Saturation") {
      if (is.na(match("P2O5", traces_file))) {
        traces_file <- c(traces_file, "P2O5")
      }
      if (is.na(match("P2O5", input$trace_elements))) {
        updateSelectizeInput(session, "trace_elements", selected = traces_file)
      }
    } else {
      if (!is.na(match("P2O5", traces_file))) {
        traces_file <- traces_file[-match("P2O5", trace_elements)]
      }
      if (!is.na(match("P2O5", trace_elements_r$data))) {
        trace_elements_r$data <- trace_elements_r$data[-match("P2O5", trace_elements)]
      }
      if (!is.na(match("P2O5", input$trace_elements))) {
        updateSelectizeInput(session, "trace_elements", selected = "load")
      }
    }
    if (any(input$trace_elements == "load")) {
      selectizeInput("trace_elements", "Trace elements",
        c(
          trace_elements_r$data,
          setdiff(traces_file, trace_elements_r$data),
          "load"
        ),
        selected = trace_elements_r$data,
        multiple = TRUE
      )
    } else {
      selectizeInput("trace_elements", "Trace elements",
        c(
          input$trace_elements,
          setdiff(traces_file, input$trace_elements),
          "load"
        ),
        selected = input$trace_elements,
        multiple = TRUE
      )
    }
  })
  # Warning message for valide kd file
  output$kd_message <- renderText({
    if (is.null(traces_message$data)) {
      return()
    }
    paste0(traces_message$data)
  })
  # Dynamically use number of bulk definitions to create the correct number of From,To,bulk inputs
  output$bulk <- renderUI({
    if (input$n_bulk_def == "load") {
      if (all(bulk_definitions_r$data == "") | is.null(bulk_definitions_r$data)) {
        def_num <- ""
      } else {
        def_num <- length(bulk_definitions_r$data)
      }
      updateTextInput(session, "n_bulk_def", value = def_num)
    } else {
      def_num <- input$n_bulk_def
    }
    if (!(is.null(def_num) | def_num == "" | def_num == 0)) {
      validate(if (is.na(suppressWarnings(as.numeric(def_num)))) {
        "Error: Number of Bulk definitions must be numeric"
      } else {
        NULL
      })
      validate(if (as.numeric(def_num) < 1) {
        "Error: Number of Bulk definitions must be greater than 0"
      } else {
        NULL
      })
      validate(if (!as.numeric(def_num) %% 1 == 0) {
        "Error: Number of Bulk definitions must be a whole number"
      } else {
        NULL
      })
      fixedRow(
        lapply(1:(as.numeric(def_num) * 3), function(i) {
          a <- 1
          ii <- i
          while (ii > 3) {
            ii <- ii - 3
            a <- a + 1
          }
          if (ii == 1) {
            chk_From <- try(unlist(strsplit(names(bulk_definitions_r$data)[a], split = "_"))[1], silent = TRUE)
            if (class(chk_From) == "try-error") {
              chk_From <- NULL
            }
            column(2, textInput(paste0("bulk_from_", a), "From", value = chk_From))
          } else {
            if (ii == 2) {
              chk_To <- try(unlist(strsplit(names(bulk_definitions_r$data)[a], split = "_"))[2], silent = TRUE)
              if (class(chk_To) == "try-error") {
                chk_To <- NULL
              }
              column(2, textInput(paste0("bulk_to_", a), "To", value = chk_To))
            } else {
              if (ii == 3) {
                if (is.null(input$major_elements)) {
                  bulk_label <- "Please choose compositional elements above"
                } else {
                  if (input$set_oxygen_fugacity) {
                    bulk_label <- paste(c(input$major_elements, input$trace_elements, "log10(fugacity)", "mass"), collapse = ",")
                  } else {
                    bulk_label <- paste(c(input$major_elements, input$trace_elements, "mass"), collapse = ",")
                  }
                }
                chk_bulk <- try(bulk_definitions_r$data[[a]], silent = TRUE)
                if (class(chk_bulk) == "try-error") {
                  chk_bulk <- NULL
                }
                column(8, textInput(paste0("bulk_", a), bulk_label, value = paste(chk_bulk, collapse = ",")))
              } else {
              }
            }
          }
        })
      )
    }
  })
  # Import bulk definitions button
  observeEvent(input$import_bulk, {
    # mod-tag: try rather updating text to have only a local load #updateTextInput(session,input$n_ph_add_def,value="load")
    reactive_message$data <- paste0(on_save())
    # Import definitions from file (csv of definitions or txt as Rcrust input file)
    # read Rcrust Input file
    if (input$file_bulk["type"] == "text/plain") {
      data_in <- readLines(paste(input$file_bulk[4]))
      input_file <- readLines(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      if (!length(grep("comp_transformations<-", data_in)) == 0) {
        if (length(grep("comp_transformations<-", input_file)) == 0) {
          input_file <- c(input_file, grep("comp_transformations<-", data_in))
        } else {
          input_file[grep("comp_transformations<-", input_file)] <- data_in[grep("comp_transformations<-", data_in)]
        }
      }
      if (!length(grep("major_elements<-", data_in)) == 0) {
        if (length(grep("major_elements<-", input_file)) == 0) {
          input_file <- c(input_file, grep("major_elements<-", data_in))
        } else {
          input_file[grep("major_elements<-", input_file)] <- data_in[grep("major_elements<-", data_in)]
        }
      }
      if (!length(grep("bulk_definitions<-", data_in)) == 0) {
        if (length(grep("bulk_definitions<-", input_file)) == 0) {
          input_file <- c(input_file, grep("bulk_definitions<-", data_in))
        } else {
          input_file[grep("bulk_definitions<-", input_file)] <- data_in[grep("bulk_definitions<-", data_in)]
        }
      }
      write(input_file, file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    } else {
      # read csv of definitions
      data_in <- as.matrix(read.csv(paste(input$file_bulk[4])))
      lines_all <- NULL
      for (line_no in 1:nrow(data_in)) {
        line_i <- paste(
          "\"{",
          data_in[line_no, 1], ";",
          data_in[line_no, 2], "}_{",
          data_in[line_no, 3], ";",
          data_in[line_no, 4], "}\"=c(\"",
          paste(data_in[line_no, c(-1, -2, -3, -4)], collapse = "\",\""),
          "\")",
          sep = ""
        )
        if (line_no == 1) {
          lines_all <- paste("bulk_definitions<-c(list(", line_i, sep = "")
        } else {
          lines_all <- paste(lines_all, line_i, sep = ",")
        }
        if (line_no == nrow(data_in)) {
          lines_all <- paste(lines_all, "))", sep = "")
        }
      }
      input_file <- readLines(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      if (length(grep("bulk_definitions<-", input_file)) == 0) {
        input_file <- c(input_file, lines_all)
      } else {
        input_file[grep("bulk_definitions<-", input_file)] <- lines_all
      }
      write(input_file, file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    }
    reactive_message$data <- paste0(on_load())
  })
  # Dyanmically use number of Phase Addition definitions to create the correct number of From,To,P,T inputs
  output$ph_add <- renderUI({
    if (input$n_ph_add_def == "load") {
      if (all(ph_add_definitions_r$data == "") | is.null(ph_add_definitions_r$data)) {
        def_num <- ""
      } else {
        def_num <- length(ph_add_definitions_r$data)
      }
      updateTextInput(session, "n_ph_add_def", value = def_num)
    } else {
      def_num <- input$n_ph_add_def
    }
    if (!(is.null(def_num) | def_num == "" | def_num == 0)) {
      validate(if (is.na(suppressWarnings(as.numeric(def_num)))) {
        "Error: Number of Phase Addition definitions must be numeric"
      } else {
        NULL
      })
      validate(if (as.numeric(def_num) < 1) {
        "Error: Number of Phase Addition definitions must be greater than 0"
      } else {
        NULL
      })
      validate(if (!as.numeric(def_num) %% 1 == 0) {
        "Error: Number of Phase Addition definitions must be a whole number"
      } else {
        NULL
      })
      fixedRow(
        # Mod-tag: could use some commenting here.
        lapply(1:(as.numeric(def_num) * 5), function(i) {
          a <- 1
          ii <- i
          while (ii > 5) {
            ii <- ii - 5
            a <- a + 1
          }
          if (ii == 1) {
            chk_From <- try(unlist(strsplit(names(ph_add_definitions_r$data)[a], split = "_"))[1], silent = TRUE)
            if (class(chk_From) == "try-error") {
              chk_From <- NULL
            }
            column(2, textInput(paste0("ph_add_from_", a), "From", value = chk_From))
          } else {
            if (ii == 2) {
              chk_To <- try(unlist(strsplit(names(ph_add_definitions_r$data)[a], split = "_"))[2], silent = TRUE)
              if (class(chk_To) == "try-error") {
                chk_To <- NULL
              }
              column(2, textInput(paste0("ph_add_to_", a), "To", value = chk_To))
            } else {
              if (ii == 3) {
                chk_Condition <- try(ph_add_definitions_r$data[[a]][["condition"]], silent = TRUE)
                if (class(chk_Condition) == "try-error") {
                  chk_Condition <- NULL
                }
                column(3, textInput(paste0("ph_add_con_", a), "Condition", value = chk_Condition))
              } else {
                if (ii == 4) {
                  chk_Phases <- try(paste(names(ph_add_definitions_r$data[[a]][-1]), collapse = ","), silent = TRUE)
                  if (class(chk_Phases) == "try-error") {
                    chk_Phases <- NULL
                  }
                  column(4, textInput(paste0("ph_add_phs_", a), "Phases", value = chk_Phases))
                } else {
                  if (ii == 5) {
                    column(12, uiOutput(paste0("ph_add_", a)))
                  }
                }
              }
            }
          }
        })
      )
    }
  })
  # auto create phase addition inputs given number of phases
  observe({
    if (all(
      !input$n_ph_add_def == "load",
      !is.null(input$n_ph_add_def),
      !input$n_ph_add_def == "",
      !input$n_ph_add_def == 0
    )) {
      lapply(1:input$n_ph_add_def, function(i) {
        eval(parse(text = paste0("output$ph_add_", i, "<-renderUI({
          a<-i
          chk_Phases<-try(eval(parse(text=paste0(\'input$ph_add_phs_\',a))),silent=TRUE)
					if(!is.null(chk_Phases)){
            if(!chk_Phases==\"\"){
              phases<-unlist(strsplit(chk_Phases,split=\",\"))
              fixedRow(
                lapply(1:(length(phases)), function(j) {
                chk_in_ph<-try(ph_add_definitions_r$data[[a]][phases[j]],silent=TRUE)
                if(class(chk_in_ph)==\"try-error\"){chk_in_ph<-NULL}
                column(3,textInput(paste0(\'ph_add_phs_\',a,\'_\',phases[j]),phases[j],value=chk_in_ph))
              }))
						}
          }
        })")))
      })
    }
  })
  # Import addition definitions button
  observeEvent(input$import_ph_add, {
    # mod-tag: try rather updating text to have only a local load #updateTextInput(session,input$n_ph_add_def,value="load")
    reactive_message$data <- paste0(on_save())
    # Import definitions from file (csv of definitions or txt as Rcrust input file)
    # read Rcrust Input file
    if (input$file_ph_add["type"] == "text/plain") {
      data_in <- readLines(paste(input$file_ph_add[4]))
      input_file <- readLines(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      if (!length(grep("ph_add<-", data_in)) == 0) {
        if (length(grep("ph_add<-", input_file)) == 0) {
          input_file <- c(input_file, grep("ph_add<-", data_in))
        } else {
          input_file[grep("ph_add<-", input_file)] <- data_in[grep("ph_add<-", data_in)]
        }
      }
      if (!length(grep("ph_add_definitions<-", data_in)) == 0) {
        if (length(grep("ph_add_definitions<-", input_file)) == 0) {
          input_file <- c(input_file, grep("ph_add_definitions<-", data_in))
        } else {
          input_file[grep("ph_add_definitions<-", input_file)] <- data_in[grep("ph_add_definitions<-", data_in)]
        }
      }
      write(input_file, file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    } else {
      # read csv of definitions
      data_in <- as.matrix(read.csv(paste(input$file_ph_add[4])))
      lines_all <- NULL
      for (line_no in 1:nrow(data_in)) {
        line_i <- paste(
          "\"{",
          data_in[line_no, 1], ";",
          data_in[line_no, 2], "}_{",
          data_in[line_no, 3], ";",
          data_in[line_no, 4], "}\"=c(condition=\"",
          data_in[line_no, 5], "\",",
          data_in[line_no, 6], ")",
          sep = ""
        )
        if (line_no == 1) {
          lines_all <- paste("ph_add_definitions<-c(list(", line_i, sep = "")
        } else {
          lines_all <- paste(lines_all, line_i, sep = ",")
        }
        if (line_no == nrow(data_in)) {
          lines_all <- paste(lines_all, "))", sep = "")
        }
      }
      input_file <- readLines(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      if (length(grep("ph_add<-", input_file)) == 0) {
        input_file <- c(input_file, "ph_add<-TRUE")
      } else {
        input_file[grep("ph_add<-", input_file)] <- "ph_add<-TRUE"
      }
      if (length(grep("ph_add_definitions<-", input_file)) == 0) {
        input_file <- c(input_file, lines_all)
      } else {
        input_file[grep("ph_add_definitions<-", input_file)] <- lines_all
      }
      write(input_file, file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    }
    reactive_message$data <- paste0(on_load())
    updateTextInput(session, input$file_ph_add, value = "")
  })
  # Import extraction definitions button
  observeEvent(input$import_ph_extr, {
    # mod-tag: try rather updating text to have only a local load #updateTextInput(session,input$n_ph_extr_def,value="load")
    reactive_message$data <- paste0(on_save())
    # Import definitions from file (csv of definitions or txt as Rcrust input file)
    # read Rcrust Input file
    if (input$file_ph_extr["type"] == "text/plain") {
      data_in <- readLines(paste(input$file_ph_extr[4]))
      input_file <- readLines(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      if (!length(grep("ph_extr<-", data_in)) == 0) {
        if (length(grep("ph_extr<-", input_file)) == 0) {
          input_file <- c(input_file, grep("ph_extr<-", data_in))
        } else {
          input_file[grep("ph_extr<-", input_file)] <- data_in[grep("ph_extr<-", data_in)]
        }
      }
      if (!length(grep("reequilibrate_steps<-", data_in)) == 0) {
        if (length(grep("reequilibrate_steps<-", input_file)) == 0) {
          input_file <- c(input_file, grep("reequilibrate_steps<-", data_in))
        } else {
          input_file[grep("reequilibrate_steps<-", input_file)] <- data_in[grep("reequilibrate_steps<-", data_in)]
        }
      }
      if (!length(grep("ph_extr_definitions<-", data_in)) == 0) {
        if (length(grep("ph_extr_definitions<-", input_file)) == 0) {
          input_file <- c(input_file, grep("ph_extr_definitions<-", data_in))
        } else {
          input_file[grep("ph_extr_definitions<-", input_file)] <- data_in[grep("ph_extr_definitions<-", data_in)]
        }
      }
      write(input_file, file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    } else {
      # read csv of definitions
      data_in <- as.matrix(read.csv(paste(input$file_ph_extr[4])))
      lines_all <- NULL
      for (line_no in 1:nrow(data_in)) {
        line_i <- paste(
          "\"{",
          data_in[line_no, 1], ";",
          data_in[line_no, 2], "}_{",
          data_in[line_no, 3], ";",
          data_in[line_no, 4], "}\"=c(condition=\"",
          data_in[line_no, 5], "\",",
          data_in[line_no, 6], ")",
          sep = ""
        )
        if (line_no == 1) {
          lines_all <- paste("ph_extr_definitions<-c(list(", line_i, sep = "")
        } else {
          lines_all <- paste(lines_all, line_i, sep = ",")
        }
        if (line_no == nrow(data_in)) {
          lines_all <- paste(lines_all, "))", sep = "")
        }
      }
      input_file <- readLines(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
      input_file[grep("ph_extr_definitions", input_file)] <- lines_all
      write(input_file, file = paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
    }
    reactive_message$data <- paste0(on_load())
    updateTextInput(session, input$file_ph_extr, value = "")
  })
  # Dyanmically use number of Phase Extraction definitions to create the correct number of From,To,P,T inputs
  output$ph_extr <- renderUI({
    # Mod-tag: could use some commenting here
    if (input$n_ph_extr_def == "load") {
      if (all(ph_extr_definitions_r$data == "") | is.null(ph_extr_definitions_r$data)) {
        def_num <- ""
      } else {
        def_num <- length(ph_extr_definitions_r$data)
      }
      updateTextInput(session, "n_ph_extr_def", value = def_num)
    } else {
      def_num <- input$n_ph_extr_def
    }
    if (!(is.null(def_num) | def_num == "" | def_num == 0)) {
      validate(if (is.na(suppressWarnings(as.numeric(def_num)))) {
        "Error: Number of Phase Extraction definitions must be numeric"
      } else {
        NULL
      })
      validate(if (as.numeric(def_num) < 1) {
        "Error: Number of Phase Extraction definitions must be greater than 0"
      } else {
        NULL
      })
      validate(if (!as.numeric(def_num) %% 1 == 0) {
        "Error: Number of Phase Extraction definitions must be a whole number"
      } else {
        NULL
      })
      fixedRow(
        lapply(1:(as.numeric(def_num) * 5), function(i) {
          a <- 1
          ii <- i
          while (ii > 5) {
            ii <- ii - 5
            a <- a + 1
          }
          if (ii == 1) {
            chk_From <- try(unlist(strsplit(names(ph_extr_definitions_r$data)[a], split = "_"))[1], silent = TRUE)
            if (class(chk_From) == "try-error") {
              chk_From <- NULL
            }
            column(2, textInput(paste0("ph_extr_from_", a), "From", value = chk_From))
          } else {
            if (ii == 2) {
              chk_To <- try(unlist(strsplit(names(ph_extr_definitions_r$data)[a], split = "_"))[2], silent = TRUE)
              if (class(chk_To) == "try-error") {
                chk_To <- NULL
              }
              column(2, textInput(paste0("ph_extr_to_", a), "To", value = chk_To))
            } else {
              if (ii == 3) {
                chk_Condition <- try(ph_extr_definitions_r$data[[a]][["condition"]], silent = TRUE)
                if (class(chk_Condition) == "try-error") {
                  chk_Condition <- NULL
                }
                column(3, textInput(paste0("ph_extr_con_", a), "Condition", value = chk_Condition))
              } else {
                if (ii == 4) {
                  phase_names <- try(names(ph_extr_definitions_r$data[[a]][-1]), silent = TRUE)
                  chk_Phases <- NULL
                  if (class(phase_names) != "try-error") {
                    for (i in 1:length(phase_names)) {
                      if (i >= 2) {
                        chk_Phases <- paste0(chk_Phases, ",")
                      }
                      # place quotes around any phase with a comma in its name
                      if (length(grep(",", phase_names[i])) > 0) {
                        chk_Phases <- paste0(chk_Phases, pasteq(phase_names[i]))
                      } else {
                        chk_Phases <- paste0(chk_Phases, phase_names[i])
                      }
                    }
                  }
                  column(4, textInput(paste0("ph_extr_phs_", a), "Phases", value = chk_Phases))
                } else {
                  if (ii == 5) {
                    column(12, uiOutput(paste0("ph_extr_", a)))
                  }
                }
              }
            }
          }
        })
      )
    }
  })
  # auto create phase extraction inputs given number of phases
  # remove phases<-unlist(strsplit(chk_Phases,split=\",\"))
  observe({
    if (all(
      !input$n_ph_extr_def == "load",
      !is.null(input$n_ph_extr_def),
      !input$n_ph_extr_def == "",
      !input$n_ph_extr_def == 0
    )) {
      lapply(1:input$n_ph_extr_def, function(i) {
        eval(parse(text = paste0("output$ph_extr_", i, "<-renderUI({
          a<-i
          chk_Phases<-try(eval(parse(text=paste0(\'input$ph_extr_phs_\',a))),silent=TRUE)
					if(!is.null(chk_Phases)){
					if(!chk_Phases==\"\"){
					     phases<-gsub(\'\"\',\"\",break_on_comma(chk_Phases))
						fixedRow(
              lapply(1:(length(phases)), function(j) {
                chk_in_ph<-try(ph_extr_definitions_r$data[[a]][phases[j]],silent=TRUE)
                if(class(chk_in_ph)==\"try-error\"){chk_in_ph<-NULL}
                column(3,textInput(paste0(\'ph_extr_phs_\',a,\'_\',sub_brackets(phases[j])),phases[j],value=chk_in_ph))
              })
            )
          }}
        })")))
      })
    }
  })
  # Dynamically use solution model file to build phase models selection
  output$solution_models <- renderUI({
    if (input$solution_models_file == "load") {
      if (all(solution_models_file_r$data == "") | is.null(solution_models_file_r$data)) {
        solution_models_file <- ""
      } else {
        solution_models_file <- solution_models_file_r$data
      }
      updateTextInput(session, "solution_models_file", value = solution_models_file)
    }
    if (!(is.null(input$solution_models_file) | input$solution_models_file == "")) {
      if (!input$solution_models_file == "load") {
        validate(if (!dir.exists(gsub("/code", "/data", getwd()))) {
          return("Error: data directory does not exist")
        } else {
          NULL
        })
        validate(if (!file.exists(paste0(gsub("/code", "/data", getwd()), "/", input$solution_models_file))) {
          return("Error: Solution model file does not exist")
        } else {
          NULL
        })
        selectizeInput("use_sol_models", "Solution models", solution_models_available_r(), selected = use_sol_models_r$data, multiple = TRUE)
      }
    }
  })
  # Dynamically use number of component packet definitions to create the correct number of Component isolation,Phases inputs
  output$cp_packet_ui <- renderUI({
    if (input$cp_components == "load") {
      if (all(cp_components_r$data == "") | is.null(cp_components_r$data)) {
        updateTextInput(session, "cp_components", value = "")
      } else {
        updateTextInput(session, "cp_components", value = paste0(names(cp_components_r$data), collapse = ","))
      }
    }
    fixedRow(
      lapply(1:(as.numeric(length(break_on_comma(input$cp_components))) * 5), function(i) {
        a <- 1
        ii <- i
        while (ii > 5) {
          ii <- ii - 5
          a <- a + 1
        }
        if (ii == 1) {
          chk_cp_isolate <- try(cp_components_r$data[[a]], silent = TRUE)
          if (class(chk_cp_isolate) == "try-error") {
            chk_cp_isolate <- NULL
          }
          column(4, textInput(
            paste0("cp_isolate_", a),
            paste0("Proportion of component to isolate in ", break_on_comma(input$cp_components)[a]),
            value = chk_cp_isolate
          ))
        } else {
          if (ii == 2) {
            chk_cp_phases <- try(paste0(names(eval(parse(text = paste0("cp_phases_", names(cp_components_r$data[a]))))), collapse = ","), silent = TRUE)
            if (class(chk_cp_phases) == "try-error") {
              chk_cp_phases <- NULL
            }
            column(4, textInput(paste0("cp_phases_", a), "Phases to partition component into", value = chk_cp_phases))
          } else {
            if (ii == 3) {
              column(12, uiOutput(paste0("cp_packet_ui_ph_", a)))
            }
          }
        }
      })
    )
  })
  # auto create component packet inputs given number of phases
  observe({
    if (all(
      !input$cp_components == "load",
      !is.null(input$cp_components),
      !input$cp_components == "",
      !input$cp_components == 0
    )) {
      lapply(1:length(break_on_comma(input$cp_components)), function(i) {
        eval(parse(text = paste0("output$cp_packet_ui_ph_", i, "<-renderUI({
          a<-i
          chk_Phases<-try(eval(parse(text=paste0(\'input$cp_phases_\',a))),silent=TRUE)
					if(!is.null(chk_Phases)){
            if(!chk_Phases==\"\"){
              phases<-unlist(strsplit(chk_Phases,split=\",\"))
              fixedRow(
                lapply(1:(length(phases)), function(j) {
                  chk_in_cp_phases<-try(cp_phases_r$data[[a]][j],silent=TRUE)
                  if(class(chk_in_cp_phases)==\"try-error\"){chk_in_cp_phases<-NULL}
                  column(4,
                    tags$div(title=\"Input as a percentage, numeric value, equation based on chemistry of phases, or type 'excess' to partition remainder in packet.\",
                    textInput(paste0(\'cp_phases_\',a,\'_\',phases[j]),phases[j],value=chk_in_cp_phases))
                  )
                }
              ))
            }
          }
        })")))
      })
    }
  }) # end of component packet UI section

  # Dyanmically use number of phases to rename to create the correct number of inputs
  output$rename_phases_inputs <- renderUI({
    if (input$n_phases_to_rename == "load") {
      if (all(phases_to_rename_r$data == "") | is.null(phases_to_rename_r$data)) {
        def_num <- ""
      } else {
        def_num <- length(phases_to_rename_r$data)
      }
      updateTextInput(session, "n_phases_to_rename", value = def_num)
    } else {
      def_num <- input$n_phases_to_rename
    }
    if (!(is.null(def_num) | def_num == "" | def_num == 0)) {
      validate(if (is.na(suppressWarnings(as.numeric(def_num)))) {
        "Error: Number of phases to rename must be numeric"
      } else {
        NULL
      })
      validate(if (as.numeric(def_num) < 1) {
        "Error: Number of phases to rename must be greater than 0"
      } else {
        NULL
      })
      validate(if (!as.numeric(def_num) %% 1 == 0) {
        "Error: Number of phases to rename must be a whole number"
      } else {
        NULL
      })
      fixedRow(
        lapply(1:(as.numeric(def_num) * 3), function(i) {
          a <- 1
          ii <- i
          while (ii > 3) {
            ii <- ii - 3
            a <- a + 1
          }
          if (ii == 1) {
            chk_old <- try(unlist(strsplit(names(phases_to_rename_r$data)[a], split = "~"))[1], silent = TRUE)
            if (class(chk_old) == "try-error") {
              chk_old <- NULL
            }
            column(2, textInput(paste0("old_name_", a), "Old name", value = chk_old))
          } else {
            if (ii == 2) {
              chk_new <- try(unlist(strsplit(names(phases_to_rename_r$data)[a], split = "~"))[2], silent = TRUE)
              if (class(chk_new) == "try-error") {
                chk_new <- NULL
              }
              column(2, textInput(paste0("new_name_", a), "New name", value = chk_new))
            } else {
              if (ii == 3) {
                chk_condition <- try(phases_to_rename_r$data[a], silent = TRUE)
                if (class(chk_condition) == "try-error") {
                  chk_condition <- NULL
                }
                column(2, textInput(paste0("condition_", a), "Condition", value = chk_condition))
              } else {
              }
            }
          }
        })
      )
    }
  })

  # Dyanmically use abundance phases to select which phases to show
  output$select_abundance_phases <- renderUI({
    selectizeInput(
      "show_abundance_phases",
      "Show phases",
      c(
        abundance_phases_available_r(),
        "Reactive Subsystem",
        "Extract Subsystem",
        "Cumulative Extract Subsystem",
        "Full System"
      ),
      multiple = TRUE,
      selected = "Reactive Subsystem"
    )
  })
  abundance_phases_available_r <- reactive({
    if (!is.null(input$axis)) {
      rownames(switch(input$axis,
        x = phase_abundance(store_r$crust_r, input$axis, as.numeric(input$path_y), input$start_x, input$end_x, input$path_label, input_pt = store_r$input_pt_r, input$proportion),
        y = phase_abundance(store_r$crust_r, input$axis, as.numeric(input$path_x), input$start_y, input$end_y, input$path_label, input_pt = store_r$input_pt_r, input$proportion)
      )[[2]])
    }
  })
  # Create Compilation button
  observeEvent(input$create_compilation, {
    # Compile PAM fields (corelation to other PAMs)
    validate(need(!input$PAM_compilation == "", "Please provide a comma seperated list of working files in PAM compilation"))
    cat("Creating Compilation Legend for", input$PAM_compilation, "\n")
    flush.console()
    # mergers must contain the current file
    mergers <- sort(union(unlist(strsplit(input$PAM_compilation, ",")), input$working_file))
    if (!dir.exists(paste0(input$projects_directory, "/Compile"))) {
      dir.create(paste0(input$projects_directory, "/Compile"))
    }
    compile_names <- NULL
    for (i in mergers) {
      attach(paste0(input$projects_directory, "/", i, "/", i, ".RData"), warn.conflicts = FALSE)
      pull_crust <- get("crust", which(search() == paste0("file:", input$projects_directory, "/", i, "/", i, ".RData")))
      detach(pos = which(search() == paste0("file:", input$projects_directory, "/", i, "/", i, ".RData")))
      compile_names <- union(compile_names, get_PAM_names(neaten_crust(pull_crust, input$phase_aliases), input$PAM_system)[[1]])
    }
    write.table(
      compile_names,
      paste0(input$projects_directory, "/Compile/", input$PAM_compilation, " compilation legend.txt"),
      sep = "\t",
      quote = F,
      col.names = FALSE
    )
    cat("Compilation successfully created for", input$PAM_compilation, "\n")
    flush.console()
  })
  # Save Data button
  observeEvent(input$save_data, {
    switch(input$output_type,
      "Data File" =
        reactive_message$data <- write_data_file(
          data_file(crust_out(),
            x_n = length(crust_out()[[1]]),
            y_n = length(crust_out()),
            input$choose_columns,
            input$choose_rows,
            input$choose_points
          ),
          input$working_file,
          input$projects_directory,
          input$file_type
        ),
      "Grid" =
        if (TRUE) {
          reactive_message$data <- "Saving"
          reactive_message$data <- draw_Grid_r()
        },
      "Phase Abundance Along Path" =
        if (TRUE) {
          reactive_message$data <- "Saving"
          reactive_message$data <- draw_abundance_r()
        },
      "PAM" =
        if (TRUE) {
          reactive_message$data <- "Saving"
          reactive_message$data <- draw_PAM_r()
        }
    )
  })
  # Send to GCDkit button
  observeEvent(input$send_gcdkit, {
    # is sourcing necessary here?
    source("Rcrust_functions.r")
    crust_to_gcdkit(store_r$crust_r, input$choose_columns, input$choose_rows, input$choose_points, GCDkitGUI = FALSE)
  })
  # Assign labels button
  observeEvent(input$assign_label, {
    source("Rcrust_functions.r")
    crust <- assign_label(crust_out(), input$from_label, input$to_label, input$label_name, input$label_value)
    store_r$crust_r <- crust
  })
  solution_models_available_r <- reactive({
    solution_models_available <- NULL
    qq <- scan(file = paste0(gsub("/code", "/data", getwd()), "/", input$solution_models_file), what = "character", sep = "\n", quiet = TRUE)
    foo <- grep("begin_model", qq)
    for (i in 1:length(foo)) {
      if (i == length(foo)) {
        toto <- substr(qq[foo[i]:length(qq)], 1, 1)
      } else {
        toto <- substr(qq[foo[i]:foo[i + 1]], 1, 1)
      }
      modrow <- grep("[A-Za-z]", toto)[2]
      line <- qq[foo[i] + modrow - 1]
      solution_models_available <- c(solution_models_available, strsplit(line, " ")[[1]][1])
    }
    return(sort(solution_models_available))
  })

  # Accepts custom definition of the form grid_data("H2O","Bulk_rs",crust_out())
  # or of the form list("K2O/Na2O",grid_data("K2O","Melt_rs",crust_out())[[2]]/grid_data("Na2O","Melt_rs",crust_out())[[2]])
  # or of the form list("K2O/Na2O",grid_data("K2O","Melt_rs",crust_out())[[2]]/grid_data("Na2O","Melt_rs",crust_out())[[2]])
  grid_data_r <- reactive({
    if (input$Grid_variable == "Custom") {
      grid_in <- eval(parse(text = input$Custom_selection))
    } else {
      grid_in <- grid_data(input$Grid_variable, input$Grid_variable_phase, crust_out(), input_pt)
      # mod tag allow this sort of reading below in the GUI
      # grid_in[[2]]<-as.matrix(read.csv("H2O saturation.csv"))
    }
    # Remove values
    if (input$remove_values != "") {
      for (val in unlist(strsplit(input$remove_values, split = ",|;"))) {
        for (i in 1:length(grid_in[[2]])) {
          if (!is.na(grid_in[[2]][i])) {
            if (length(unlist(strsplit(val, split = "-"))) > 1) {
              if (grid_in[[2]][i] >= unlist(strsplit(val, split = "-"))[1] & grid_in[[2]][i] <= unlist(strsplit(val, split = "-"))[2]) {
                grid_in[[2]][i] <- NA
              }
            } else {
              if (grid_in[[2]][i] == val) {
                grid_in[[2]][i] <- NA
              }
            }
          }
        }
      }
    }
    # apply matrix rotations
    if (input$rotation != 0) {
      for (i in 1:input$rotation) {
        grid_in[[2]] <- rotate(grid_in[[2]])
      }
    }
    # apply matrix reflections
    if (any(input$reflection == "Horizontal")) {
      grid_in[[2]] <- flip_y(grid_in[[2]])
    }
    if (any(input$reflection == "Vertical")) {
      grid_in[[2]] <- flip_x(grid_in[[2]])
    }
    validate(need(nrow(grid_in[[2]]) > 1, "Cannot plot a grid with only 1 row"))
    validate(need(ncol(grid_in[[2]]) > 1, "Cannot plot a grid with only 1 column"))
    if (any(input$Grid_axes == "bottom")) {
      bottom <- create_Axes(input$Grid_bottom_axis, input$Grid_bottom_axis_grid_phase, input$Grid_bottom_axis_increments, "x", crust_out(), input_pt)
    } else {
      bottom <- NULL
    }
    if (any(input$Grid_axes == "left")) {
      left <- create_Axes(input$Grid_left_axis, input$Grid_left_axis_grid_phase, input$Grid_left_axis_increments, "y", crust_out(), input_pt)
    } else {
      left <- NULL
    }
    if (any(input$Grid_axes == "top")) {
      top <- create_Axes(input$Grid_top_axis, input$Grid_top_axis_grid_phase, input$Grid_top_axis_increments, "x", crust_out(), input_pt, y_n)
    } else {
      top <- NULL
    }
    if (any(input$Grid_axes == "right")) {
      right <- create_Axes(input$Grid_right_axis, input$Grid_right_axis_grid_phase, input$Grid_right_axis_increments, "y", crust_out(), input_pt, x_n)
    } else {
      right <- NULL
    }
    # mod-tag: look into rotations
    # Apply rotations
    # if(any(input$Grid_axes=="bottom")){colnames(grid_in[[2]])<-paste0(create_Axes(input$Grid_bottom_axis,input$Grid_bottom_axis_grid_phase,length(crust[[1]])+1,"x",crust_out(),input_pt)[[2]][-(length(crust[[1]])+1)],"_{",1:length(crust[[1]]),";1}")}
    # if(any(input$Grid_axes=="left")){rownames(grid_in[[2]])<-paste0(create_Axes(input$Grid_left_axis,input$Grid_left_axis_grid_phase,length(crust)+1,"y",crust_out(),input_pt)[[2]][-(length(crust)+1)],"_{1;",1:length(crust),"}")}
    return(list(grid_in[[1]], grid_in[[2]], bottom, left, top, right))
  })
  phase_abundance_r <- reactive({
    if (!is.null(input$axis)) {
      switch(input$axis,
        x = phase_abundance(crust_out(), input$axis, as.numeric(input$path_y), input$start_x, input$end_x, input$path_label, input_pt = store_r$input_pt_r, input$proportion),
        y = phase_abundance(crust_out(), input$axis, as.numeric(input$path_x), input$start_y, input$end_y, input$path_label, input_pt = store_r$input_pt_r, input$proportion)
      )
    }
  })
  PAM_r <- reactive({
    PAM_calc(crust_out(), input$PAM_system, input$compile_PAM, input$PAM_compilation)
  })
  neaten_crust <- function(crust_neat = NULL, phase_aliases = NULL) {
    # rename using aliases if given
    if (!phase_aliases == "") {
      split_aliases <- strsplit(strsplit(phase_aliases, c(","))[[1]], "=")
      phase_aliases <- NULL
      phase_aliases_names <- NULL
      for (i in 1:length(split_aliases)) {
        phase_aliases <- c(phase_aliases, split_aliases[[i]][1])
        phase_aliases_names <- c(phase_aliases_names, split_aliases[[i]][2])
      }
      names(phase_aliases) <- phase_aliases_names
      # seperate merge commands
      merge_aliases <- grep("&", phase_aliases, value = TRUE)
      if (length(merge_aliases) != 0) {
        phase_aliases <- phase_aliases[-grep("&", phase_aliases)]
      }
      merge_aliases <- strsplit(merge_aliases, "&")
      if (length(merge_aliases) == 0) {
        merge_aliases <- ""
      }
      if (length(phase_aliases) == 0) {
        phase_aliases <- ""
      }
      if (all(phase_aliases != "") | all(merge_aliases != "")) {
        for (x_i in 1:length(crust[[1]])) {
          for (y_i in 1:length(crust)) {
            delete <- NULL
            grab <- rownames(crust_neat[[y_i]][[x_i]])
            if (length(grab) > 0) {
              grab_split <- strsplit(grab, "_")
              grab_sys <- grab_name <- grab
              for (ph in 1:length(grab)) {
                grab_name[ph] <- grab_split[[ph]][1]
                # if first part of system is numeric extract up to end of _
                if (!is.na(suppressWarnings(as.numeric(grab_split[[ph]][-1][1])))) {
                  grab_sys[ph] <- paste(grab_split[[ph]][c(-1, -2)], collapse = "_")
                } else {
                  grab_sys[ph] <- paste(grab_split[[ph]][-1], collapse = "_")
                }
              }
              # rename using aliases
              if (all(phase_aliases != "")) {
                for (ph in 1:length(grab)) {
                  if (any(grab_name[ph] == phase_aliases)) {
                    grab_name[ph] <- names(which(grab_name[ph] == phase_aliases)[1])
                    # label unwanted phases (phases of zero mass or phases labelled as "hide")
                    if (crust_neat[[y_i]][[x_i]][ph, "mass"] == 0 | grab_name[ph] == "hide") {
                      delete <- c(delete, ph)
                    }
                  }
                }
                # number duplicates within systems
                renamed <- paste(grab_name, grab_sys, sep = "_")
                num <- 1
                while (any(duplicated(renamed))) {
                  renamed[which(duplicated(renamed))] <- paste(grab_name[which(duplicated(renamed))], num, grab_sys[which(duplicated(renamed))], sep = "_")
                  num <- num + 1
                }
                rownames(crust_neat[[y_i]][[x_i]]) <- renamed
              }
              # merge if required
              if (all(merge_aliases != "")) {
                systems <- c("rs", "es", "cumul")
                for (sys in systems) {
                  for (merge_try in 1:length(merge_aliases)) {
                    merge_nos <- NULL
                    for (i in 1:length(merge_aliases[[merge_try]])) {
                      if (any(paste(grab_name, grab_sys, sep = "_") == paste(merge_aliases[[merge_try]][i], sys, sep = "_"))) {
                        merge_nos <- union(merge_nos, which(grab_name == merge_aliases[[merge_try]][i]))
                      }
                    }
                    merger <- NULL
                    if (!is.null(merge_nos)) {
                      for (i in merge_nos) {
                        merger <- rbind(merger, crust_neat[[y_i]][[x_i]][i, , drop = FALSE])
                      }
                      merged <- .wtd.add(merger, avname = names(merge_aliases)[merge_try])
                      crust_neat[[y_i]][[x_i]][merge_nos[1], ] <- merged
                      rownames(crust_neat[[y_i]][[x_i]])[merge_nos[1]] <- paste(names(merge_aliases)[merge_try], sys, sep = "_")
                      if (length(merge_nos) > 1) {
                        delete <- c(delete, merge_nos[-1])
                      }
                    }
                  }
                }
              }
              # remove unwanted phases
              if (!is.null(delete)) {
                crust_neat[[y_i]][[x_i]] <- crust_neat[[y_i]][[x_i]][-delete, ]
              }
            }
          }
        }
      }
    }
    return(crust_neat)
  }
  crust_out <- reactive({
    return(neaten_crust(store_r$crust_r, input$phase_aliases))
  })
  create_Axes <- function(axes_variable, axes_variable_phase, axes_increments, axes_direction, crust = NULL, input_pt = NULL, axes_line = 1) {
    # Given axes_variable,axes_variable_phase,axes_increments,axes_direction create axis values
    if (axes_increments == "Increments") {
      axes_increments <- 11
    } else {
      axes_increments <- as.numeric(axes_increments)
    }
    if (axes_direction == "y") {
      i_n <- length(crust) / (axes_increments - 1) * (0:(axes_increments - 1)) + 1
    } else {
      i_n <- length(crust[[1]]) / (axes_increments - 1) * (0:(axes_increments - 1)) + 1
    }
    if (axes_variable == "y_i" | axes_variable == "x_i") {
      axis_values <- i_n[-length(i_n)]
    } else {
      axis_values <- NULL
      for (i in i_n[-length(i_n)]) {
        if (axes_direction == "y") {
          if (axes_variable == "Temperature (C)") {
            axis_values <- c(axis_values, input_pt[[i]][[axes_line]][, "Temperature"])
          } else if (axes_variable == "Pressure (kbar)") {
            axis_values <- c(axis_values, input_pt[[i]][[axes_line]][, "Pressure"])
          } else if (axes_variable == "Pressure (Mpa)") {
            axis_values <- c(axis_values, input_pt[[i]][[axes_line]][, "Pressure"] * 100)
          } else if (axes_variable == "Temperature (K)") {
            axis_values <- c(axis_values, input_pt[[i]][[axes_line]][, "Temperature"] + 273)
          } else {
            chk <- try(crust[[i]][[axes_line]][axes_variable_phase, axes_variable], silent = TRUE)
            if (class(chk) == "try-error") {
              axis_values <- c(axis_values, 0)
            } else {
              axis_values <- c(axis_values, crust[[i]][[axes_line]][axes_variable_phase, axes_variable])
            }
          }
        } else {
          if (axes_variable == "Temperature (C)") {
            axis_values <- c(axis_values, input_pt[[axes_line]][[i]][, "Temperature"])
          } else if (axes_variable == "Pressure (kbar)") {
            axis_values <- c(axis_values, input_pt[[axes_line]][[i]][, "Pressure"])
          } else if (axes_variable == "Pressure (Mpa)") {
            axis_values <- c(axis_values, input_pt[[axes_line]][[i]][, "Pressure"] * 100)
          } else if (axes_variable == "Temperature (K)") {
            axis_values <- c(axis_values, input_pt[[axes_line]][[i]][, "Temperature"] + 273)
          } else {
            chk <- try(crust[[axes_line]][[i]][axes_variable_phase, axes_variable], silent = TRUE)
            if (class(chk) == "try-error") {
              axis_values <- c(axis_values, 0)
            } else {
              axis_values <- c(axis_values, crust[[axes_line]][[i]][axes_variable_phase, axes_variable])
            }
          }
        }
      }
    }
    axis_values <- c(axis_values, axis_values[length(axis_values)] + axis_values[length(axis_values)] - axis_values[length(axis_values) - 1])
    if (axes_variable == "y_i" |
      axes_variable == "x_i" |
      axes_variable == "Temperature (C)" |
      axes_variable == "Pressure (kbar)" |
      axes_variable == "Pressure (Mpa)" |
      axes_variable == "Temperature (K)"
    ) {
      axis_title <- axes_variable
    } else {
      axis_title <- paste(axes_variable_phase, axes_variable)
    }
    return(list(axis_title, signif(axis_values, digits = 4)))
  }
  draw_Grid_r <- reactive({
    library(graphics)
    library(grDevices)
    validate(need(nrow(grid_data_r()[[2]]) > 1, "Cannot plot a grid with only 1 row"))
    validate(need(ncol(grid_data_r()[[2]]) > 1, "Cannot plot a grid with only 1 column"))
    bottom <- grid_data_r()[[3]]
    left <- grid_data_r()[[4]]
    top <- grid_data_r()[[5]]
    right <- grid_data_r()[[6]]
    action <- "View"
    if (!is.null(reactive_message$data)) {
      if (reactive_message$data == "Saving") {
        action <- "Save"
      }
    }
    if (action == "Save") {
      outfile_path <- paste0(input$projects_directory, "/", input$working_file, "/Outputs/", input$working_file, "_grid_", gsub("\\/", " per ", grid_data_r()[[1]]), input$file_type)
      if (input$file_type == ".txt") {
        write.table(grid_data_r()[[2]], outfile_path, sep = "\t", quote = F)
      }
      if (input$file_type == ".csv") {
        write_test <- try(write.csv(grid_data_r()[[2]], outfile_path), silent = TRUE)
        if (class(write_test) == "try-error") {
          cat("Error cannot write to ", outfile_path, ", please close all programs that may be accessing the file then try again\n")
          return(paste0("Error could not save Grid File: ", grid_data_r()[[1]], ", file may be open in another program, please close all programs that may be accessing the file then try again\n"))
        }
      }
      if (input$file_type == ".ps") {
        outfile_path <- gsub("%", "%%", outfile_path)
        postscript(file = outfile_path, onefile = TRUE, horizontal = TRUE)
      }
    }
    if (action == "View" | (action == "Save" & input$file_type == ".ps")) {
      if (input$Grid_colours == "RColorBrewer") {
        col_levels <- as.numeric(input$Grid_RColorBrewer_levels)
        col_brew <- RColorBrewer::brewer.pal(n = col_levels - 1, name = input$Grid_RColorBrewer_colours)
        if (input$Grid_RColorBrewer_reverse) {
          col_brew <- rev(col_brew)
        }
        col_pal <- colorRampPalette(col_brew)(col_levels)
        # 		mod tag - edit below line to allow GUI to set hard limits on z value
        # 		filled.contour(t(flip_y(grid_data_r()[[2]])),zlim=c(0,2.6),
        if (input$Grid_RColorBrewer_zlimits != "") {
          filled.contour(t(flip_y(grid_data_r()[[2]])),
            zlim = as.numeric(unlist(strsplit(input$Grid_RColorBrewer_zlimits, ","))),
            plot.axes = {
              if (!is.null(bottom)) {
                axis(1, (0:(length(bottom[[2]]) - 1)) / (length(bottom[[2]]) - 1), bottom[[2]])
              } else {
                NULL
              }
              if (!is.null(left)) {
                axis(2, (0:(length(left[[2]]) - 1)) / (length(left[[2]]) - 1), left[[2]])
              } else {
                NULL
              }
              if (!is.null(top)) {
                axis(3, (0:(length(top[[2]]) - 1)) / (length(top[[2]]) - 1), top[[2]])
              } else {
                NULL
              }
              if (!is.null(right)) {
                axis(4, (0:(length(right[[2]]) - 1)) / (length(right[[2]]) - 1), right[[2]])
              } else {
                NULL
              }
            },
            nlevels = col_levels, col = col_pal
          )
        } else {
          filled.contour(t(flip_y(grid_data_r()[[2]])),
            plot.axes = {
              if (!is.null(bottom)) {
                axis(1, (0:(length(bottom[[2]]) - 1)) / (length(bottom[[2]]) - 1), bottom[[2]])
              } else {
                NULL
              }
              if (!is.null(left)) {
                axis(2, (0:(length(left[[2]]) - 1)) / (length(left[[2]]) - 1), left[[2]])
              } else {
                NULL
              }
              if (!is.null(top)) {
                axis(3, (0:(length(top[[2]]) - 1)) / (length(top[[2]]) - 1), top[[2]])
              } else {
                NULL
              }
              if (!is.null(right)) {
                axis(4, (0:(length(right[[2]]) - 1)) / (length(right[[2]]) - 1), right[[2]])
              } else {
                NULL
              }
            },
            nlevels = col_levels, col = col_pal
          )
        }
      } else {
        filled.contour(t(flip_y(grid_data_r()[[2]])),
          plot.axes = {
            if (!is.null(bottom)) {
              axis(1, (0:(length(bottom[[2]]) - 1)) / (length(bottom[[2]]) - 1), bottom[[2]])
            } else {
              NULL
            }
            if (!is.null(left)) {
              axis(2, (0:(length(left[[2]]) - 1)) / (length(left[[2]]) - 1), left[[2]])
            } else {
              NULL
            }
            if (!is.null(top)) {
              axis(3, (0:(length(top[[2]]) - 1)) / (length(top[[2]]) - 1), top[[2]])
            } else {
              NULL
            }
            if (!is.null(right)) {
              axis(4, (0:(length(right[[2]]) - 1)) / (length(right[[2]]) - 1), right[[2]])
            } else {
              NULL
            }
          },
          color.palette = eval(parse(text = input$Grid_colours))
        )
      }
      # Colour options = gray.colors,heat.colours,terrain.colours,rainbow,topo.colours
      # look into http://geog.uoregon.edu/datagraphics/color_scales.htm
      # colorRampPalette(rev(brewer.pal(n = 7, name = "BuPu")))(100)
      # graphics::filled.contour(t(flip_y(grid_data_r()[[2]])),nlevels=10,col=colorRampPalette(rev(brewer.pal(n = 9, name = "YlOrRd")))(10))
      if (!is.null(bottom)) {
        mtext(bottom[[1]], 1, line = 3)
      }
      if (!is.null(left)) {
        mtext(left[[1]], 2, line = 3)
      }
      if (!is.null(top)) {
        mtext(top[[1]], 3, line = 3)
      }
      if (!is.null(right)) {
        mtext(right[[1]], 4, line = 3)
      }
      if (action == "Save") {
        title(paste("Grid of ", gsub("\\/", " per ", grid_data_r()[[1]]), " for ", input$working_file))
        dev.off()
      }
    }
    detach("package:graphics")
    detach("package:grDevices")
    if (action == "Save") {
      cat("File written to ", outfile_path, "\n")
      return(paste0(paste("Grid of", gsub("\\/", " per ", grid_data_r()[[1]])), " saved to ", input$projects_directory, "/", input$working_file, "/Outputs/\n"))
    }
  })
  draw_abundance_r <- reactive({
    # Plot phase abundance versus path cumulative column graph
    library(graphics)
    library(RColorBrewer)
    library(grDevices)
    # All RColorBrewer palettes display.brewer.all()
    # Display a specific pallette display.brewer.pal(12,"Set3")
    action <- "View"
    if (!is.null(reactive_message$data)) {
      if (reactive_message$data == "Saving") {
        action <- "Save"
      }
    }
    if (action == "Save") {
      outfile_path <- paste0(input$projects_directory, "/", input$working_file, "/Outputs/", input$working_file, "_phase_abundance_", gsub("\\/", " per ", phase_abundance_r()[[1]]), input$file_type)
      if (input$file_type == ".txt") {
        write_phase_abundance(phase_abundance_r(), input$working_file, input$projects_directory, input$file_type)
      }
      if (input$file_type == ".csv") {
        write_test <- try(write_phase_abundance(phase_abundance_r(), input$working_file, input$projects_directory, input$file_type), silent = TRUE)
        if (class(write_test) == "try-error") {
          cat("Error cannot write file please close all programs that may be accessing the file then try again\n")
          return(paste0("Error could not save Phase abundance File: file may be open in another program, please close all programs that may be accessing the file then try again\n"))
        }
      }
      if (input$file_type == ".ps") {
        cat("Error Functionality not available yet\n")
        return(paste0("Error Functionality not available yet\n"))
        # mod-tag: look into saving phase abundance directly to .ps
        # outfile_path<-gsub("%","%%",outfile_path)
        # postscript(file=outfile_path,onefile=TRUE,horizontal=TRUE)
      }
    }
    if (action == "View" | (action == "Save" & input$file_type == ".ps")) {
      phase_abundance_data <- phase_abundance_r()[[2]]
      # subset phases
      if (!is.null(input$show_abundance_phases)) {
        phases_to_show <- input$show_abundance_phases
        phases_present <- rownames(phase_abundance_data)
        if (any(phases_to_show == "Reactive Subsystem")) {
          phases_to_show <- union(phases_to_show, phases_present[setdiff(which(substrRight(phases_present, 3) == "_rs"), which(phases_present == "Bulk_rs"))])
          phases_to_show <- phases_to_show[-grep("Reactive Subsystem", phases_to_show)]
        }
        if (any(phases_to_show == "Extract Subsystem")) {
          phases_to_show <- union(phases_to_show, phases_present[setdiff(which(substrRight(phases_present, 3) == "_es"), which(phases_present == "Bulk_es"))])
          phases_to_show <- phases_to_show[-grep("Extract Subsystem", phases_to_show)]
        }
        if (any(phases_to_show == "Cumulative Extract Subsystem")) {
          phases_to_show <- union(phases_to_show, phases_present[setdiff(which(substrRight(phases_present, 9) == "_es_cumul"), which(phases_present == "Bulk_es_cumul"))])
          phases_to_show <- phases_to_show[-grep("Cumulative Extract Subsystem", phases_to_show)]
        }
        if (any(phases_to_show == "Full System")) {
          phases_to_show <- union(phases_to_show, phases_present)
          phases_to_show <- phases_to_show[-grep("Full System", phases_to_show)]
        }
        phase_abundance_data <- phase_abundance_data[match(phases_to_show, rownames(phase_abundance_data)), , drop = FALSE]
      }
      # Normalise to percentage
      phase_abundance_data <- apply(phase_abundance_data, 2, function(x) {
        x * 100 / sum(x, na.rm = T)
      })
      # Set colour pallette
      cols <- RColorBrewer::brewer.pal(min(nrow(phase_abundance_data), 12), "Set3")
      # Plot
      barplot(
        phase_abundance_data,
        space = 0,
        col = cols,
        border = NA,
        legend.text = (input$legend != "None"),
        args.legend = list(x = input$legend, bty = "n"),
        xlab = input$path_label,
        ylab = paste("Phase abundance", input$proportion)
      )
      if (action == "Save") {
        title(paste("Grid of ", gsub("\\/", " per ", grid_data_r()[[1]]), " for ", input$working_file))
        dev.off()
      }
    }
    detach("package:graphics")
    detach("package:grDevices")
    if (action == "Save") {
      cat("File written to ", outfile_path, "\n")
      return(paste0(paste("Grid of", gsub("\\/", " per ", grid_data_r()[[1]])), " saved to ", input$projects_directory, "/", input$working_file, "/Outputs/\n"))
    }
  })
  draw_PAM_r <- reactive({
    # Bind all_pres to PAM legend
    PAM_legend <- c(PAM_r()[[3]], names(PAM_r()[[2]]))
    names(PAM_legend) <- c("\"All fields are +\"", as.numeric(PAM_r()[[2]]))
    compiled_legend <- c(PAM_r()[[3]], names(PAM_r()[[6]]))
    names(compiled_legend) <- c("\"All fields are +\"", as.numeric(PAM_r()[[6]]))
    action <- "View"
    if (!is.null(reactive_message$data)) {
      if (reactive_message$data == "Saving") {
        action <- "Save"
      }
    }
    if (action == "Save") {
      outfile_path <- paste0(input$projects_directory, "/", input$working_file, "/Outputs/", input$working_file, " PAM", input$file_type)
      outfile_legend_path <- paste0(input$projects_directory, "/", input$working_file, "/Outputs/", input$working_file, " PAM legend", input$file_type)
      outfile_compile_legend_path <- paste0(
        input$projects_directory,
        "/Compile/",
        paste(sort(union(unlist(strsplit(input$PAM_compilation, ",")), input$working_file)), collapse = ","),
        " compile legend",
        input$file_type
      )
      if (input$file_type == ".txt") {
        write.table(PAM_r()[[1]], outfile_path, sep = "\t", quote = F, row.names = TRUE)
        write.table(PAM_legend, outfile_legend_path, sep = "\t", quote = F, row.names = TRUE, col.names = FALSE)
        if (input$compile_PAM) {
          write.table(compiled_legend, outfile_compile_legend_path, sep = "\t", quote = F, row.names = TRUE, col.names = FALSE)
        }
      }
      if (input$file_type == ".csv") {
        write_test <- try(write.csv(PAM_r()[[1]], outfile_path, row.names = TRUE), silent = TRUE)
        if (class(write_test) == "try-error") {
          cat("Error cannot write to ", outfile_path, ", please close all programs that may be accessing the file then try again\n")
          return(paste0("Error could not save phase assemblage map: ", outfile_path, ", file may be open in another program, please close all programs that may be accessing the file then try again\n"))
        }
        write_test <- try(write.csv(PAM_legend, outfile_legend_path, row.names = TRUE, col.names = FALSE), silent = TRUE)
        if (class(write_test) == "try-error") {
          cat("Error cannot write to ", outfile_legend_path, ", please close all programs that may be accessing the file then try again\n")
          return(paste0("Error could not save phase assemblage map legend: ", outfile_legend_path, ", file may be open in another program, please close all programs that may be accessing the file then try again\n"))
        }
        if (input$compile_PAM) {
          write_test <- try(write.csv(compiled_legend, outfile_compile_legend_path, row.names = TRUE, col.names = FALSE), silent = TRUE)
          if (class(write_test) == "try-error") {
            cat("Error cannot write to ", outfile_compile_legend_path, ", please close all programs that may be accessing the file then try again\n")
            return(paste0("Error could not save phase assemblage map legend: ", outfile_compile_legend_path, ", file may be open in another program, please close all programs that may be accessing the file then try again\n"))
          }
        }
      }
      if (input$file_type == ".ps") {
        library(grDevices)
        postscript(file = outfile_path, onefile = TRUE, horizontal = TRUE)
      }
    }
    if (action == "View" | (action == "Save" & input$file_type == ".ps")) {
      library(graphics)
      par(mar = c(5, 4, 4, 5) + .1)
      library(grDevices)
      library(raster)
      library(rgeos)
      max_pol <- max(length(PAM_r()[[2]]), length(PAM_r()[[6]]))
      raster::plot(PAM_r()[[4]], col = gray(1:max_pol / max_pol)[PAM_r()[[5]]])
      if (input$PAM_labels == "Numbers") {
        raster::text(gCentroid(PAM_r()[[4]], byid = TRUE), labels = PAM_r()[[5]], col = "Black")
      } else {
        if (input$compile_PAM) {
          raster::text(gCentroid(PAM_r()[[4]], byid = TRUE), labels = compiled_legend[as.character(PAM_r()[[5]])], col = "Black")
        } else {
          raster::text(gCentroid(PAM_r()[[4]], byid = TRUE), labels = PAM_legend[as.character(PAM_r()[[5]])], col = "Black")
        }
      }
      # Create axes
      if (length(input$PAM_axes) > 0) {
        for (j in 1:length(input$PAM_axes)) {
          # get variable and increment choice
          var_choice <- eval(parse(text = paste("input$PAM_", input$PAM_axes[j], "_axis", sep = "")))
          if (!(var_choice == "y_i" | var_choice == "x_i" | var_choice == "Temperature" | var_choice == "Pressure")) {
            grid_phase_choice <- eval(parse(text = paste("input$PAM_", input$PAM_axes[j], "_axis_grid_phase", sep = "")))
          }
          increment_choice <- eval(parse(text = paste("input$PAM_", input$PAM_axes[j], "_axis_increments", sep = "")))
          if (increment_choice == "Increments") {
            increment_choice <- 11
          } else {
            increment_choice <- as.numeric(increment_choice)
          }
          if (input$PAM_axes[j] == "left" | input$PAM_axes[j] == "right") {
            i_n <- length(crust) / (increment_choice - 1) * (0:(increment_choice - 1)) + 1
          } else {
            i_n <- length(crust[[1]]) / (increment_choice - 1) * (0:(increment_choice - 1)) + 1
          }
          if (var_choice == "y_i" | var_choice == "x_i") {
            axis_values <- i_n[-length(i_n)]
          } else {
            axis_values <- NULL
            for (i in i_n[-length(i_n)]) {
              if (input$PAM_axes[j] == "left" | input$PAM_axes[j] == "right") {
                if (var_choice == "Temperature" | var_choice == "Pressure") {
                  axis_values <- c(axis_values, input_pt[[i]][[1]][, var_choice])
                } else {
                  chk <- try(crust[[i]][[1]][grid_phase_choice, var_choice], silent = TRUE)
                  if (class(chk) == "try-error") {
                    axis_values <- c(axis_values, 0)
                  } else {
                    axis_values <- c(axis_values, crust[[i]][[1]][grid_phase_choice, var_choice])
                  }
                }
              } else {
                if (var_choice == "Temperature" | var_choice == "Pressure") {
                  axis_values <- c(axis_values, input_pt[[1]][[i]][, var_choice])
                } else {
                  chk <- try(crust[[1]][[i]][grid_phase_choice, var_choice], silent = TRUE)
                  if (class(chk) == "try-error") {
                    axis_values <- c(axis_values, 0)
                  } else {
                    axis_values <- c(axis_values, crust[[1]][[i]][grid_phase_choice, var_choice])
                  }
                }
              }
            }
          }
          axis_values <- c(axis_values, axis_values[length(axis_values)] + axis_values[length(axis_values)] - axis_values[length(axis_values) - 1])
          # round to 2 significant figures
          axis_values <- signif(axis_values, digits = 4)
          side_no <- switch(input$PAM_axes[j],
            "bottom" = 1,
            "left" = 2,
            "top" = 3,
            "right" = 4
          )
          axis(side_no, (0:(increment_choice - 1)) / (increment_choice - 1), axis_values)
          if (var_choice == "y_i" | var_choice == "x_i" | var_choice == "Temperature" | var_choice == "Pressure") {
            mtext(var_choice, side = side_no, line = 3)
          } else {
            mtext(paste(grid_phase_choice, var_choice), side = side_no, line = 3)
          }
        }
      }
      if (!PAM_r()[[3]] == "") {
        title(, paste("All fields are +", PAM_r()[[3]]))
      }
      # Add contour
      if (!input$PAM_contour == "None") {
        grid_out_mat <- matrix(0, y_n, x_n)
        if (input$PAM_contour == "Pressure" | input$PAM_contour == "Temperature") {
          for (x_i in 1:x_n) {
            for (y_i in 1:y_n) {
              grid_out_mat[y_i, x_i] <- input_pt[[y_i]][[x_i]][, input$PAM_contour]
            }
          }
        } else {
          for (x_i in 1:x_n) {
            for (y_i in 1:y_n) {
              chk <- try(crust_out()[[y_i]][[x_i]][input$PAM_contour_grid_phase, input$PAM_contour], silent = TRUE)
              if (class(chk) == "try-error") {
                grid_out_mat[y_i, x_i] <- 0
              } else {
                grid_out_mat[y_i, x_i] <- chk
              }
            }
          }
        }
        x <- raster::raster(flip_y(grid_out_mat))
        if (input$PAM_contour_increments == "Default Increments") {
          raster::contour(x, add = TRUE, col = "red")
        } else {
          if (input$PAM_contour_increments == "In/Out") {
            raster::contour(x, add = TRUE, levels = 0.0000000000000000001, col = "red")
          } else {
            raster::contour(x, add = TRUE, nlevels = as.numeric(input$PAM_contour_increments), col = "red")
          }
        }
      }
      if (action == "Save") {
        title(paste("Phase Assemblage Map for", input$working_file), )
        dev.off()
        write.table(PAM_legend, outfile_legend_path, sep = "\t", quote = F, row.names = TRUE, col.names = FALSE)
        if (input$compile_PAM) {
          write.table(compiled_legend, outfile_compile_legend_path, sep = "\t", quote = F, row.names = TRUE, col.names = FALSE)
        }
      }
      detach("package:graphics")
      detach("package:rgeos")
      detach("package:raster")
      detach("package:sp")
      detach("package:grDevices")
    }
    if (action == "Save") {
      cat("File written to ", outfile_path, "\n")
      cat("File written to ", outfile_legend_path, "\n")
      return(paste0("Phase Assemblage Map saved to ", input$projects_directory, "/", input$working_file, "/Outputs/\n"))
    }
  })
  output$output_header <- renderText(
    switch(input$output_type,
      "Data File" = if (!is.null(store_r$crust_r)) {
        "Compilation data file"
      },
      "Grid" = if (!is.null(grid_data_r())) {
        paste0(grid_data_r()[[1]], " on (X,Y) grid")
      },
      "Phase Abundance Along Path" = if (!is.null(phase_abundance_r())) {
        phase_abundance_r()[[1]]
      },
      "PAM" = if (!is.null(store_r$crust_r)) {
        "Phase Assemblage Map"
      }
    )
  )
  # Dyanmically create output for viewing
  output$output_view <- renderUI({
    if (!is.null(input$output_form)) {
      switch(input$output_form,
        "Data" = tableOutput("table"),
        "Legend" = tableOutput("table"),
        switch(input$output_type,
          "PAM" = plotOutput("plot"),
          "Grid" = plotOutput("plot"),
          "Phase Abundance Along Path" = plotOutput("plot")
        )
      )
    }
  })
  output$table <- renderTable(
    switch(input$output_type,
      "Data File" = data_file(crust_out(), x_n = length(crust_out()[[1]]), y_n = length(crust_out()), input$choose_columns, input$choose_rows, input$choose_points, input$assign_label),
      "Grid" = if (!is.null(grid_data_r())) {
        grid_data_r()[[2]]
      },
      "Phase Abundance Along Path" = phase_abundance_r()[[2]],
      "PAM" = matrix(names(PAM_r()[[2]]), , 1, byrow = FALSE, dimnames = list(PAM_r()[[2]], paste("All fields are +", PAM_r()[[3]])))
    ),
    rownames = TRUE
  )
  output$plot <- renderPlot({
    switch(input$output_type,
      "PAM" = draw_PAM_r(),
      "Grid" = draw_Grid_r(),
      "Phase Abundance Along Path" = draw_abundance_r()
    )
  })
  # Dyanmically create output form selections
  output$output_form_selection <- renderUI({
    form_selection <- switch(input$output_type,
      "Data File" = "Data",
      "PAM" = c("Plot", "Legend"),
      c("Plot", "Data")
    )
    radioButtons("output_form", "View", form_selection, inline = TRUE)
  })
  # Fix for plots disappearing.
  # Mod-tag: Need a better way to refresh plots that has better code continuity
  observeEvent(input$output_type, {
    if (!is.null(input$output_form)) {
      inputOption <- ""
      Option1 <- ""
      Option2 <- ""
      unidentifiedPlot <- FALSE
      # setting possible options
      # You need a select input that does not re-render the plot.
      switch(input$output_type,
        "Grid" = {
          inputOption <- "Grid_colours"
          Option1 <- "gray.colors"
          Option2 <- "heat.colors"
        },
        "Phase Abundance Along Path" = {
          inputOption <- "axis"
          Option1 <- "x"
          Option2 <- "y"
        },
        "PAM" = {
          inputOption <- "PAM_labels"
          Option1 <- "Phases"
          Option2 <- "Numbers"
        },
        unidentifiedPlot <- TRUE
        # Template, add your new plot update.
        # , "plotName" = {
        # inputOption <- "inputOption"
        # Option1 <- "Option1"
        # Option2 <- "Option2"
        # }
      )
      # This code should not need to be altered, unless the plot has no Select Input option.
      # will not refresh if the plot has not been added to this fix function.
      if (!unidentifiedPlot) {
        theOption <- Option2
        if (eval(parse(text = paste0("input$", inputOption))) == Option2) {
          theOption <- Option1
          Option1 <- Option2
        }
        eval(parse(text = paste0("updateSelectInput(session, \"", inputOption, "\", selected = \"", theOption, "\")")))
        eval(parse(text = paste0("updateSelectInput(session, \"", inputOption, "\", selected = \"", Option1, "\")")))
      }
    }
  })
  # Dyanmically create output selections
  output$output_selection <- renderUI({
    # Refresh reactive outputs if they exist
    if (exists("crust")) {
      store_r$crust_r <- crust
    }
    if (exists("input_pt")) {
      store_r$input_pt_r <- input_pt
    }
    if (exists("input_bulk")) {
      store_r$input_bulk_r <- input_bulk
    }
    if (exists("major_elements")) {
      store_r$major_elements_r <- major_elements
    }
    if (exists("trace_elements")) {
      store_r$trace_elements_r <- trace_elements
    }
    if (exists("all_elements")) {
      store_r$all_elements_r <- all_elements
    }
    if (is.null(store_r$crust_r)) {
      return("To select ouputs: first run calculation or load previously saved calculation
             *(variables 'crust' and 'major_elements' must be present")
    }
    all_columns <- NULL
    all_phases <- NULL
    for (y_i in 1:length(crust_out())) {
      for (x_i in 1:length(crust_out()[[1]])) {
        all_columns <- union(all_columns, colnames(crust_out()[[y_i]][[x_i]]))
        all_phases <- union(all_phases, rownames(crust_out()[[y_i]][[x_i]]))
      }
    }
    # mod-tag: figure out better way of parsing multiple panels
    conditionalPanel(
      "true",
      conditionalPanel(
        "input.output_type == 'Data File'",
        selectizeInput(
          "choose_columns",
          "Select Columns",
          c(
            "All" = "",
            "Brief",
            "ID",
            "Phase",
            "y_i",
            "x_i",
            "Pressure(kbar)",
            "Temperature(C)",
            "wt%",
            "vol%",
            all_elements,
            "mass",
            "G(J)",
            "V(J/bar)",
            "H(J)",
            "Gruneisen_T",
            "Ks(bar)",
            "Mu(bar)",
            "V0(km/s)",
            "Vp(km/s)",
            "Vs(km/s)",
            "Vp/Vs",
            "Rho(kg/m3)",
            "Cp(J/K)",
            "alpha(1/K)",
            "beta(1/bar)",
            "S(J/K)",
            "N(g)",
            "Cp/Cv"
          ),
          multiple = TRUE
        ),
        # mod-tag should populate all potential columns here
        selectizeInput("choose_rows", "Select System/Phase", c("All" = "", "Reactive subsystem", "Extract subsystem", all_phases), multiple = TRUE),
        textInput("choose_points", "Select Points", "{1;1}"),
        textInput("from_label", "From"),
        textInput("to_label", "To"),
        textInput("label_name", "Label Name"),
        textInput("label_value", "Label Value"),
        actionButton("assign_label", "Assign Label")
      ),
      conditionalPanel(
        "input.output_type == 'Grid'",
        selectInput("Grid_variable", "Variable", c("Pressure", "Temperature", all_columns, "A/CNK", "Custom"), selected = "Pressure", selectize = TRUE),
        conditionalPanel(
          condition = "input.Grid_variable.indexOf('Custom') != -1",
          textInput("Custom_selection", "Custom Selection")
        ),
        conditionalPanel(
          "['Pressure','Temperature','Custom'].indexOf(input.Grid_variable) == -1",
          selectInput("Grid_variable_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
        ),
        # Mod-tag: Make remove phase here. Condition must be that Grid_variable is wt% or vol% or mol%
        # mod-tag: allow incrments setting
        # selectizeInput('Grid_variable_increments', NULL, c("Increments",3:y_n),"Increments") ,
        selectInput("Grid_axes", "Labelled axes", c("bottom", "left", "top", "right"), selected = c("bottom", "left"), selectize = TRUE, multiple = TRUE),
        conditionalPanel(
          condition = "input.Grid_axes.indexOf('bottom') != -1",
          selectInput(
            "Grid_bottom_axis",
            "Bottom Axis",
            c(
              "x_i",
              "Pressure (kbar)",
              "Temperature (C)",
              "Pressure (Mpa)",
              "Temperature (K)",
              colnames(crust[[1]][[1]]),
              "Custom"
            ),
            selected = "x_i",
            selectize = TRUE
          ),
          conditionalPanel(
            "['x_i','Pressure','Temperature'].indexOf(input.Grid_bottom_axis) == -1",
            selectInput("Grid_bottom_axis_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
          ),
          selectizeInput("Grid_bottom_axis_increments", NULL, c("Increments", 3:x_n), "Increments")
        ),
        conditionalPanel(
          condition = "input.Grid_axes.indexOf('left') != -1",
          selectInput(
            "Grid_left_axis",
            "Left Axis",
            c(
              "y_i",
              "Pressure (kbar)",
              "Temperature (C)",
              "Pressure (Mpa)",
              "Temperature (K)",
              colnames(crust[[1]][[1]])
            ),
            selected = "y_i",
            selectize = TRUE
          ),
          conditionalPanel(
            "['y_i','Pressure','Temperature'].indexOf(input.Grid_left_axis) == -1",
            selectInput("Grid_left_axis_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
          ),
          selectizeInput("Grid_left_axis_increments", NULL, c("Increments", 3:y_n), "Increments")
        ),
        conditionalPanel(
          condition = "input.Grid_axes.indexOf('top') != -1",
          selectInput(
            "Grid_top_axis",
            "Top Axis",
            c(
              "x_i",
              "Pressure (kbar)",
              "Temperature (C)",
              "Pressure (Mpa)",
              "Temperature (K)",
              colnames(crust[[1]][[1]])
            ),
            selected = "x_i",
            selectize = TRUE
          ),
          conditionalPanel(
            "['x_i','Pressure','Temperature'].indexOf(input.Grid_top_axis) == -1",
            selectInput("Grid_top_axis_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
          ),
          selectizeInput("Grid_top_axis_increments", NULL, c("Increments", 3:x_n), "Increments")
        ),
        conditionalPanel(
          condition = "input.Grid_axes.indexOf('right') != -1",
          selectInput(
            "Grid_right_axis",
            "Right Axis",
            c(
              "y_i",
              "Pressure (kbar)",
              "Temperature (C)",
              "Pressure (Mpa)",
              "Temperature (K)",
              colnames(crust[[1]][[1]])
            ),
            selected = "y_i",
            selectize = TRUE
          ),
          conditionalPanel(
            "['y_i','Pressure','Temperature'].indexOf(input.Grid_right_axis) == -1",
            selectInput("Grid_right_axis_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
          ),
          selectizeInput("Grid_right_axis_increments", NULL, c("Increments", 3:y_n), "Increments")
        ),
        textInput("remove_values", "Remove Values"),
        selectInput(
          "Grid_colours",
          "Colour Scheme",
          c(
            "gray.colors",
            "heat.colors",
            "terrain.colors",
            "rainbow",
            "topo.colors",
            "RColorBrewer"
          ),
          selected = "gray.colors",
          selectize = TRUE
        ),
        conditionalPanel(
          "input.Grid_colours == 'RColorBrewer'",
          selectInput("Grid_RColorBrewer_colours", "RColorBrewer Pallette", rev(rownames(RColorBrewer::brewer.pal.info)), selected = "YlOrRd", selectize = TRUE),
          textInput("Grid_RColorBrewer_levels", "Number of levels", 10),
          textInput("Grid_RColorBrewer_zlimits", "Z limits", NULL),
          checkboxInput("Grid_RColorBrewer_reverse", "Reverse colours?", value = FALSE)
        ),
        radioButtons("rotation", "Rotation", c("0" = 0, "90" = 1, "180" = 2, "270" = 3), selected = 0, inline = TRUE, width = "200px"),
        checkboxGroupInput("reflection", "Reflection", c("Horizontal", "Vertical"))
      ),
      conditionalPanel(
        "input.output_type == 'Phase Abundance Along Path'",
        selectInput("axis", "Axis", c("x", "y"), selected = "x", selectize = TRUE),
        selectInput("proportion", "Proportion", c("mass", "wt%", "vol%"), selected = "mass", selectize = TRUE),
        conditionalPanel(
          "input.axis == 'x'",
          selectInput("path_y", "Path", 1:length(crust), selected = 1, selectize = TRUE),
          selectInput("start_x", "Start Point", 1:length(crust[[1]]), selected = 1, selectize = TRUE),
          selectInput("end_x", "End Point", 1:length(crust[[1]]), selected = 1, selectize = TRUE)
        ),
        conditionalPanel(
          "input.axis == 'y'",
          selectInput("path_x", "Path", 1:length(crust[[1]]), selected = 1, selectize = TRUE),
          selectInput("start_y", "Start Point", 1:length(crust), selected = 1, selectize = TRUE),
          selectInput("end_y", "End Point", 1:length(crust), selected = 1, selectize = TRUE)
        ),
        uiOutput("select_abundance_phases"),
        selectInput("path_label", "Path Label", c("Point", "Pressure(kbar)", "Temperature(C)"), selected = "Point", selectize = TRUE),
        selectInput("legend", "Legend", c("None", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"), selected = "topright", selectize = TRUE)
      ),
      conditionalPanel(
        "input.output_type == 'PAM'",
        selectInput("PAM_system", "System", c("Reactive Subsystem", "Extract Subsystem", "Full System"), selected = "Reactive Subsystem", selectize = TRUE),
        selectInput("PAM_labels", "Field Labels", c("Phases", "Numbers"), selected = "Phases", selectize = TRUE),
        selectInput("PAM_axes", "Labelled axes", c("bottom", "left", "top", "right"), selected = c("bottom", "left"), selectize = TRUE, multiple = TRUE),
        conditionalPanel(
          condition = "input.PAM_axes.indexOf('bottom') != -1",
          selectInput("PAM_bottom_axis", "Bottom Axis", c("x_i", "Pressure", "Temperature", colnames(crust[[1]][[1]])), selected = "x_i", selectize = TRUE),
          conditionalPanel(
            "['x_i','Pressure','Temperature'].indexOf(input.PAM_bottom_axis) == -1",
            selectInput("PAM_bottom_axis_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
          ),
          selectizeInput("PAM_bottom_axis_increments", NULL, c("Increments", 3:x_n), "Increments")
        ),
        conditionalPanel(
          condition = "input.PAM_axes.indexOf('left') != -1",
          selectInput("PAM_left_axis", "Left Axis", c("y_i", "Pressure", "Temperature", colnames(crust[[1]][[1]])), selected = "y_i", selectize = TRUE),
          conditionalPanel(
            "['y_i','Pressure','Temperature'].indexOf(input.PAM_left_axis) == -1",
            selectInput("PAM_left_axis_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
          ),
          selectizeInput("PAM_left_axis_increments", NULL, c("Increments", 3:y_n), "Increments")
        ),
        conditionalPanel(
          condition = "input.PAM_axes.indexOf('top') != -1",
          selectInput("PAM_top_axis", "Top Axis", c("x_i", "Pressure", "Temperature", colnames(crust[[1]][[1]])), selected = "x_i", selectize = TRUE),
          conditionalPanel(
            "['x_i','Pressure','Temperature'].indexOf(input.PAM_top_axis) == -1",
            selectInput("PAM_top_axis_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
          ),
          selectizeInput("PAM_top_axis_increments", NULL, c("Increments", 3:x_n), "Increments")
        ),
        conditionalPanel(
          condition = "input.PAM_axes.indexOf('right') != -1",
          selectInput("PAM_right_axis", "Right Axis", c("y_i", "Pressure", "Temperature", colnames(crust[[1]][[1]])), selected = "y_i", selectize = TRUE),
          conditionalPanel(
            "['y_i','Pressure','Temperature'].indexOf(input.PAM_right_axis) == -1",
            selectInput("PAM_right_axis_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
          ),
          selectizeInput("PAM_right_axis_increments", NULL, c("Increments", 3:y_n), "Increments")
        ),
        selectInput("PAM_contour", "Contour", c("None", "Pressure", "Temperature", colnames(crust[[1]][[1]])), selected = "None", selectize = TRUE),
        conditionalPanel(
          "['None','Pressure','Temperature'].indexOf(input.PAM_contour) == -1",
          selectInput("PAM_contour_grid_phase", NULL, sort(all_phases), selected = "Bulk_rs", selectize = TRUE)
        ),
        selectizeInput("PAM_contour_increments", NULL, c("Default Increments", "In/Out", "Set Levels", "Set Contours", 3:y_n), "Default Increments"),
        textInput("PAM_compilation", "PAM Compilation"),
        actionButton("create_compilation", "Compile/Refresh Legend"),
        checkboxInput("compile_PAM", "Apply compilation?", value = FALSE)
      ),
      selectInput("file_type", "File type", c(".csv", ".txt", ".ps"), selected = ".csv", selectize = TRUE),
      actionButton("save_data", "Save To File")
      # mod-tag: allow this functionality
      , actionButton("send_gcdkit", "Send To GCDkit")
    )
  })
  # if working_file is not blank on first opening, load file
  observe({
    if (!exists("first_load")) {
      first_load <<- TRUE
    }
    if (first_load) {
      projects_directory <- input$projects_directory
      working_file <- input$working_file
      # error handling
      reactive_message$data <- error_handling(working_file, projects_directory)
      if (reactive_message$data == "error handling passed") {
        # load
        if (file.exists(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))) {
          reactive_message$data <- paste0(on_load())
        } else {
          reactive_message$data <- paste0("No input file found at ", paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
        }
        # Load workspace if it exists (previous calculation results)
        if (file.exists(paste0(projects_directory, "/", working_file, "/", working_file, ".RData"))) {
          load(paste0(projects_directory, "/", working_file, "/", working_file, ".RData"), envir = .GlobalEnv)
        }
        # Refresh reactive outputs if they exist
        if (exists("crust")) {
          store_r$crust_r <- crust
        }
        if (exists("input_pt")) {
          store_r$input_pt_r <- input_pt
        }
        if (exists("input_bulk")) {
          store_r$input_bulk_r <- input_bulk
        }
        if (exists("major_elements")) {
          store_r$major_elements_r <- major_elements
        }
        if (exists("trace_elements")) {
          store_r$trace_elements_r <- trace_elements
        }
      } else {
        reactive_message$data
      }
    }
    first_load <<- FALSE
  })
})
