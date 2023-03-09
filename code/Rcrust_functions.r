# Repository for Rcrust functions
crust_to_gcdkit <- function(crust, choose_columns = NULL, choose_rows = NULL, choose_points = "All", GCDkitGUI = FALSE, source.first = TRUE, source.plugins = TRUE) {
  # crust_to_gcdkit(crust,choose_points="{1;1}")
  # sends crust object to gcdkit
  # note this function has to assign globally so that GCDkit can see it
  data_crust <- as.data.frame(data_file(crust, x_n, y_n, choose_columns, choose_rows, choose_points))
  rownames(data_crust) <- data_crust[, 1]
  require(GCDkitDevelop)
  # library(GCDkit)
  stopApp(returnValue = invisible())
  # source(paste(gcdx.dir,"/GCDkit.r",sep=""))
  # setwd(gcdx.dir)
  accessVar(data_crust, source.first = source.first, source.plugins = source.plugins, poke.data = FALSE, GUI = GCDkitGUI)
  print("To launch the Rcrust GUI type Rcrust() into the console and press enter")
  # Provide again the launch with GUI function in order to get back
  Rcrust <<- function() {
    # If working directory is x\Projects\y then set to x\code
    if (length(grep("Rcrust/Projects/", getwd())) == 1) {
      setwd(paste0(strsplit(getwd(), split = "Projects")[[1]][1], "code"))
    }
    runApp()
  }
}

assign_label <- function(crust, from_label, to_label, label_name, label_value) {
  # adds column to crust and populates it with value over range
  # example
  # label_name<-"Protolith"
  # from_label<-"{1;1}"
  # to_label<-"{1;1}"
  # label_value<-"Sample1"
  a <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", from_label)), split = ";"))
  b <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", to_label)), split = ";"))
  for (x_i in a[1]:b[1]) {
    for (y_i in a[2]:b[2]) {
      if (!any(colnames(crust[[y_i]][[x_i]]) == label_name)) {
        lab_column <- matrix(label_value, nrow(crust[[y_i]][[x_i]]), 1)
        colnames(lab_column) <- label_name
        crust[[y_i]][[x_i]] <- cbind(crust[[y_i]][[x_i]], lab_column)
      } else {
        crust[[y_i]][[x_i]][, lab_column] <- label_value
      }
    }
  }
  return(crust)
}
# Modified from GCDmodel, written by Jean-francois Moyen 2019
.sanitize2 <- function(obj, normalize.TO = 0) {
  if (class(obj)[1] == "matrix") {
    nn <- colnames(obj)
    obj <- as.vector(obj)
    names(obj) <- nn
  }
  if (class(obj)[1] == "data.frame") {
    nn <- colnames(obj)
    obj <- as.numeric(obj)
    names(obj) <- nn
  }
  mis <- is.na(obj)
  obj <- obj[!mis]
  if (normalize.TO > 0) {
    obj <- obj / sum(obj) * normalize.TO
  }
  return(obj)
}
# Modified from GCDmodel, written by Jean-francois Moyen 2019
BatchPM2 <- function(kd, c0, pm, cmins = matrix(), min.props, melt.arg = list(), dont = character(0)) {
  c0 <- .sanitize2(c0)
  min.props <- .sanitize2(min.props, normalize.TO = 1)
  which.elems <- intersect(names(c0), colnames(kd))
  which.mins <- names(min.props)
  missing <- setdiff(names(min.props), rownames(kd))
  if (length(missing) != 0) {
    stop(paste(
      "Missing Kd value for mineral(s)", missing,
      "\n"
    ))
  }
  c0 <- c0[which.elems]
  kd <- kd[which.mins, which.elems, drop = F]
  DD <- min.props %*% kd
  FF <- pm / 100
  cL <- c0 / (DD + FF * (1 - DD))
  ee <- sapply(names(min.props), function(i) {
    z <- cL[1, which.elems] * kd[i, which.elems]
    return(z)
  }, simplify = TRUE)
  cmins <- t(ee)
  cmins <- cmins[which.mins, which.elems, drop = F]
  cS <- cL[1, which.elems] * DD[1, which.elems]
  invisible(list(
    c0 = c0,
    cL = .sanitize2(cL),
    cmins = cmins,
    cS = cS,
    min.props = min.props,
    FF = pm,
    kd = kd,
    DD = .sanitize2(DD)
  ))
}

# exists_and_true
exists_and_true <- function(x) {
  chk <- try(x, silent = TRUE)
  if (class(chk) == "try-error") {
    return(FALSE)
  } else {
    return(x)
  }
}
#############################################
#
# Ancillary function to merge several lines (by mass)
#
############################################
.wtd.add <- function(thelines, prop = "mass", avname = "Averaged") {
  if (nrow(thelines) == 1) {
    foo <- thelines
  } else {
    # Two sets of cols
    if (exists_and_true(calc_mol)) {
      extensive.cn <- c("wt%", "vol%", "mol%", "mol") # Extensive properties (mass dependant) -- add the others if required
    } else {
      extensive.cn <- c("wt%")
    }
    intensive.cn <- setdiff(colnames(thelines), c(prop, extensive.cn))
    foo <- matrix(rep(0, length(colnames(thelines))), nrow = 1)
    colnames(foo) <- colnames(thelines)
    # replace NA with 0
    thelines[is.na(thelines)] <- 0
    # Intensive
    foo[, intensive.cn] <- thelines[, prop] %*% thelines[, intensive.cn, drop = F] / sum(thelines[, prop])
    # Extensive
    foo[, prop] <- sum(thelines[, prop])
    if (!length(intersect(extensive.cn, thelines)) == 0) {
      foo[, extensive.cn] <- colSums(thelines[, extensive.cn])
    }
  }
  rownames(foo) <- avname
  return(foo)
}
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
# fix-tag remove working directory and projects directory and replace with gsub argument to auto locate relative to opened file
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
phase_abundance <- function(crust, axis, path = 1, p_a = 1, p_b = p_a, path_label = "Point", input_pt = NULL) {
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
      chk_phase <- try(crust[[eval(parse(text = pnt_y))]][[eval(parse(text = pnt_x))]][ph, "mass"], silent = TRUE)
      if (class(chk_phase) == "try-error") {
        abundace_phase <- c(abundace_phase, 0)
      } else {
        abundace_phase <- c(abundace_phase, crust[[eval(parse(text = pnt_y))]][[eval(parse(text = pnt_x))]][ph, "mass"])
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
    validate(need(file.exists(paste0(sub("/code", "/Projects", getwd()), "/Compile/", PAM_compilation, " compilation legend.txt")), paste0(PAM_compilation, " compilation legend.txt not found in ", sub("/code", "/Projects", getwd()), "/Compile/", "\nPlease compile legend first")))
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
data_file <- function(crust, x_n = length(crust[[1]]), y_n = length(crust), choose_columns = NULL, choose_rows = NULL, choose_points = "All") {
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
    choose_columns <- union(choose_columns[-which(choose_columns == "Brief")], c("ID", "Phase", "y_i", "x_i", "Pressure(kbar)", "Temperature(C)", "wt%", comps, "mass"))
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
#############################################
#
# Expression evaluator
#
############################################
# function-def: eval_expr<-function(expr,calc_phases=calc_phases,crust=crust)
# Evaluate expression a given calc_phases and crust
# () are for solution models, {} are for function terms and bodmas
# evaluate any { and word before it up until }
eval_expr <- function(expr, calc_phases = calc_phases, crust = crust) {
  # wrap outside for evaluation
  a <- paste0("{", expr, "}")
  while (!unlist(gregexpr("[{]", a))[1] == -1) {
    letters <- gregexpr("[a-z]", a)
    left_bracs <- unlist(gregexpr("[{]", a))
    names(left_bracs) <- rep("left", length(left_bracs))
    right_bracs <- unlist(gregexpr("[}]", a))
    names(right_bracs) <- rep("right", length(right_bracs))
    bracs <- sort(c(left_bracs, right_bracs))
    # Find inner brackets
    i <- 1
    while (!(names(bracs)[i] == "left" & names(bracs)[i + 1] == "right")) {
      i <- i + 1
    }
    # evaluate inner brackets
    # If require function call find name
    j <- 1
    while ((bracs[i] - j) %in% letters[[1]]) {
      j <- j + 1
    }
    # if j==1 dont require function call, else j+1 is first letter of function name
    if (!j == 1) {
      funct_name <- substr(a, bracs[i] - j + 1, bracs[i] - 1)
      # apply function
      if (funct_name == "ph") {
        # arguments of the form ph{phase;unit;x_i;y_i} where unit can be any column name in calc_phases and x_i and y_i are the current point by default
        ph_args <- strsplit(substr(a, bracs[i] + 1, bracs[i + 1] - 1), "[;,]")[[1]]
        if (length(ph_args) == 2) {
          # current_variable
          chk_var <- try(eval(parse(text = "calc_phases[ph_args[1],ph_args[2]]")), silent = TRUE)
          if (class(chk_var) == "try-error") {
            out <- 0
          } else {
            out <- chk_var
          }
        } else {
          # past_variable
          # default look in rs
          if (length(grep("_rs", ph_args[1])) == 0 & length(grep("_es", ph_args[1])) == 0) {
            ph_args[1] <- paste0(ph_args[1], "_rs")
          }
          chk_var <- try(eval(parse(text = "crust[[eval(parse(text=ph_args[[4]]))]][[eval(parse(text=ph_args[[3]]))]][ph_args[1],ph_args[2]]")), silent = TRUE)
          if (class(chk_var) == "try-error") {
            out <- 0
          } else {
            out <- chk_var
          }
        }
      }
      if (funct_name == "delta") {
        skip_delta <- FALSE
        # arguments of the form delta{ph;x_#;y_#;unit} where unit can be wt% or mass
        delta_args <- strsplit(substr(a, bracs[i] + 1, bracs[i + 1] - 1), "[;,]")[[1]]
        delta_phs <- delta_args[1]
        if (substring(delta_args[2], 1, 8) == "prev_ext") {
          find_ph <- substring(delta_args[2], 10)
          pnt <- "find"
          l <- 1
          while (pnt == "find") {
            # find earliest of phase extraction or phase absent point, if neither show warning and take pnt=1
            # find previous increase in _es_cumul
            chk_pnt <- try(crust[[y_i]][[x_i - l]][paste0(find_ph, "_es"), "mass"], silent = TRUE)
            if (class(chk_pnt) == "try-error") {
              chk_pnt <- try(crust[[y_i]][[x_i - l]][paste0(find_ph, "_rs"), "mass"], silent = TRUE)
              if (class(chk_pnt) == "try-error") {
                pnt <- x_i - l
                cat(paste0("Delta calculated to previous phase absent point at x_i = ", pnt, "\n"))
              }
            } else {
              pnt <- x_i - l
              cat(paste0("Delta calculated to previous extract at x_i = ", pnt, "\n"))
            }
            if ((x_i - l) == 1 & pnt == "find") {
              pnt <- x_i - l
              cat(paste0("Delta calculated to first point at x_i = ", pnt, "\n"))
            }
            l <- l + 1
          }
          delta_index_x <- pnt
        } else {
          delta_index_x <- eval(parse(text = delta_args[2]))
        }
        if (!(delta_index_x <= x_n & delta_index_x >= 1)) {
          cat("Warning delta_index_x not <x_n and >=1, skipping this delta calculation\n")
          skip_delta <- TRUE
        }
        delta_index_y <- eval(parse(text = delta_args[3]))
        if (!(delta_index_y <= y_n & delta_index_y >= 1)) {
          cat("Warning delta_index_y not <y_n and >=1, skipping this delta calculation\n")
          skip_delta <- TRUE
        }
        delta_unit <- delta_args[4]
        if (!(delta_unit == "mass" | delta_unit == "wt%")) {
          cat("Error delta only currently accepted as mass or wt%\n")
          stop()
        }
        if (!skip_delta) {
          current_mode <- 0
          previous_mode <- 0
          # Check for plus sign (e.g. aluminosilicate given as sill+ky+and)
          delta_phases <- strsplit(delta_args[1], "+", fixed = TRUE)[[1]]
          for (each_ph in delta_phases) {
            # current_mode
            chk_var <- try(eval(parse(text = "calc_phases[each_ph,delta_unit]")), silent = TRUE)
            if (!class(chk_var) == "try-error") {
              current_mode <- current_mode + chk_var
            }
            # previous_mode
            delta_phase_rs <- paste(each_ph, "_rs", sep = "")
            chk_var <- try(eval(parse(text = "crust[[delta_index_y]][[delta_index_x]][delta_phase_rs,delta_unit]")), silent = TRUE)
            if (!class(chk_var) == "try-error") {
              previous_mode <- previous_mode + chk_var
            }
          }
          # delta
          if (current_mode > previous_mode) {
            out <- current_mode - previous_mode
          } else {
            out <- 0
          }
        } else {
          out <- 0
        }
      }
      if (funct_name == "retain") {
        # c# If Retention Mode - Calculate mass of retention phases to extract
        # c# Extract till retention amount of retention phases is left
        # arguments of the form retain{amount;unit;ph} where ph can be omitted to take on the current ph
        retain_args <- strsplit(substr(a, bracs[i] + 1, bracs[i + 1] - 1), "[;,]")[[1]]
        ret <- as.numeric(retain_args[1])
        retention_unit <- retain_args[2]
        ph <- names(expr)
        if (!is.na(retain_args[3])) {
          ph <- retain_args[3]
        }
        bulk_no <- which(rownames(calc_phases) == "Bulk_rs")
        ph_no <- which(rownames(calc_phases) == ph)
        if (retention_unit == "mass") {
          system_less_ret <- calc_phases[-c(bulk_no, ph_no), "mass"]
          if (calc_phases[ph, "mass"] <= ret) {
            out <- 0
          } else {
            out <- calc_phases[ph, "mass"] - ret
          }
        }
        if (retention_unit == "vol%") {
          system_less_ret <- calc_phases[-c(bulk_no, ph_no), "vol%"]
          if (calc_phases[ph, "vol%"] <= ret) {
            out <- 0
          } else {
            new_ret_vol <- (ret * (sum(system_less_ret))) / (100 - ret)
            new_ret_mass <- new_ret_vol * calc_phases[ph, "mass"] / calc_phases[ph, "vol%"]
            out <- calc_phases[ph, "mass"] - new_ret_mass
          }
        }
        if (retention_unit == "wt%") {
          system_less_ret <- calc_phases[-c(bulk_no, ph_no), "wt%"]
          if (calc_phases[ph, "wt%"] <= ret) {
            out <- 0
          } else {
            new_ret_wt <- (ret * (sum(system_less_ret))) / (100 - ret)
            new_ret_mass <- new_ret_wt * calc_phases[ph, "mass"] / calc_phases[ph, "wt%"]
            out <- calc_phases[ph, "mass"] - new_ret_mass
          }
        }
      }
      if (funct_name == "return") {
        # arguments of the form return{phase;amount} where amount can be % or mass
        return_args <- strsplit(substr(a, bracs[i] + 1, bracs[i + 1] - 1), "[;,]")[[1]]
        # check if phase is in extract cumul subsystem
        chk_var <- try(eval(parse(text = paste0("cumul_extract_pnt[\"", return_args[1], "_es_cumul\",\"mass\"]"))), silent = TRUE)
        if (class(chk_var) != "try-error") {
          if (!is.null(chk_var)) {
            percentage <- FALSE
            # if percentage tag to calculate
            if (length(grep("%", return_args[2])) != 0) {
              percentage <- TRUE
              return_args[2] <- gsub("%", "", return_args[2])
            }
            # evaluate for number
            chk_num <- try(eval(parse(text = return_args[2])), silent = TRUE)
            if (class(chk_num) == "try-error") {
              cat("Error phase addition could not evaluate isolated function correctly")
              stop()
            }
            # calculate
            if (percentage) {
              take <- chk_num / 100 * chk_var
              leave <- (100 - chk_num) / 100 * chk_var
            } else {
              take <- chk_num
              leave <- as.numeric(chk_var) - chk_num
            }
            if (leave < 0) {
              take <- chk_var
              leave <- 0
            }
            # transfer
            cumul_extract_pnt[paste0(return_args[1], "_es_cumul"), "mass"] <<- leave
            # leave calc
            a <- paste(c(cumul_extract_pnt[paste0(return_args[1], "_es_cumul"), comps], take), collapse = ",")
            break
          } else {
            cat(paste0("Phase ", return_args[1], " not found in extract cumul\n"))
            a <- paste(rep(0, length(comps) + 1), collapse = ",")
            break
            out <- ""
          }
        } else {
          cat(paste0("Phase ", return_args[1], " not found in extract cumul\n"))
          a <- paste(rep(0, length(comps) + 1), collapse = ",")
          break
          out <- ""
        }
      }
      # replace inner brackets with function output
      a <- paste(substr(a, 1, bracs[i] - j), out, substr(a, bracs[i + 1] + 1, nchar(a)), sep = "")
    } else {
      # evaluate inner brackets for bodmas output
      out <- eval(parse(text = substr(a, bracs[i] + 1, bracs[i + 1] - 1)))
      a <- paste(substr(a, 1, bracs[i] - j), out, substr(a, bracs[i + 1] + 1, nchar(a)), sep = "")
    }
  }
  return(a)
}
# rename kf
renameFsp <- function(all_elements, calc_phases) {
  if (length(intersect(toupper(all_elements), c("CAO", "K2O"))) == 2) {
    for (ph in grep("Fsp", rownames(calc_phases))) {
      CaO_pos <- which(toupper(names(calc_phases[ph, ])) == "CAO")
      K2O_pos <- which(toupper(names(calc_phases[ph, ])) == "K2O")
      if (calc_phases[ph, K2O_pos] <= 0) {
        calc_phases[ph, K2O_pos] <- 0.0001
      }
      if (calc_phases[ph, CaO_pos] / calc_phases[ph, K2O_pos] > 1) {
        rownames(calc_phases)[ph] <- "Pl"
      } else {
        rownames(calc_phases)[ph] <- "Kf"
      }
    }
  }
  return(calc_phases)
}
# rename phases as specified by the user
renamePhases <- function(phases_to_rename, calc_phases) {
  if (!phases_to_rename[1] == "") {
    # Remove duplicate numbering in calc_phases
    phases <- strsplit(rownames(calc_phases), "_")
    phases <- sapply(phases, "[[", 1)
    rownames(calc_phases) <- phases
    # Replace the Bulk_rs label as we'll need it later
    rownames(calc_phases)[which(rownames(calc_phases) == "Bulk")] <- "Bulk_rs"
    for (i in 1:length(phases_to_rename)) {
      rename_inputs <- unlist(strsplit(names(phases_to_rename[i]), "~"))
      # if calculated phase matches phase name after removing underscore
      for (j in which(phases == rename_inputs[1])) {
        # If the phase meets the criteria
        if (eval_expr(phases_to_rename[i], calc_phases[j, , drop = FALSE])) {
          # rename the phase
          rownames(calc_phases)[j] <- rename_inputs[2]
        }
      }
    }
  }
  return(calc_phases)
}
# normalise totals in calc_phases
# mass, wt%, vol%, mol%
normTotals <- function(calc_phases, c0, roundDigits = 5) {
  br <- which(rownames(calc_phases) == "Bulk_rs") - 1
  # mass
  if (any(colnames(calc_phases) == "mass")) {
    calc_phases[1:br, "mass"] <- calc_phases[1:br, "mass"] / sum(calc_phases[1:br, "mass"]) * c0["mass"]
    calc_phases["Bulk_rs", "mass"] <- c0["mass"]
    calc_phases[1:br, "mass"] <- round(calc_phases[1:br, "mass"], roundDigits)
  }
  # wt%
  if (any(colnames(calc_phases) == "wt%")) {
    calc_phases[1:br, "wt%"] <- calc_phases[1:br, "mass"] / sum(calc_phases[1:br, "mass"]) * 100
    calc_phases["Bulk_rs", "wt%"] <- 100
    calc_phases[1:br, "wt%"] <- round(calc_phases[1:br, "wt%"], roundDigits)
  }
  # vol%
  if (any(colnames(calc_phases) == "vol%")) {
    volume <- calc_phases[1:br, "mass"] / calc_phases[1:br, "Density(kg/m3)"]
    calc_phases[1:br, "vol%"] <- volume / sum(volume) * 100
    calc_phases["Bulk_rs", "vol%"] <- 100
    calc_phases[1:br, "vol%"] <- round(calc_phases[1:br, "vol%"], roundDigits)
  }
  # mol%
  if (any(colnames(calc_phases) == "mol%")) {
    calc_phases[1:br, "mol%"] <- calc_phases[1:br, "mol"] / sum(calc_phases[1:br, "mol"]) * 100
    calc_phases["Bulk_rs", "mol%"] <- 100
    calc_phases[1:br, "mol%"] <- round(calc_phases[1:br, "mol%"], roundDigits)
  }
  return(calc_phases)
}
# General function to return the molar A/CNK ratio of comp.
# Function accepts a vector object that must be named as component oxides (including Al2O3, CaO, Na2O, K2O), wt.% units
calcACNK <- function(comp) {
  mwAl2O3 <- 101.96128
  mwCaO <- 56.0774
  mwNa2O <- 61.97894
  mwK2O <- 94.19600
  # Resolving case-sensitive components.
  Al <- comp[[which(toupper(names(comp)) == "AL2O3")]]
  Ca <- comp[[which(toupper(names(comp)) == "CAO")]]
  Na <- comp[[which(toupper(names(comp)) == "NA2O")]]
  K <- comp[[which(toupper(names(comp)) == "K2O")]]
  ACNK <- (Al / mwAl2O3) / ((Ca / mwCaO) + (Na / mwNa2O) + (K / mwK2O))
  return(ACNK)
}
