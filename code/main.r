###############################
## Rcrust (main.r)
###############################
# Get dependencies function
source("dependencies.r")
# Handle Default Configurations -----------------------------------------------

if (!exists("project_name")) {
  project_name <- working_file
}
# .First ----------------------------------------------------------------------
.First <- function(run_gui = TRUE) {
  first_load <<- FALSE
  cat("Loading Dependencies\n")
  load_dependencies()
  # Launch GUI
  Rcrust <<- function() {
    # Load dependencies
    cat("Loading Dependencies\n")
    load_dependencies() # TODO removeme?
    # If working directory is x\Projects\y then set to x\code
    if (length(grep("Rcrust/Projects/", getwd())) == 1) {
      setwd(paste0(strsplit(getwd(), split = "Projects")[[1]][1], "code"))
    }
    if (run_gui == TRUE) {
      runApp(launch.browser = TRUE) # Runs server.R and ui.R
    }
  }
  # Launch without GUI
  manual_load <<- function(working_file, projects_directory = paste0(substring(getwd(), 1, nchar(getwd()) - 4), "Projects")) {
    source(paste0(projects_directory, "/", project_name, "/Inputs/", working_file, ".txt"))
    source("main.r")
  }
  # If working directory is x\Projects\y then set to x\code
  if (length(grep("Rcrust/Projects/", getwd())) == 1) {
    setwd(paste0(strsplit(getwd(), split = "Projects")[[1]][1], "code"))
    first_load <<- TRUE
  }
  if (run_gui == TRUE) {
    runApp(launch.browser = TRUE) # Runs server.R and ui.R
  }
}
# End .First ------------------------------------------------------------------
# CALEB TODO Why is this section here? Repeat of .First segment? --------------
# Launch without GUI
manual_load <<- function(working_file, projects_directory = paste0(substring(getwd(), 1, nchar(getwd()) - 4), "Projects")) {
  source(paste0(projects_directory, "/", working_file, "/Inputs/", working_file, ".txt"))
  source("main.r")
}
# END TODO --------------------------------------------------------------------

#############################################
#
# Additional Settings (Advanced users)
#
#############################################
#############################################
# Errors
run_errors <<- NULL
exit_calc <<- FALSE
pause_on_error <<- FALSE
# Calculation settings
request_readline <- FALSE
stop_errors <- NULL
# Determine calculation mode   #normal,parallel (still to come)
calc_mode <- "normal"
silent_calc <- TRUE
cumulate_on_rs <- TRUE
# Override number of cores to use in parallel computation, if enabled (still to come)
# Core count is limited to number of cores on machine
# override_cores<-TRUE
# override_core_count<-2
# TODO Caleb perhaps set number cores to 1 if normal mode?
# PT settings
pt_def <- "input"
PT_restrictions <- NULL
# PT_restrictions<-c("Pressure_>_0","Temperature_>_0")
# PT_restrictions<-c("Pressure_>_2.5","Pressure_<_20","Temperature_>_600","Temperature_<_1100")
# Reaction buffering
reaction_buffering <<- FALSE
reaction_buffer_steps <<- 1
start_in_reaction <<- FALSE
merge_duplicates <<- FALSE
# Sec to clock function
sec_to_clock <- function(total_sec) {
  hr <- floor(total_sec / 3600)
  mn <- floor((total_sec - hr * 3600) / 60)
  if (mn < 10) {
    mn_0 <- "0"
  } else {
    mn_0 <- NULL
  }
  sec <- total_sec - hr * 3600 - mn * 60
  if (sec < 10) {
    sec_0 <- "0"
  } else {
    sec_0 <- NULL
  }
  clock <- paste(hr, ":", mn_0, mn, ":", sec_0, sec, sep = "")
  return(clock)
}
#############################################
#
# Set defaults
#
#############################################
#############################################
if (class(try(reequilibrate_steps, silent = TRUE)) == "try-error") {
  reequilibrate_steps <- TRUE
}
if (class(try(calculate_activities, silent = TRUE)) == "try-error") {
  calculate_activities <- FALSE
}
# Sean-tag:
if (class(try(component_packet, silent = TRUE)) == "try-error") {
  component_packet <- FALSE
}
if (class(try(calc_choice, silent = TRUE)) == "try-error") {
  calc_choice <- "read.meemum" # options = lars.wrap,read.meemum
}
#############################################
#
# Initialize datastructures and inputs
#
#############################################
# Saves the operating system to global variable.
# Values: 'Windows', TODO fill in rest here
assign("OperatingSystem", Sys.info()["sysname"], envir = .GlobalEnv)
cat("Running on Operating System: ", OperatingSystem, "\n")

# cat("Loading Dependencies\n")
# load_dependencies()

#-- Workaround to run from cmd line
args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L) {
  # If no args, run interface mode
  if (length(args) == 0) {
    .First() # Set up everything, and run GUI (runApp call)
  } else {
    .First(run_gui = FALSE) # Set up everything, but skip runApp call
    manual_load(args[1])
  }
} else { # Run main.r code
  # manual_load, runApp both sources main.r, so this prevents double runs
  #--
  cat("\n")
  if (exists("working_file")) {
    cat(paste0("Starting up ", working_file, " for project ", project_name, ".\n"))
  }
  source("run.Rcrust.R") # The function doing the job
  chk_valid <- try({
    # TODO CALEB Linux is case-sensitive. Any reason why some files named with R and some with r?
    source("init_bulk.r") # Bulk composition
    source("init_wrapper.r") # Initialise wrapper function
    source("meemum_connect.r") # Initialise run and read meemum functions
    source("init_meem.r") # Create meemum build file
    source("init_pt.r") # PT conditions
    source("init_ph_add.R") # phase addition
    source("init_ph_extr.r") # phase extraction
    source("init_dependence.R") # Determine dependence relations
    if (calculate_activities) {
      source("init_activities.r") # Gibbs enthalpies for activity calculations
    }
    # Sean-tag:
    source("apSaturation.r") # Saturation calculations involving apatite
    source("Rcrust_functions.r")
  })
  # Set number of cores for parallel computation (TODO Caleb)
  # core_count<-detectCores()
  # if(override_cores==TRUE) {
  #   # Trying to assign more cores than available. Warn and set max
  #   if (override_core_count > core_count) {
  #     cat("WARNING: Attempting to assign more cores (", override_core_count,
  #     ") than available (", core_count, ").\n")
  #     cat("Assigning maximum number of cores possible.")
  #   } else {
  #     core_count = override_core_count
  #     cat("Core count overridden to (", core_count, ")\n")
  #   }
  # }
  # Input validation
  input_valid <- TRUE
  if (class(chk_valid) == "try-error") {
    cat("Initiation failed\n")
    input_valid <- FALSE
  }
  if (input_valid) {
    # Clear meemum_outputs if they exist
    if (export_meemum_output) {
      meemum_output_path <- paste(projects_directory, working_file, "Outputs", paste0(working_file, "_output_meemum.txt"), sep = "/")
      close(file(meemum_output_path, open = "w"))
    }
    if (request_readline) {
      cat("Initiation succesful:\n   Please read the above lines and make sure this is what you wanted.\n")
      rd <- readline(prompt = "Choose \"n\" to abort or press [enter] to continue\n")
    } else {
      cat("Initiation succesful:\n   Computation beginning\n")
      rd <- "pass"
    }
    if (!rd == "n") {
      # Normal calc by tiers
      # get values and pass them to run.Rcrust
      # Create calculation structures
      # object-def: crust[[y_i]][[x_i]]
      crust <- rep(list(rep(list(NULL), x_n)), y_n)
      cumul_extract <- rep(list(rep(list(NULL), x_n)), y_n)
      cumul_add <- rep(list(rep(list(NULL), x_n)), y_n)
      calculation_matrix <- matrix(0, y_n, x_n)
      strt <<- proc.time()
      # Refresh input values for groups (j) and their members (i)
      for (j in 1:length(calc_order)) {
        for (i in 1:length(calc_order[[j]])) {
          # get and assign variables for each point in this calc group
          x_i <- calc_order[[j]][[i]]$x_i
          y_i <- calc_order[[j]][[i]]$y_i
          press <- input_pt[[y_i]][[x_i]][1]
          temp <- input_pt[[y_i]][[x_i]][2]
          # Evaluate c0 till numeric
          # Get c0
          c0 <- input_bulk[[y_i]][[x_i]]
          # replace rs,es,as,fs tuples with values
          for (k in 1:length(c0)) {
            c0_k <- c0[k]
            left_brac <- gregexpr("\\{", c0_k)[[1]]
            right_brac <- gregexpr("\\}", c0_k)[[1]]
            if (left_brac[1] > 0) {
              subsystem <- NULL
              for (kk in 1:length(left_brac)) {
                subsystem <- c(subsystem, substr(c0[k], (left_brac[kk] - 2), left_brac[kk] - 1))
              }
              tuples <- NULL
              for (tup_i in 1:length(left_brac)) {
                tuples <- c(tuples, substr(c0_k, left_brac[tup_i], right_brac[tup_i]))
              }
              split_tuples <- strsplit(gsub("\\{", "", gsub("\\}", "", tuples)), split = ";")
              # Grab dependent point if cumulating extracts - caution this only works if dependent entirely on one point
              if (k == 1 & ph_extr & cumulate_on_rs) {
                dep_pnt <- c(eval(parse(text = split_tuples[[1]][1])), eval(parse(text = split_tuples[[1]][2])))
                crust_rows <- rownames(crust[[dep_pnt[2]]][[dep_pnt[1]]])
                parse_cumul_extract <- crust[[dep_pnt[2]]][[dep_pnt[1]]][crust_rows[grep("_es_cumul", crust_rows)], , drop = FALSE]
                if (!length(parse_cumul_extract) == 0) {
                  cumul_extract[[y_i]][[x_i]] <- parse_cumul_extract
                }
              }
              # Grab dependent point if cumulating additions - caution this only works if dependent entirely on one point
              if (k == 1 & ph_add & cumulate_on_rs) {
                dep_pnt <- c(eval(parse(text = split_tuples[[1]][1])), eval(parse(text = split_tuples[[1]][2])))
                crust_rows <- rownames(crust[[dep_pnt[2]]][[dep_pnt[1]]])
                parse_cumul_addition <- crust[[dep_pnt[2]]][[dep_pnt[1]]][crust_rows[grep("_as_cumul", crust_rows)], , drop = FALSE]
                if (!length(parse_cumul_addition) == 0) {
                  cumul_add[[y_i]][[x_i]] <- parse_cumul_addition
                }
              }
              # substitute in value
              c0_k <- gsub("rs", "", gsub("es", "", gsub("as", "", gsub("fs", "", c0_k))))
              c0_k <- gsub("\\{", "", gsub("\\}", "", c0_k))
              for (split_i in 1:length(split_tuples)) {
                outval <- try(crust[[eval(parse(text = split_tuples[[split_i]][2]))]][[eval(parse(text = split_tuples[[split_i]][1]))]][paste0("Bulk_", subsystem[split_i]), names(c0_k)], silent = TRUE)
                if (class(outval) == "try-error") {
                  outval <- 0
                }
                if (is.null(outval)) {
                  c0_k_try <- 0
                } else {
                  c0_k_try <- try(sub(paste0(unlist(split_tuples[split_i]), collapse = ";"), outval, c0_k))
                }
                if (class(c0_k_try) != "try-error") {
                  c0_k <- c0_k_try
                }
              }
            }
            c0[k] <- eval(parse(text = c0_k))
          }
          # Make numeric
          nam <- names(c0)
          c0 <- as.numeric(c0)
          names(c0) <- nam
          # Save to input_bulk
          input_bulk[[y_i]][[x_i]] <- c0
        }
        if (calc_mode == "normal") {
          # run Rcrust in singular for each point in the calc group
          for (i in 1:length(calc_order[[j]])) {
            x_i <- calc_order[[j]][[i]]$x_i
            y_i <- calc_order[[j]][[i]]$y_i
            if (silent_calc) {
              # Calculate times
              pull_time <- proc.time() - strt
              total_sec <- round(pull_time[3])
              run_time <- sec_to_clock(total_sec)
              cat("Computing Point", "x_i=", x_i, " ; y_i=", y_i, "... Simulation ", round(sum(calculation_matrix != 0) / (nrow(calculation_matrix) * ncol(calculation_matrix)) * 100, 2), "% complete\n")
              cat("Total run time:", run_time, "\n")
              flush.console()
            }
            # Check if point should be calculated
            # object-def: 0=remaining,1=calculated,2=aborted
            chk <- try(PT_restrictions, silent = TRUE)
            if (class(chk) == "try-error") {
              PT_restrictions <- NULL
            }
            if (!is.null(PT_restrictions[1])) {
              PT_split <- strsplit(PT_restrictions, split = "_")
              for (i in 1:length(PT_restrictions)) {
                chk <- paste(input_pt[[y_i]][[x_i]][, PT_split[[i]][1]], PT_split[[i]][2], PT_split[[i]][3], collapse = "")
                if (!eval(parse(text = chk))) {
                  calculation_matrix[y_i, x_i] <- 2
                  crust[[y_i]][[x_i]] <- matrix(c0, 1, )
                  rownames(crust[[y_i]][[x_i]]) <- "Bulk_rs"
                  colnames(crust[[y_i]][[x_i]]) <- names(c0)
                }
              }
            }
            if (calculation_matrix[y_i, x_i] == 0) {
              crust[[y_i]][[x_i]] <- run.Rcrust(
                comps = comps, c0 = input_bulk[[y_i]][[x_i]],
                press = input_pt[[y_i]][[x_i]][1],
                temp = input_pt[[y_i]][[x_i]][2],
                ph_extr_pnt = input_ph_extr[[y_i]][[x_i]],
                cumul_extract_pnt = cumul_extract[[y_i]][[x_i]],
                ph_add_pnt = input_ph_add[[y_i]][[x_i]],
                cumul_add_pnt = cumul_add[[y_i]][[x_i]]
              )
              calculation_matrix[y_i, x_i] <- 1
            }
          }
        } else {
          # Troy
          library(snow)
          # Run rcrust in parallel
          for (i in 1:length(calc_order[[j]])) {
            x_i <- calc_order[[j]][[i]]$x_i
            y_i <- calc_order[[j]][[i]]$y_i
            # cl<-makeCluster(3,type="SOCK")
            # square<-function(x){x^2}
            # clusterApply(cl,c(1:3),square)
            # cl
            # one_out<-run.Rcrust(comps,c0,press,temp)


            # 	clusterApply(cl,c(1:3),run.Rcrust)

            # 	one<-c(comps,c0,c0=input_bulk[[y_i]][[x_i]],press=input_pt[[y_i]][[x_i]][1],temp=input_pt[[y_i]][[x_i]][2])
            # 			one_out<-run.Rcrust(comps,c0,press,temp)

            crust[[y_i]][[x_i]] <- run.Rcrust(
              comps = comps, c0 = input_bulk[[y_i]][[x_i]],
              press = input_pt[[y_i]][[x_i]][1],
              temp = input_pt[[y_i]][[x_i]][2],
              ph_extr_pnt = input_ph_extr[[y_i]][[x_i]],
              cumul_extract_pnt = cumul_extract[[y_i]][[x_i]],
              ph_add_pnt = input_ph_add[[y_i]][[x_i]],
              cumul_add_pnt = cumul_add[[y_i]][[x_i]]
            )
          }
        }
      }
    }
    # Save data
    save.image(file = paste0(projects_directory, "/", project_name, "/", working_file, ".RData"))
    cat("\n\nDone with calculations:\nResults saved to ", paste0(projects_directory, "/", project_name, "/", working_file, ".RData"), "\n\nSelect outputs through the Rcrust GUI or press esc to edit data in the Rconsole\n")
    flush.console()
  }
}
