#############################################
#
#  Create a meemum build file
#
############################################
# temporary setting to FALSE
if (class(try(set_oxygen_fugacity, silent = TRUE)) == "try-error") {
  set_oxygen_fugacity <- FALSE
}

cat("Creating meemum build file...\n")
thermo_in <- readLines(paste0(gsub("/code", "/data", getwd()), "/", thermodynamic_data_file))
comp_start <- grep("begin_components", thermo_in)
comp_end <- grep("end_components", thermo_in)
number_components <- comp_end - comp_start - 1
# fixtag must allow user to specify mol or wt and should be read to perplex option file
molar_vs_wt <- 1
# new_name must be less than 6 characters, MUST BE ALL CAPS
# can only build components with total of 11 or less components
# replace components each time
transformations <- NULL
n_comp_trans <- 0
if (exists("comp_transformations")) {
  if (length(comp_transformations) >= 1 & !all(comp_transformations == "")) {
    n_comp_trans <- length(comp_transformations)
  }
}
if (n_comp_trans > 0) {
  # modtag check that avilable components refreshed in GUI
  available_components <- NULL
  suppressWarnings(qq_try <- try(scan(file = gsub("Rcrust/code", paste0("Rcrust/data/", thermodynamic_data_file), getwd()), what = "character", sep = "\n", quiet = TRUE), silent = TRUE))
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
  for (i in 1:length(comp_transformations)) {
    new_comp <- strsplit(names(comp_transformations)[i], split = "_")[[1]][2]
    old_comp <- strsplit(names(comp_transformations)[i], split = "_")[[1]][1]
    transformations <- c(transformations, paste(new_comp, paste(rep(" ", 6 - nchar(new_comp)), collapse = ""), which(available_components == old_comp), " component transformation\n", paste(strsplit(comp_transformations[[i]], ",")[[1]], collapse = "   "), sep = ""))
  }
}
##################
# Create components matrix
comp_mat <- NULL
for (ox in major_elements) {
  comp_mat <- c(comp_mat, paste(ox, "     1  1.00000  0.00000  ", if (molar_vs_wt == 1) {
    "weight amount"
  } else {
    "molar amount"
  }, "\n", sep = ""))
}
# Sean-tag
if (component_packet) {
  if (!isTRUE(new_components == "")) {
    for (i in 1:length(new_components)) {
      comp_mat <- comp_mat[-grep(new_components[i], comp_mat)]
    }
  }
}
# Write dummy meemum build file
dummy <- paste(
  gsub("/code", "/data", getwd()), "/", thermodynamic_data_file, "     thermodynamic data file
no_print | print generates print output
no_plot     | no_plot suppresses plot output
",
  gsub("/code", "/data", getwd()), "/", solution_models_file, "     solution model file, blank = none
parse_meem
",
  gsub("/code", "/data", getwd()), "/", perplex_option_file, "     computational option file
    5 calculation type: 0 - composition, 1 - Schreinemakers, 3 - Mixed, 4 - gwash, 5 - gridded min, 7 - 1d fract, 8 - gwash 9 - 2d fract, 10 - 7 w/file input
    0 unused place holder, post 06
    0 unused place holder, post 06
    0 unused place holder, post 06
    0 unused place holder, post 06
    0 unused place holder, post 06
    0 unused place holder, post 06
    0 unused place holder, post 06
    0 unused place holder, post 06
    0 unused place holder, post 06
    ", n_comp_trans, " number component transformations\n    ",
  number_components, " number of components in the data base\n",
  paste(transformations, collapse = "\n"), "\n",
  molar_vs_wt, " component amounts, 0 - molar, 1 weight
    0 unused place holder, post 06
    0 unused place holder, post 06
    0 unused place holder, post 05
    5 ifug EoS for saturated phase
    2 gridded minimization dimension (1 or 2)
    0 special dependencies: 0 - P and T independent, 1 - P(T), 2 - T(P)
 0.00000      0.00000      0.00000      0.00000      0.00000     Geothermal gradient polynomial coeffs.

begin thermodynamic component list
", paste(comp_mat, collapse = ""),
  "end thermodynamic component list

begin saturated component list
", paste(unlist(strsplit(saturated_components, split = ",")), collapse = "\n"), "
end saturated component list


begin saturated phase component list
", paste(unlist(strsplit(saturated_phase_components, split = ",")), collapse = "\n"), "
end saturated phase component list


begin independent potential/fugacity/activity list
", if (set_oxygen_fugacity) {
    paste0("O2    f_O2     O2      ", "\n")
  },
  paste(unlist(strsplit(independent_potential_fugacity_activity, split = ",")), collapse = "\n"), "
end independent potential list


begin excluded phase list
", paste(exclude_phases, collapse = "\n"), "
end excluded phase list


begin solution phase list
", paste(use_sol_models, collapse = "\n"),
  "\nend solution phase list

 0.0000      0.0000     0.00000000  0.0000      0.0000     max p, t, xco2, u1, u2
 0.0000      0.0000     0.00000000  0.0000      0.0000     min p, t, xco2, u1, u2
 0.0000      0.0000     0.00000000  0.0000      0.0000     unused place holder post 06

 2  1  4  5  3   indices of 1st & 2nd independent & sectioning variables
",
  sep = ""
)
write(dummy, file = paste0(gsub("/code", "/data", getwd()), "/parse_meem.dat"))
cat("Created meemum build file as ", paste0(gsub("/code", "/data", getwd()), "/parse_meem.dat"), "\n", sep = "")
cat("..............................................\n")
