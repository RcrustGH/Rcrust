#############################################
#
#  Initialize Gibbs enthalpy of pure phases for activity calculations
#
############################################
cat("Initializing Gibbs enthalpy of pure phases for activity calculations...\n")
# functions needed
# mod-tag: Surely there's a better way to make a matrix numeric?
make_numeric_mat <- function(x) {
  y <- matrix(as.numeric(as.vector(as.matrix(x))), nrow(x), )
  rownames(y) <- rownames(x)
  colnames(y) <- colnames(x)
  return(y)
}
#############################################
#
#  G of pure phases from input files
#
############################################
data_directory <- gsub("code", "data", getwd())
G_files <- unlist(strsplit(G_pure_phases, ";"))
# Extract component names from G_list
G_list <- unlist(strsplit(G_files, "_"))
activities_to_calculate <- G_list[seq(2, length(G_list), 3)]
input_G <- vector("list", length(G_files))
for (G in 1:length(G_files)) {
  G_table <- read.table(file = paste(data_directory, G_files[G], sep = "/"), sep = "\t")
  # convert to matrix
  G_table <- as.matrix(G_table)
  # error validation
  # Check headers
  if (!all(G_table[1, 1] == "Pressure(kbar)", G_table[1, 2] == "Temperature(C)", G_table[1, 3] == "G(J)")) {
    cat(paste0("Error in Gibbs enthalpy inititation for activity calculations, incorrect column names in ", G, "\n"))
    stop()
    flush.console()
  }
  # remove headers
  G_table <- G_table[-1, ]
  # make numeric
  G_table <- make_numeric_mat(G_table)
  input_G[[G]] <- G_table
  cat(paste0("Gibbs enthalpy successfully read from ", G_files[G], "\n"))
}
cat("Done with activity calculation preperations\n")
cat("..............................................\n")
