#############################################
#
#  Initialize phase extraction options
#
############################################
# handles definitions of the form: ph_extr_definitions<-list("{1,1}_{6,1}"=c(condition="melt(HP)>=8wt.%",Pl="20%",ky="3","melt(HP)"="rt(4,wt.%)"),"{4,1}_{9,1}"=c(condition="melt(HP)>=8wt.%",Pl="20%",ky="3","melt(HP)"="rt(4,wt.%)"))
if (ph_extr) {
  cat("Setting phase extraction options...\n")
  # fix-tag: error handling
  # Create data structure
  input_ph_extr <- rep(list(rep(list(NULL), x_n)), y_n)
  #############################################
  #
  #  Phase extraction from input definitions
  #
  ############################################
  cat("Setting phase extractions from definitions in configuration file\n")
  pnts <- unlist(strsplit(names(ph_extr_definitions), "_"))
  for (h in 1:length(ph_extr_definitions)) {
    a <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", pnts[h * 2 - 1])), split = ";"))
    b <- unlist(strsplit(gsub("\\{", "", gsub("\\}", "", pnts[h * 2])), split = ";"))
    for (x_i in a[1]:b[1]) {
      for (y_i in a[2]:b[2]) {
        if (is.null(input_ph_extr[[y_i]][[x_i]])) {
          input_ph_extr[[y_i]][[x_i]] <- list(ph_extr_definitions[[h]])
        } else {
          input_ph_extr[[y_i]][[x_i]][[length(input_ph_extr[[y_i]][[x_i]]) + 1]] <- ph_extr_definitions[[h]]
        }
      }
    }
  }
} else {
  cat("No phase extraction.\n")
}
cat("Done with phase extraction options\n")
cat("..............................................\n")
