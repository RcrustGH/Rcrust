#############################################
#
#  Initialize phase addition options
#
############################################
#handles definitions of the form: ph_add_definitions<-list("{1,1}_{10,1}"=c(2.95,0.618,15.3,7.11,5.17,0.15,0.183,6.44,63.5,5),"{1,1}_{10,1}"=c(2.95,0.618,15.3,7.11,5.17,0.15,0.183,6.44,63.5,2))
if(ph_add){
  cat("Setting phase addition options...\n")
  # fix-tag: error handling
  # Create data structure
  input_ph_add<-rep(list(rep(list(NULL),x_n)),y_n)
  #############################################
  #
  #  Phase addition from input definitions
  #
  ############################################
  cat("Setting phase additions from definitions in configuration file\n")
  pnts<-unlist(strsplit(names(ph_add_definitions),"_"))
  for(h in 1:length(ph_add_definitions)){
    a<-unlist(strsplit(gsub("\\{","",gsub("\\}","",pnts[h*2-1])),split=";"))
    b<-unlist(strsplit(gsub("\\{","",gsub("\\}","",pnts[h*2])),split=";"))
    if(length(ph_add_definitions[[h]])==length(c(all_elements,"mass"))){
    names(ph_add_definitions[[h]])<-c(all_elements,"mass")}
    for(x_i in a[1]:b[1]){
      for(y_i in a[2]:b[2]){
        if(is.null(input_ph_add[[y_i]][[x_i]])){
          input_ph_add[[y_i]][[x_i]]<-list(ph_add_definitions[[h]])
        }else{
          input_ph_add[[y_i]][[x_i]][[length(input_ph_add[[y_i]][[x_i]])+1]]<-ph_add_definitions[[h]]
        }
      }
    }
  }
}else{
  cat("No phase addition.\n")
}
cat("Done with phase addition options\n")
cat("..............................................\n")    