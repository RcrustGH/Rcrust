#############################################
#
#  Wrapper function containing compiled form of meemum from Perple_X (Give citation to Lars and Jamie here)
#
############################################
#function-def:wrapper<-function(comps,c0,press,temp,calc_choice="read.meemum")
wrapper<-function(comps,c0,press,temp,calc_choice="read.meemum"){
#############################################
#
# Wrapper function, gives some exception handling
# So an error in the calc does not necessarily kill all the previous results !
#
############################################
if(!(calc_choice=="read.meemum"|calc_choice=="lars.wrap")){
cat("\nError: ",calc_choice,"Is not an appropriate calculation choice (calc_choice)\n")
stop()
}
if(calc_choice=="read.meemum"){
#control comps order
input_file<-readLines(paste0(gsub("Projects","data",projects_directory),"/parse_meem.dat"))
comp_mat<-NULL
 for(ox in major_elements){
 comp_mat<-c(comp_mat,paste0(ox,"     1  ",c0[which(major_elements==ox)],"  0.00000  ",if(molar_vs_wt==1){"weight amount"}else{"molar amount"},sep=""))
 }
start_write<-grep("begin thermodynamic component list",input_file)+1
end_write<-grep("end thermodynamic component list",input_file)-1
input_file[start_write:end_write]<-comp_mat
write(input_file, file = paste0(gsub("Projects","data",projects_directory),"/parse_meem.dat"))
#fix-tag normalise c0  here but should ensure mass balnace at extraction
c0[major_elements]<-c0[major_elements]/sum(c0[major_elements])*100
# mod-tag: have to bring c0 into global, look for alternative
c0<<-c0
meem_calc<-run.meemum(comps,c0,press,temp,meemum_path)
#Write meemum output to file?
if(!exists("export_meemum_output")){export_meemum_output<-FALSE}
if(export_meemum_output){sink("output_meem_mohit.txt", append = TRUE)}
#Write meemum output to console?
if(!exists("print_meem")){print_meem<-FALSE}
if(print_meem){print(meem_calc)}
calc_out<-read.meemum(meem_calc)
}
# mod-tag: Do we still need "lars.wrap"?
if(calc_choice=="lars.wrap"){
#mod-tag: check if need to normalise c0 or if perplex does this
#a (Options: P, T, bulk composition in given order)
  a <- .Call("R_phaseq", as.numeric(press*1000), as.numeric(temp+273.15), as.numeric(c0[comps]))
  #'a' is a list, with following components:
#    1 - Return value (should be zero in case of successful minimization
#    2 - Amounts of stable phases (wt% or mol, depending on PerpleX options)
#    3 - Names of the stable phases
#    4 - Compositions of stable phases (a[[4]][[i]][j] for phase a[[3]][[i]],
#        element comps[j] (wt% or mol, depending on PerpleX options)* PerpleX options not editable
#	 5 - Bulk properties
#	 6 - System properties of stable phases (a[[6]][[i]][j] for phase a[[3]][[i]], system property [j]
#Check for error
error<<-FALSE
chk_call<-try(length(a[[3]])==length(a[[4]]),silent=TRUE)
if(class(chk_call)=="try-error"){
error<<-TRUE
}
else{
if(!chk_call){
error<<-TRUE
}
}
out_colnames<<-c("wt%","vol%",comps,"mass","V(J/bar)","H(J)","Gruneisen_T","Ks(bar)","Mu(bar)","V0(km/s)","Vp(km/s)","Vs(km/s)","Vp/Vs","Rho(kg/m3)","Cp(J/K)","alpha(1/K)","beta(1/bar)","S(J/K)","N(g)","Cp/Cv")
if(!error){
  #compile phases
phases<-NULL
for(ii in 1:length(a[[3]])){
one<-matrix(a[[4]][[ii]],1,length(a[[4]][[ii]]))
prop_one<-matrix(a[[6]][[ii]],1,length(a[[6]][[ii]]))
phases<-rbind(phases,cbind(one,prop_one))
}
colnames(phases)<-c(comps,"V(J/bar)","H(J)","Gruneisen_T","Ks(bar)","Mu(bar)","V0(km/s)","Vp(km/s)","Vs(km/s)","Vp/Vs","Rho(kg/m3)","?","Cp(J/K)","alpha(1/K)","beta(1/bar)","S(J/K)","??","N(g)",17,18,19,20,21,22,23,24,25,26,"Cp/Cv")
# Create mass and wt columns
wt<-matrix(a[[2]],length(a[[3]]),1)
colnames(wt)<-"wt%"
mass<-matrix(a[[2]]/100*c0["mass"],length(a[[3]]),1)
colnames(mass)<-"mass"
# Create vol% columns
vol<-(wt/phases[,"Rho(kg/m3)",drop=FALSE])/(100/a[[5]][10])*100
colnames(vol)<-"vol%"
#Create Bulk
bulk_row<-matrix(c(100,100,c0[],a[[5]]),1)
rownames(bulk_row)<-"Bulk_rs"
# Compile calc_out
calc_out<-cbind(wt,vol,phases[,names(c0[1:length(comps)])],mass,phases[,(length(comps)+1):ncol(phases)])
rownames(calc_out)<-a[[3]]
calc_out<-rbind(calc_out,bulk_row)
calc_out<-calc_out[,out_colnames,drop=FALSE]
}else{
# Return Blank Comp
colmns<-c("wt%",out_colnames)
calc_out<-matrix(0,2,length(out_colnames))
rownames(calc_out)<-c("Error","Bulk_error")
colnames(calc_out)<-out_colnames
}
}
return(calc_out)
}
