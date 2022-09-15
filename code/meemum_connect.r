# Launch & read meemum

########################################################
#
# Run meemum with parameters
#
########################################################
#function-def:run.meemum<-function(meemum.path="",build.file="",meemum.order="",press=press*1000,temp=temp+273.15,bulk="",components="",pt_comp=pt_comp)
run.meemum<-function(comps=comps,c0=c0,press=press,temp=temp,meemum_path=""){
					if(set_oxygen_fugacity){
					parms<-c(paste0(gsub("/code","/data",getwd()),"/parse_meem"),"n",c(temp+273.15,press*1000,as.numeric(input_bulk[[y_i]][[x_i]]["log10(fugacity)"])),"0","0","0")
					}else{
					parms<-c(paste0(gsub("/code","/data",getwd()),"/parse_meem"),"n",c(temp+273.15,press*1000),"0","0")
					}
					if (OperatingSystem == 'Windows') {
						# WINDOWS-ONLY VERSION (throws warnings on Linux)
						res<-system(shQuote(paste0(gsub("/code","/data",getwd()),"/",meemum_path)),invisible=T,input=parms,intern=T)
					} else {
						# LINUX-COMPATIBLE VERSION
						res<-system(shQuote(paste0(gsub("/code","/data",getwd()),"/",meemum_path)),input=parms,intern=T)
					}
return(res)
}
#####################################
# Parse meemum .prn file
#####################################
read.meemum<-function(meemum_in=meemum_in){
#Grab data matrix after line given line text
data_matrix_from_line<-function(line_text){
ln_given<-grep(line_text,meemum_in)
meemum_select<-meemum_in[-(1:ln_given)]
while(meemum_select[1]==""){
meemum_select<-meemum_select[-1]
}
nm_phases<-which(nchar(meemum_select)==0)[1]-2
text_given<-meemum_select[1:(nm_phases+1)]
text_altered<-trimws(text_given)
if(line_text=="Other Bulk Properties"){
text_cut<-unlist(strsplit(text_altered,"="))
mat_out<-trimws(text_cut)
}else{
text_altered<-gsub(" %","%",text_altered)
# mod-tag: remove system - fluid for now, may be useful in future to have access to this data
text_altered<-gsub("System - fluid","System-fluid",text_altered)
text_altered<-gsub("Poisson ratio","Poisson_ratio",text_altered)
text_altered<-strsplit(text_altered,split="\\s{1,}")

if(line_text=="Chemical Potentials"){text_altered[[2]]<-c("Bulk",text_altered[[2]])}
col_names<-c("Phase",text_altered[[1]])
mat_out<-matrix(unlist(text_altered)[-c(1:(length(col_names)-1))],nm_phases,byrow=TRUE)
colnames(mat_out)<-col_names
#number duplicates
for(ph in mat_out[which(duplicated(mat_out[,1])),1]){
      mat_out[which(mat_out[,1]==ph)[-1],1]<-paste0(ph,"_",1:length(which(mat_out[,1]==ph)[-1]))
}
}
return(mat_out)
}
#merge data tables by common first field with rows in same order
merge_matrix<-function(a,b){
full_rows<-union(a[,1],b[,1])
#add extra row of NA to any table that needs it
max_row<-max(nrow(a),nrow(b))
if(nrow(a)<max_row){
a<-rbind(a,matrix(NA,max_row-nrow(a),ncol(a)))
}
if(nrow(b)<max_row){
b<-rbind(b,matrix(NA,max_row-nrow(b),ncol(b)))
}
#rename first row
a[,1]<-full_rows
#delete key row in second table
b<-b[,-1]
#merge
return(cbind(a,b))
}
#grab crust output by merging tables
#reorder major_elements
ph_comp_mat<-data_matrix_from_line("Phase Compositions")
ph_comp_mat<-ph_comp_mat[,c(setdiff(colnames(ph_comp_mat),major_elements),major_elements),drop=FALSE]
#merge info
merge_1<-merge_matrix(ph_comp_mat,data_matrix_from_line("Molar Properties and Density"))
merge_2<-merge_matrix(merge_1,data_matrix_from_line("Seismic Properties"))
#make numeric and set phase column as row names
meemum_out<-merge_2[,-1]
meemum_out<-matrix(as.numeric(meemum_out),nrow(meemum_out))
rownames(meemum_out)<-merge_2[,1]
colnames(meemum_out)<-colnames(merge_2)[-1]
# mod-tag: remove system - fluid properties for now, may be useful in future to have access to this data
chk_call<-try(meemum_out["System-fluid",],silent=TRUE)
if(class(chk_call)!="try-error"){
meemum_out<-meemum_out[-grep("System-fluid",rownames(meemum_out)),]
}
#merge phases with the same name
if(merge_duplicates){
rownames(meemum_out)<-gsub("_[0-9]","",rownames(meemum_out))
while(any(duplicated(rownames(meemum_out)))){
ph<-rownames(meemum_out)[which(duplicated(rownames(meemum_out)))[1]]
thelines<-meemum_out[which(rownames(meemum_out)==ph),]
prop<-"wt%"
new<-.wtd.add(thelines,prop,ph)
meemum_out<-rbind(new,meemum_out[-which(rownames(meemum_out)==ph),])
}
}
#Create Bulk_rs with system
new_names<-rownames(meemum_out)
new_names[which(new_names=="System")]<-"Bulk_rs"
rownames(meemum_out)<-new_names
#meem_comps<-as.numeric(data_matrix_from_line("Bulk Composition")[,"wt%"])
#names(meem_comps)<-data_matrix_from_line("Bulk Composition")[,"Phase"]
#meemum_out["Bulk_rs",comps]<-meem_comps[comps]
meemum_out["Bulk_rs",comps]<-c0[comps]
#Set % to 100 for Bulk_rs
chk_call<-try(meemum_out["Bulk_rs","vol%"],silent=TRUE)
if(class(chk_call)!="try-error"){
meemum_out["Bulk_rs","wt%"]<-100
}
chk_call<-try(meemum_out["Bulk_rs","vol%"],silent=TRUE)
if(class(chk_call)!="try-error"){
meemum_out["Bulk_rs","vol%"]<-100
}
chk_call<-try(meemum_out["Bulk_rs","vol%"],silent=TRUE)
if(class(chk_call)!="try-error"){
meemum_out["Bulk_rs","mol%"]<-100
}
mass<-matrix(meemum_out[,"wt%"]/100*c0["mass"],nrow(meemum_out),1)
colnames(mass)<-"mass"
meemum_out<-cbind(mass,meemum_out)
#Add system properties to Bulk_rs
add_system_properties_to_bulk_rs<-TRUE
if(add_system_properties_to_bulk_rs){
Bulk_prop<-matrix(data_matrix_from_line("Other Bulk Properties"),2,)
colnames(Bulk_prop)<-Bulk_prop[1,]
mat_add<-matrix(0,nrow(meemum_out),ncol(Bulk_prop))
colnames(mat_add)<-Bulk_prop[1,]
meemum_out<-cbind(meemum_out,mat_add)
meemum_out["Bulk_rs",colnames(Bulk_prop)]<-as.numeric(Bulk_prop[2,colnames(Bulk_prop)])
}

#Nearest points average: grid interpolation {npav(grid_input,column1_val,column2_val)}
#finds points above and below in the x and y directions and returns average of 4 points (therefore resolution sets x:y aspect ratio)
npav<-function(grid_input,column1_val,column2_val){
#Error handling, press and temp must be contained within G table
contained<-TRUE
if(!press<=max(grid_input[,1])){contained<-FALSE}
if(!press>=min(grid_input[,1])){contained<-FALSE}
if(!temp<=max(grid_input[,2])){contained<-FALSE}
if(!temp>=min(grid_input[,2])){contained<-FALSE}
if(contained){
#Find closest 4 points as ur = upper right (x_above;y_above), lr = lower right (x_above;y_below), ul = upper left (x_below;y_above), ll = lower left (x_below;y_below)
rownames(grid_input)<-1:length(grid_input[,1])
#find equal or greater x
x_above_pnts<-which(grid_input[,2]>=temp)
#find least greater
x_min_above<-min(grid_input[x_above_pnts,2])
x_min_above_pnts<-which(grid_input[,2]==x_min_above)
	#ur
		#find equal or greater y
		x_a_y_a_pnts<-which(grid_input[x_min_above_pnts,1]>=press)
		#find least greater
		y_min_above<-min(grid_input[names(x_a_y_a_pnts),1])
		y_min_above_pnts<-which(grid_input[names(x_a_y_a_pnts),1]==y_min_above)
		ur<-grid_input[y_min_above_pnts[1],]
	#lr
		#find equal or lesser y
		x_a_y_b_pnts<-which(grid_input[x_min_above_pnts,1]<=press)
		#find least lesser
		y_max_below<-max(grid_input[names(x_a_y_b_pnts),1])
		y_max_below_pnts<-which(grid_input[names(x_a_y_b_pnts),1]==y_max_below)
		lr<-grid_input[y_max_below_pnts[1],]
#find equal or lesser x
x_below_pnts<-which(grid_input[,2]<=temp)
#find least lesser
x_max_below<-max(grid_input[x_below_pnts,2])
x_max_below_pnts<-which(grid_input[,2]==x_max_below)
	#ul
		#find equal or greater y
		x_b_y_a_pnts<-which(grid_input[x_max_below_pnts,1]>=press)
		#find least greater
		y_min_above_x_below<-min(grid_input[names(x_b_y_a_pnts),1])
		y_min_above_x_below_pnts<-which(grid_input[names(x_b_y_a_pnts),1]==y_min_above_x_below)
		ul<-grid_input[y_min_above_x_below_pnts[1],]
	#ll
		#find equal or lesser y
		x_b_y_b_pnts<-which(grid_input[x_max_below_pnts,1]<=press)
		#find least lesser
		y_max_below_x_below<-max(grid_input[names(x_b_y_b_pnts),1])
		y_max_below_x_below_pnts<-which(grid_input[names(x_b_y_b_pnts),1]==y_max_below_x_below)
		ll<-grid_input[y_max_below_x_below_pnts[1],]
#take average of 4 corners
return(mean(ur[3],lr[3],ul[3],ll[3]))
}else{
cat(paste0("Error activity cannot be calculated for ",activities_to_calculate[G_i]," as press or temp lie outside of the gibbs enthalpy provided\n"))
return(0)
}
}


#Calculate activities
#modtag - add activities
#a=e^((u-u0)/RT)
#u0 = chemical potential of the pure substance at P,T since G=u for a pure substance (J/mol)
#u = chemical potential of the component (i.e. TiO2) at P,T (J/mol)
#T = Temperature in Kelvin
#R = Molar gas constant = 8.3144598 m2 kg s-2 K-1 mol-1
if(calculate_activities){
#fixtag - can't use this temp and press statement as it erronously gives another value, fix this
press<-input_pt[[y_i]][[x_i]][1]
temp<-input_pt[[y_i]][[x_i]][2]
Chemical_potentials<-data_matrix_from_line("Chemical Potentials")
#modtag - convert to upper to ensure case matching
colnames(Chemical_potentials)<-toupper(colnames(Chemical_potentials))
activities_to_calculate<-toupper(activities_to_calculate)
for(G_i in 1:length(activities_to_calculate)){
if(FALSE){
#find u0 from pure phase
u0<-meemum_out["ru","G(J)"]   #(J) as molar property
#Interpolate G of pure phase from input grid
u0<-npav(input_G[[G_i]],press,temp)
which(input_G[[G_i]][,1]==press)
}else{
#find G in table
#fixtag - equivalence statements on press did not find equality between press=12.6 and 12.6, had to add plaster to round to significant figures of its original value to allow equivalance
f <- function(x){length(gregexpr("[[:digit:]]", as.character(x))[[1]])}
press_pnts<-which(input_G[[G_i]][,1]==signif(press,f(press)))
temp_pnts<-which(input_G[[G_i]][,2]==temp)
G_pnt<-intersect(press_pnts,temp_pnts)
if(length(G_pnt)!=1){
cat("Error non unique G found in pure phase enthalpy table")
stop()}
u0<-input_G[[G_i]][G_pnt,3]
}
u<-as.numeric(Chemical_potentials[1,activities_to_calculate[G_i]])  #(J/mol)
R<-8.3144598
aA<-exp((u-u0)/(R*(temp+273.15)))
#aB<-exp(-((u0-u)/(R*(temp+273.15))))
activities<-matrix(0,nrow(meemum_out),1)
colnames(activities)<-paste0("a",activities_to_calculate[G_i])
meemum_out<-cbind(meemum_out,activities)
meemum_out["Bulk_rs",paste0("a",activities_to_calculate[G_i])]<-aA
}
}
return(meemum_out)
}
