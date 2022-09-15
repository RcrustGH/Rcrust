#Break on comma unless enclosed in quotes
break_on_comma<-function(x){
#find commas
commas<-which(strsplit(x,"")[[1]]==",")
#name commas
names(commas)<-rep("break",length.out=length(commas))
#find quotes
quotes<-which(strsplit(x,"")[[1]]=="\"")
#number quotes
names(quotes)<-rep(c("odd","even"),length.out=length(quotes))
#if comma is between an odd and an even quotes set to not break
if(length(quotes>0)){
new_names<-NULL
for(comma in commas){
new_name<-"break"
if(comma>min(quotes)&comma<max(quotes)){
below<-which(quotes<comma)[length(which(quotes<comma))]
above<-which(quotes>comma)[1]
if(names(below)=="odd"&names(above)=="even"){
new_name<-"dont"
}
}
new_names<-c(new_names,new_name)
}
names(commas)<-new_names
}
#break by commas tagged as break
if(any(names(commas)=="break")){
breaks<-c(-1,commas[which(names(commas)=="break")]-1,nchar(x))
new_x<-NULL
for(i in 1:(length(breaks)-1)){
new_x<-c(new_x,substr(x,breaks[i]+2,breaks[i+1]))
}
x<-new_x
}
return(x)
}


#Plot phase boundaries
if(TRUE){
      all_columns<-NULL
	  all_phases<-NULL
      for(y_i in 1:length(crust_out())){
        for(x_i in 1:length(crust_out()[[1]])){
		  all_columns<-union(all_columns,colnames(crust_out()[[y_i]][[x_i]]))
          all_phases<-union(all_phases,rownames(crust_out()[[y_i]][[x_i]]))
        }
      }
library(RColorBrewer)
library(grDevices)
grDevices::postscript(file="phase_boundaries.ps",onefile=TRUE,horizontal=TRUE)
cols<-RColorBrewer::brewer.pal(length(all_phases), "Spectral")
cols<-rep(cols,length.out=length(all_phases))
names(cols)<-all_phases 
graphics::contour(t(flip_y(grid_data("wt%","Bulk_rs",crust_out())[[2]])))	    
for(ph in all_phases){
grid_data<-t(flip_y(grid_data("wt%",ph,crust_out())[[2]]))
graphics::contour(grid_data,levels=c(0,0.001),labels=ph,add=TRUE,col=cols[ph])
}	  
dev.off()
}

outfile_path<-paste0(sub("/code","/Projects",getwd()),"/",working_file,"_melt cont.ps")
postscript(file=outfile_path,onefile=TRUE,horizontal=TRUE)
plot(cont_mat,col=cont_mat[,3],xlim=c(640,920),ylim=c(2,12))

	  
	  
contour(grid_data,levels=c(0,0.001),labels=ph,add=TRUE,col=cols[ph])	  
	  
contours_phs<-NULL	  
for(ph in all_phases){
grid_data<-t(flip_y(grid_data("wt%",ph,crust_out())[[2]]))
if(all(grid_data!=0)){
a<-list(level=ph,"x"=0,"y"=0)
}else{	  
a<-grDevices::contourLines(grid_data,levels=c(0,0.001))
for(i in 1:length(a)){
a[[i]][[1]]<-ph
}
}
contours_phs<-c(contours_phs,a)
}


contour(t(flip_y(grid_data("wt%","Bulk_rs",crust_out())[[2]])))
lines(contours_phs)

for(ph in 1:length(all_phases)){
browser()
lines(contours_phs[[ph]]$x,contours_phs[[ph]]$y)
}



for(ph in all_phases){
grid_data<-t(flip_y(grid_data("wt%",ph,crust_out())[[2]]))
contour()
}



contour(t(flip_y(grid_data_r()[[2]])),



graphics::plot.window(c(0,nrow(grid_data)),c(0,ncol(grid_data)))
graphics::plot.new()
lines(contours_phs[[1]])
contours_phs<-NULL	  
for(ph in all_phases){
grid_data<-t(flip_y(grid_data("wt%",ph,crust_out())[[2]]))
if(all(grid_data!=0)){
a<-list(level=ph,"x"=0,"y"=0)
}else{	  
a<-grDevices::contourLines(grid_data,levels=c(0,0.001))
for(i in 1:length(a)){
a[[i]][[1]]<-ph
}
}
contours_phs<-c(contours_phs,a)
}

	  
t(flip_y(grid_data("wt%",all_phases[1],crust_out())))	  




t(flip_y(grid_data("H2O","Bt_rs",crust_out())[[2]]))


axis(1,(0:(length(bottom[[2]])-1))/(length(bottom[[2]])-1),bottom[[2]])}else


#Convert upper case to mixed case for major elements
		convert_case<-function(x){
		mix_case<-c("SiO2","K2O","TiO2","FeO","Fe2O3","O2","CaO","MgO","Na2O","Al2O3","H2O","MnO","NiO","ZrO2","Cl2","CO2","CuO","Cr2O3","S2","F2")
		converted_case<-NULL
		for(i in 1:length(x)){
		converted_case<-c(converted_case,mix_case[which(toupper(mix_case)==toupper(x[i]))])
		}
		return(converted_case)
		}


#Compare trends
gradients<-NULL
for(y_i in c(-1,1:y_n)){
#Select group
subgroup<-WR[which(WR[,"y_i"]==y_i),]
#Calculate linear regression for each path
model<-lm(K2O~SiO2,data=data.frame(subgroup))
gradients<-c(gradients,coef(model)[2])
}


for(y_i in 2:length(gradients)){
abline(coef=c(0,gradients[y_i]), col = "gray60",add=TRUE)
}
y_i<-1
abline(coef=c(0,gradients[y_i]), col = "red")

lines<-matrix(c(gradients,rep(0,y_n+1)),y_n+1,2)



#Compare trends
cols<-heat.colors(y_n+1)
gradients<-NULL
for(y_i in c(-1,1:y_n)){
#Select group
subgroup<-WR[which(WR[,"y_i"]==y_i),]
#Calculate linear regression for each path
model<-lm(K2O~SiO2,data=data.frame(subgroup))
abline(model,col=cols[y_i])
}




#Compare trends
comps<-c("K2O","Al2O3","FeO")
y_n<-72
gradients<-matrix(c(-1,1:y_n),y_n+1,1)
colnames(gradients)<-"path"
for(i in 1:length(comps)){
gradients_comp<-matrix(,y_n+1,1)
for(y_i in c(-1,1:y_n)){
#Select group
subgroup<-WR[which(WR[,"y_i"]==y_i),]
#Calculate linear regression for each path
model<-lm(eval(parse(text=paste0(comps[i],"~SiO2"))),data=data.frame(subgroup))
gradients_comp[y_i,]<-coef(model)[2]
colnames(gradients_comp)<-comps[i]
}
gradients<-cbind(gradients,gradients_comp)
}
gradients_normalised<-gradients
for(i in 1:(length(comps)+1)){
gradients_normalised[,i]<-gradients[,i]/gradients[1,i]
}
gradients_distance<-gradients
for(i in 1:(length(comps)+1)){
gradients_distance[,i]<-gradients_normalised[,i]-1
}
gradients_distance<-abs(gradients_distance)
gradients_distance[,1]<-gradients[,1]
#sort by comp1
sort_comp1<-gradients_distance[sort.list(gradients_distance[,2]),]
top_comp1<-sort_comp1[2:11,1]
#sort by comp2
sort_comp2<-gradients_distance[sort.list(gradients_distance[,3]),]
top_comp2<-sort_comp2[2:11,1]
#sort by comp3
sort_comp3<-gradients_distance[sort.list(gradients_distance[,4]),]
top_comp3<-sort_comp3[2:11,1]
intersect(intersect(top_comp1,top_comp2),top_comp3)



#Plot phase abundance versus path cumulative column graph
library(RColorBrewer)
library(graphics)
#All RColorBrewer palettes display.brewer.all()
#Display a specific pallette display.brewer.pal(12,"Set3")

#Normalise to percentage
phase_abundance_data<-apply(phase_abundance_r()[[2]], 2, function(x){x*100/sum(x,na.rm=T)})
#Set colour pallette
cols<-RColorBrewer::brewer.pal(min(nrow(phase_abundance_data),12), "Set3") 
#Plot
barplot(phase_abundance_data,space=0,col=cols,border=NA,legend.text=TRUE)




#select = 20

limit<-30
melts_below_limit<-NULL
below_limits<-NULL
for(y_i in 1:y_n){
for(x_i in 1:x_n){
#if less than x% of system consumed
chk_H2O_extr_cumul<-try(crust[[y_i]][[x_i]]["Melt_es_cumul","mass"],silent=TRUE)
chk_melt_extr_cumul<-try(crust[[y_i]][[x_i]]["H2O_es_cumul","mass"],silent=TRUE)
if(class(chk_H2O_extr_cumul)!="try-error"){
if(class(chk_melt_extr_cumul)!="try-error"){
total_melt_loss_anhydrous<-crust[[y_i]][[x_i]]["Melt_es_cumul","mass"]/(100-crust[[y_i]][[x_i]]["H2O_es_cumul","mass"])*100
if(total_melt_loss_anhydrous<limit){
#if melt extract exists
chk_melt_extract<-try(crust[[y_i]][[x_i]]["Melt_es","mass"],silent=TRUE)
if(class(chk_melt_extract)!="try-error"){
#grab comp
melts_below_limit<-rbind(melts_below_limit,crust[[y_i]][[x_i]]["Melt_es",])
below_limits<-rbind(below_limits,matrix(c(x_i,y_i),1,2))
}
}
}
}
}
}
write.csv(below_limits,"test.csv",row.names = TRUE)



below_limits<-as.matrix(read.csv("C:/Rcrust/code/test.csv"))
WR_select<-NULL
melt_names<-NULL
for(xrow in 1:nrow(WR)){
if(any(WR[xrow,"y_i"]==below_limits[,3])){
from_limits<-below_limits[which(below_limits[,3]==WR[xrow,"y_i"]),2]
from_WR<-WR[which(WR[,"y_i"]==below_limits[xrow,3]),"x_i"]
keep<-NULL
for(i in 1:length(from_WR)){
if(any(from_WR[i]==from_limits)){
keep<-c(keep,from_WR[i])
}
}
melt_names<-c(melt_names,names(keep))
}
#WR_select<-rbind(WR_select,)
}
WR_rows<-NULL
for(i in 1:length(melt_names)){
WR_rows<-c(WR_rows,which(rownames(WR)==melt_names[i]))
}
WR_select<-WR[WR_rows,]
rownames(WR_select)<-paste0(rownames(WR_select),"_select")
write.csv(WR_select,"C:/Rcrust/code/out.csv",row.names = TRUE)



rownames(full_data)<-full_data[,"ID"]

below_limits<-as.matrix(read.csv("C:/Rcrust/code/test.csv"))
WR_select<-NULL
melt_names<-NULL
for(xrow in 1:nrow(WR)){
if(any(WR[xrow,"y_i"]==below_limits[,3])){
from_limits<-below_limits[which(below_limits[,3]==WR[xrow,"y_i"]),2]
from_WR<-WR[which(WR[,"y_i"]==below_limits[xrow,3]),"x_i"]
keep<-NULL
for(i in 1:length(from_WR)){
if(any(from_WR[i]==from_limits)){
keep<-c(keep,from_WR[i])
}
}
melt_names<-c(melt_names,names(keep))
}
#WR_select<-rbind(WR_select,)
}
WR_rows<-NULL
for(i in 1:length(melt_names)){
WR_rows<-c(WR_rows,which(rownames(WR)==melt_names[i]))
}
WR_select<-WR[WR_rows,]
rownames(WR_select)<-paste0(rownames(WR_select),"_select")
write.csv(WR_select,"C:/Rcrust/code/out.csv",row.names = TRUE)


clusters<-groups
border<-NULL
fill<-FALSE
precision<-50

#center of polygon
 on.exit(options(show.error.messages = TRUE))
    if (is.numeric(clusters)) {
        winDialog(type = "ok", paste("clusters not defined!"))
        return()
    }
    labs <- as.character(levels(factor(clusters)))
    labs <- abbreviate(labs, minlength = 5)
    if (is.null(border)) {
        border <- .getAvgColour(clusters)
    }
    names(border) <- labs
    gr <- abbreviate(clusters, minlength = 5)
    ee <- lapply(labs, function(i) {
        xx <- x.data[gr == i]
        yy <- y.data[gr == i]
        xx <- xx[!is.na(xx) & !is.na(yy)]
        yy <- yy[!is.na(xx) & !is.na(yy)]
        if (length(xx) > 0 & length(yy) > 0) {
            picked <- intersect(names(xx), names(yy))
            grid <- kde2d(xx[picked], yy[picked], n = 500, lims = c(sheet$d$c$xlim, 
                sheet$d$c$ylim))
            contour(grid, levels = max(unlist(grid$z)/precision), 
                col = border[i], labels = i, drawlabels = TRUE, 
                add = TRUE, ...)
        }
    })
    invisible()



#Write data file
write_data_file(data_file(crust,x_n,y_n,c("ID","Phase","y_i","x_i","Pressure(kbar)","Temperature(C)","wt%",comps,"mass"),"All","All"),working_file,projects_directory,".csv")

#Load Grid
#matt tag
grid_temp_def<-"600+2*(x_i-1)"
grid_press_def<-"2.5+0.07*(y_i-1)"
outfile<-paste0(gsub("code","Projects/",getwd()),working_file,"/Outputs/",working_file,"_grid.rds")
store_r$crust_r<-readRDS(outfile)
store_r$crust_r<<-readRDS(outfile)
crust<-readRDS(outfile)
crust<<-readRDS(outfile)
x_n<-length(crust[[1]])
x_n<<-length(crust[[1]])
y_n<-length(crust)
y_n<<-length(crust)
input_pt<-rep(list(rep(list(c(Pressure=NULL,Temperature=NULL)),x_n)),y_n) 
for(y_i in 1:y_n){
for(x_i in 1:x_n){
input_pt[[y_i]][[x_i]]<-matrix(c(eval(parse(text=grid_press_def)),eval(parse(text=grid_temp_def))),nrow=1,dimnames = list("PT",c("Pressure","Temperature")))
}
}
store_r$input_pt_r<-input_pt
store_r$input_pt_r<<-input_pt
input_pt<<-input_pt

if(exists("crust")){store_r$crust_r<-crust}
if(exists("input_pt")){store_r$input_pt_r<-input_pt}


#Function  - is.real tests if number is a finite real number between the numbers "lower" and "upper"
#NA NULL NaN error Inf -Inf logical(0) interger(0) numeric(0)
#doesnt pick up c(NULL,1)
is.real<-function(x,lower=-Inf,upper=Inf){
chk.real<-function(x,lower=-Inf,uppper=Inf){
if(length(x)==0){return(FALSE)}else{
for(i in 1:length(x)){
if(length(x[i])!=1){return(FALSE)}
if(suppressWarnings(is.null(x[i]))){return(FALSE)}
if(suppressWarnings(is.na(x[i]))){return(FALSE)}
if(suppressWarnings(is.na(suppressWarnings(as.numeric(x[i]))))){return(FALSE)}
if(!x[i]<=uppper){return(FALSE)}
if(!x[i]>=lower){return(FALSE)}
}
}
return(TRUE)
}
chk<-try(chk.real(x,lower,upper),silent=TRUE)
if(class(chk)=="try-error"){return(FALSE)}else{return(chk)}
}

#Function  - num returns a finite numbeic form of x or 0
num<-function(x){
if(is.null(x)){return(0)}
if(is.na(x)){return(0)}
suppressWarnings(if(is.na(as.numeric(x))){return(0)})
return(x)
}

if(FALSE){
#asp_fac<-((1.582/2.054)*10+250)/10
#19.935/1.38393
#14,40463:1
#0.7995115/14.40463
}

#Given T and P defs plot points
#modtag - allow this as "view" from Rcrust GUI
if(FALSE){
T<-"cos((-8.53076560994813+0.0782789097942497*(y_i-1))*pi/180)*(x_i-1)+650"
P<-"sin((-8.53076560994813+0.0782789097942497*(y_i-1))*pi/180)*(x_i-1)+10"
x_n<-451
y_n<-201
pt_mat<-NULL
for(y_i in 1:y_n){
#for(x_i in 1:x_n){
for(x_i in c(1,seq(10,x_n,30))){
pt_mat<-rbind(pt_mat,matrix(c(eval(parse(text=paste(T))),eval(parse(text=paste(P)))),1,2))
}
}
graphics::plot(pt_mat,xlim=c(0,1100),ylim=c(0,25))

#Lines
graphics::plot(c(0,0),xlim=c(0,1100),ylim=c(0,25))
lines(pt_mat)
}


#Plot Paths in P-T space
if(FALSE){
T0<-650
P0<-10
T<-"cos((-8.53076560994813+0.0782789097942497*(y_i-1))*pi/180)*(x_i-1)+650"
P<-"sin((-8.53076560994813+0.0782789097942497*(y_i-1))*pi/180)*(x_i-1)+10"
x_i<-451
y_n<-201
graphics::plot(c(T0,P0),xlim=c(0,1100),ylim=c(0,25))
for(y_i in 1:y_n){
graphics::segments(T0,P0,eval(parse(text=paste(T))),eval(parse(text=paste(P))))
}
}


if(FALSE){
x_scaling<-70
#Plot paths with H2O line
T0<-650
P0<-5
T<-"cos((-2.86240522611175+0.0779363487601036*(y_i-1))*pi/180)*(x_i-1)+650"
P<-"sin((-2.86240522611175+0.0779363487601036*(y_i-1))*pi/180)*(x_i-1)+5"
x_i<-1101
y_n<-111



#H2O plot
H2O_at_solidus<-matrix(NA,y_n,1)
for(y_i in 1:y_n){
for(x_i in 1:x_n){
chk<-try(crust[[y_i]][[x_i]]["Melt_rs","mass"],silent=TRUE)
if(class(chk)!="try-error"){
H2O_at_solidus[y_i,1]<-crust[[y_i]][[x_i]]["Bulk_rs","H2O"]
break}
}
}
plot_points<-NULL
for(y_i in 1:y_n){
plot_y<-eval(parse(text=paste(P)))
plot_x<-eval(parse(text=paste(T)))+x_scaling*H2O_at_solidus[y_i,1]
plot_points<-rbind(plot_points,matrix(c(plot_x,plot_y),1,2))
}
graphics::plot(plot_points,xlim=c(0,1300),ylim=c(0,30))
graphics::plot(plot_points)
lines(plot_points)

#PT lines
for(y_i in 1:y_n){
graphics::segments(T0,P0,eval(parse(text=paste(T))),eval(parse(text=paste(P))))
}
#horizontal lines
hor_length<-x_scaling*max(H2O_at_solidus[,1],na.rm=TRUE)
for(y_i in 1:y_n){
graphics::segments(eval(parse(text=paste(T))),eval(parse(text=paste(P))),eval(parse(text=paste(T)))+hor_length,eval(parse(text=paste(P))))
}
}

#H2O content at solidus
if(FALSE){
x_intercept<-1100
x_scaling<-70

H2O_at_solidus<-matrix(NA,y_n,1)
for(y_i in 1:y_n){
for(x_i in 1:x_n){
chk<-try(crust[[y_i]][[x_i]]["Melt_rs","mass"],silent=TRUE)
if(class(chk)!="try-error"){
H2O_at_solidus[y_i,1]<-crust[[y_i]][[x_i]]["Bulk_rs","H2O"]
break}
}
}
plot_points<-NULL
for(y_i in 1:y_n){
dTdP<-input_pt[[y_i]][[x_n]][2]/input_pt[[y_i]][[x_n]][1]
plot_y<-x_intercept/(dTdP)
plot_x<-x_intercept+x_scaling*H2O_at_solidus[y_i,1]
plot_points<-rbind(plot_points,matrix(c(plot_x,plot_y),1,2))
}
graphics::plot(plot_points,xlim=c(0,1300),ylim=c(0,30))
lines(plot_points)
}


#Given paths and T intercept draw H2O curve
if(FALSE){
T_x<-1100
gap<-20
H2O_fac<-75
#fixtag should take derivative
#sin((y_i-1+6)*pi/180)*10/140*(x_i-1)/cos((y_i-1+6)*pi/180)*280/140*(x_i-1)
plot_points<-NULL
for(y_i in 1:y_n){
dTdP<-input_pt[[y_i]][[x_n]][2]/input_pt[[y_i]][[x_n]][1]
plot_y<-T_x/(dTdP)
#fixtag - must find point of max H2O change (should be holding on to last change)
plot_x<-T_x+gap+H2O_fac*crust[[y_i]][[431]]["Bulk_rs","H2O"]
plot_points<-rbind(plot_points,matrix(c(plot_x,plot_y),1,2))
}
graphics::plot(plot_points,xlim=c(0,1300),ylim=c(0,30))

#seg_points<-matrix(c(plot_points[,1][-nrow(plot_points)],plot_points[,2][-nrow(plot_points)],plot_points[,1][-1],plot_points[,2][-1]),,4)
#colnames(seg_points)<-c("y0","x0","y1","x1")
#graphics::segments(seg_points[,"x0"],seg_points[,"y0"],seg_points[,"x1"],seg_points[,"y1"])
#need aspect ratio such that degrees equate to C/kbar rounded
#15C/kbar
#tan(theta)=y/x
}

#Points calculated in P-T space
if(FALSE){
cont_mat<-NULL
for(y_i in 1:y_n){
#for(x_i in seq(1,x_n,10)){
for(x_i in 1:x_n){
if(calculation_matrix[y_i,x_i]==1){
cont_mat<-rbind(cont_mat,matrix(c(input_pt[[y_i]][[x_i]][2],input_pt[[y_i]][[x_i]][1]),1,2))
}
}
}
graphics::plot(cont_mat)
graphics::plot(cont_mat,xlim=c(0,1100),ylim=c(0,25))
graphics::plot(cont_mat,xlim=c(600,900),ylim=c(0,25),asp=1)
}


#Plot dT/dP given starting P-T
if(FALSE){
#plot(640,3,xlim=c(650,950),ylim=c(3,20),xaxp=c(650,950,4),yaxp=c(3,20,6))
plot(0,0,xlim=c(0,950),ylim=c(0,20),xaxp=c(0,950,19),yaxp=c(0,20,20))
start_T<-0
start_P<-0
end_T<-950
#dTdP<-c(89.4722222222222,76.25,105.111111111111,80.2222222222222,45.3611111111111,105.888888888889,55.4722222222222,81.6111111111111,58.5833333333333,62.6944444444444,100.5,91.5,77.1944444444444,84.3333333333333,110.166666666667,78.9722222222222,78.8611111111111,83.6666666666667,70.8333333333333,108.305555555556,94.0277777777778,59.2777777777778,69.8333333333333,45.3333333333333,72.3888888888889,48.75,49.2222222222222)
dTdP<-read.csv("AMG.csv")[,5]
end_P<-start_P+(end_T-start_T)/dTdP
graphics::segments(start_T,start_P,end_T,end_P,col=read.csv("AMG.csv")[,1])
}

#Irregular grid contouring
if(FALSE){
#Solidus
#mm<-flip_y(grid_data("wt%","Bio_rs",crust_out())[[2]])
mm<-flip_y(grid_data("Min_formula","Bio_rs",select=5,crust=crust_out())[[2]])
#mm<-flip_y(grid_data("mass","Bi_rs")[[2]])
#melt, pheng, crd wt%
cont_incr<-c(0,1,seq(10,100,10))
cont_incr<-c(0,seq(0.71,1.8,0.05))

cont_incr<-seq(3.1,3.8,0.1)
#cont_incr<-c(0,100)

#remove non numeric
for(y_i in 1:y_n){
for(x_i in 1:x_n){
suppressWarnings(class(mm[y_i,x_i])<-"numeric")
}
}


for(y_i in 1:y_n){
for(x_i in 1:x_n){
if(!is.na(mm[y_i,x_i])){if(mm[y_i,x_i]==0){mm[y_i,x_i]<-NA}
}
}
}

cont_mat<-NULL
#if from left to right cross any contour value then that contour value is assinged left coordinate
#if(min(mm,na.rm=TRUE)<min(cont_incr)){cont_incr<-c(min(mm,na.rm=TRUE),cont_incr)}
#if(max(mm,na.rm=TRUE)>max(cont_incr)){cont_incr<-c(cont_incr,max(mm,na.rm=TRUE))}
for(y_i in 1:y_n){
for(x_i in 1:(x_n-1)){
pnt_1<-mm[y_i,x_i]<=cont_incr
pnt_2<-mm[y_i,x_i+1]<=cont_incr
if(!identical(pnt_1,pnt_2)){
if(is.na(pnt_1!=pnt_2)[1]){rbind(cont_mat,matrix(c(input_pt[[y_i]][[x_i]][2],input_pt[[y_i]][[x_i]][1],0),1,3))}else{
cont_mat<-rbind(cont_mat,matrix(c(input_pt[[y_i]][[x_i]][2],input_pt[[y_i]][[x_i]][1],max(which(pnt_1!=pnt_2))),1,3))
}
}
}
}
grDevices::dev.new()
graphics::plot(cont_mat,col=cont_mat[,3],xlim=c(600,1100),ylim=c(2.5,20))
#contour vals
conts<-sort(intersect(cont_mat[,3],cont_mat[,3]))
draw<-matrix(c(rep(1,length(conts)),1:length(conts)),,2)
plot(draw,col=conts)
labels<-cont_incr[sort(intersect(cont_mat[,3],cont_mat[,3]))]



lev<-seq(3.3,3.8,0.05)
col_num<-0:length(lev)/length(lev)
col_num[length(col_num)]<-0.93
col_use<-rev(grDevices::gray(col_num))
graphics::filled.contour(matrix(cont_incr,,3),levels=lev,col=col_use)


plot(matrix(c(rep(1,length(cont_incr)),1:length(cont_incr)),,2),col=cont_incr)
val_num<-4
x<-1
y<-1:4
val<-cont_incr[1:length(cont_incr)]


outfile_path<-paste0(sub("/code","/Projects",getwd()),"/",working_file,"_melt cont.ps")
postscript(file=outfile_path,onefile=TRUE,horizontal=TRUE)
plot(cont_mat,col=cont_mat[,3],xlim=c(640,920),ylim=c(2,12))
dev.off()


plot(cont_mat,col=cont_mat[,3],axes=FALSE,
plot.axes={
		graphics::axis(1,seq(600,1200,50)) 
		#axis(2,(0:(length(left)-1))/(length(left)-1),left)
}
)


liq_coords<-apply(mm,1, function(x) head(which(!x==0),1))
liq_coords2<-list(x=liq_coords,y=1:length(liq_coords))
plot(liq_coords2)
}

#Irregular grid contouring
if(FALSE){
#Solidus
mm<-flip_y(grid_data("wt%","Melt_rs",crust_out())[[2]])
#cont_incr<-c(1,2,3)
cont_incr<-c(1,seq(10,100,10))
cont_mat<-NULL
#if from left to right cross any contour value then that contour value is assinged left coordinate
if(min(mm,na.rm=TRUE)<min(cont_incr)){cont_incr<-c(min(mm,na.rm=TRUE),cont_incr)}
if(max(mm,na.rm=TRUE)>max(cont_incr)){cont_incr<-c(cont_incr,max(mm,na.rm=TRUE))}
for(y_i in 1:y_n){
for(x_i in 1:(x_n-1)){
pnt_1<-mm[y_i,x_i]<=cont_incr
pnt_2<-mm[y_i,x_i+1]<=cont_incr
if(!identical(pnt_1,pnt_2)){
if(is.na(pnt_1!=pnt_2)[1]){rbind(cont_mat,matrix(c(input_pt[[y_i]][[x_i]][2],input_pt[[y_i]][[x_i]][1],0),1,3))}else{
cont_mat<-rbind(cont_mat,matrix(c(input_pt[[y_i]][[x_i]][2],input_pt[[y_i]][[x_i]][1],max(which(pnt_1!=pnt_2))),1,3))
}
}
}
}
dev.new()
plot(cont_mat,col=cont_mat[,3])


outfile_path<-paste0(sub("/code","/Projects",getwd()),"/",working_file,"_melt cont.ps")
grDevices::postscript(file=outfile_path,onefile=TRUE,horizontal=TRUE)
plot(cont_mat,col=cont_mat[,3],xlim=c(640,920),ylim=c(2,12))
grDevices::dev.off()






plot(cont_mat,col=cont_mat[,3],axes=FALSE,
plot.axes={
		graphics::axis(1,seq(600,1200,50)) 
		#axis(2,(0:(length(left)-1))/(length(left)-1),left)
}
)


liq_coords<-apply(mm,1, function(x) head(which(!x==0),1))
liq_coords2<-list(x=liq_coords,y=1:length(liq_coords))
plot(liq_coords2)
}


#P-T trajectory versus H2O in Bt	
if(FALSE){	
lev<-c(2.95,seq(3,3.5,0.05))
bottom<-seq(1,x_n,length.out=6)
left<-(seq(1,y_n,length.out=6)-1)*3
graphics::filled.contour(t(flip_y(grid_data("H2O","Bt_rs",crust_out())[[2]])),levels=lev,col=rev(gray(0:length(lev)/length(lev))),
plot.axes={
points(lapply(liq_coords2$x,function(i){i/x_n}),lapply(liq_coords2$y,function(i){i/y_n}));
#points(lapply(1:y_n,function(i){i/x_n}),lapply(1:y_n,function(i){i/y_n}));
		if(!is.null(bottom)){axis(1,(0:(length(bottom)-1))/(length(bottom)-1),bottom)}else{NULL}; 
		if(!is.null(left)){axis(2,(0:(length(left)-1))/(length(left)-1),left)}else{NULL};
})

contour(t(flip_y(grid_data("mass","Liq_rs",crust_out())[[2]])),levels=c(0,0.1,1:43),
plot.axes={
		if(!is.null(bottom)){axis(1,(0:(length(bottom)-1))/(length(bottom)-1),bottom)}else{NULL}; 
		if(!is.null(left)){axis(2,(0:(length(left)-1))/(length(left)-1),left)}else{NULL};
})
contour(t(flip_y(grid_data("mass","Liq_rs",crust_out())[[2]])),levels=c(0,0.1,1:43))
lapply(liq_coords2$y,function(i){i/y_n})
#PT paths
plot(unlist(input_pt)[seq(2,length(unlist(input_pt)),2)],unlist(input_pt)[seq(1,length(unlist(input_pt)),2)])


#plot(unlist(input_pt)[seq(2,length(unlist(input_pt)),2)],unlist(input_pt)[seq(1,length(unlist(input_pt)),2)],

#bottom2<-c(1,"",3,"",5,"",7,"")
bottom<-seq(650,1150,50)
left<-c(4,5,6)
plot(c(1,2,3),axes=FALSE,
plot.axes={
		if(!is.null(bottom)){axis(1,(0:(length(bottom)-1))/(length(bottom)-1),bottom)}else{NULL}; 
		if(!is.null(left)){axis(2,(0:(length(left)-1))/(length(left)-1),left)}else{NULL};
}
)

PT_grad<-NULL
for(i in 1:y_n){
#PT_grad<-c(PT_grad,input_pt[[i]][[x_n]][2]-input_pt[[i]][[1]][2]/((input_pt[[i]][[x_n]][1]-input_pt[[i]][[1]][1])*18/5))
PT_grad<-c(PT_grad,(input_pt[[i]][[x_n]][1]-input_pt[[i]][[1]][1])/(input_pt[[i]][[x_n]][2]-input_pt[[i]][[1]][2]))
}

#Solidus
mm<-flip_y(grid_data("mass","Melt_rs",crust_out())[[2]])
liq_coords<-apply(mm,1, function(x) head(which(!x==0),1))
liq_coords2<-list(x=liq_coords,y=1:length(liq_coords))
graphics::plot(liq_coords2)

pt_solidus<-NULL
for(y_i in 1:y_n){
if(length(liq_coords2$x[[y_i]])<=0){
pt_solidus<-c(pt_solidus,0,0)
}else{
pt_solidus<-c(pt_solidus,input_pt[[liq_coords2$y[[y_i]]]][[liq_coords2$x[[y_i]]]])
}
}
graphics::points(pt_solidus[seq(2,length(pt_solidus),2)],pt_solidus[seq(1,length(pt_solidus),2)],col="red")
}

#Ti vs H plot
if(FALSE){
grid_dat<-matrix(0,1,4)
colnames(grid_dat)<-c("H","Ti","path","temp")
for(y_i in seq(1,y_n,10)){
# 20 shows change in slope at high T
#for(y_i in 60:70){
for(x_i in 1:x_n){
grid_dat<-rbind(grid_dat,matrix(c(Bt_formula(y_i,x_i)[6]*2,Bt_formula(y_i,x_i)[12],y_i,input_pt[[y_i]][[x_i]][2]),1,4))
}}

plot(grid_dat,col=grid_dat[,3],xlim=c(0.32,0.39),ylim=c(0.02,0.06),cex=(grid_dat[,4]/100-6))
plot(grid_dat,col=gray(grid_dat[,3]/y_n),xlim=c(0.32,0.39),ylim=c(0.02,0.06),cex=(grid_dat[,4]/100-6))
}

#P-T Grid off of TIN surface     PT_tin(grid_variable,grid_phase,grid_x_n,grid_y_n,grid_temp_def,grid_press_def,crust=crust_out())
PT_tin<-function(grid_variable,grid_phase,grid_x_n,grid_y_n,grid_temp_def,grid_press_def,PT_aspect_ratio=1,return_plot=FALSE,crust=crust_out(),select=5){
#vals<-PT_tin("mass","Melt_rs",5,3,"590+100*grid_x_i","2+5*grid_y_i",1/0.03884347562821677263094156827127,FALSE,crust_out())
#vals<-PT_tin("Min_formula","Bio_rs",5,3,"590+100*grid_x_i","2+5*grid_y_i",1/0.03884347562821677263094156827127,FALSE,crust_out(),select=5)
#vals<-PT_tin("wt%","Bio_rs",5,3,"590+100*grid_x_i","2+5*grid_y_i",1/0.03884347562821677263094156827127,FALSE,crust_out(),select=5)
#vals<-PT_tin("mass","Melt_rs",50,35,"590+10*grid_x_i","2+0.5*grid_y_i",1/0.03884347562821677263094156827127,FALSE,crust_out())
#grid_phase<-"Bio_rs"
#grid_variable<-"Min_formula"
#grid_x_n<-5
#grid_y_n<-3
#grid_temp_def<-"750+2*(grid_x_i-1)"
#grid_press_def<-"2.5+0.07*(grid_y_i-1)"
#crust<-crust_out()
#PT_aspect_ratio<-1/0.03884347562821677263094156827127
#add in aspect ratio

all_temps<-unlist(input_pt)[seq(2,length(unlist(input_pt)),2)]
all_pres<-unlist(input_pt)[seq(1,length(unlist(input_pt)),2)]*PT_aspect_ratio
pnts<-matrix(c(all_temps,all_pres),,2)
vals<-matrix(,grid_y_n,grid_x_n)
for(grid_y_i in 1:grid_y_n){
for(grid_x_i in 1:grid_x_n){
a<-c(eval(parse(text=grid_temp_def)),eval(parse(text=grid_press_def))*PT_aspect_ratio)
#calc distances
dists<-apply(pnts,1,function(x)sqrt(sum((x-a)^2)))
#find closest 3
names(dists)<-1:length(dists)
closest<-sort(dists)[1:3]
value<-NULL
for(i in 1:length(closest)){
y_i<-floor(as.numeric(names(closest)[i])/x_n)
x_i<-as.numeric(names(closest)[i])-y_i*x_n
if(is.null(value)){value<-closest[i]/sum(closest)*num(get_val(y_i,x_i,grid_phase,grid_variable,crust,select=select))
}else{value<-value+closest[i]/sum(closest)*num(get_val(y_i,x_i,grid_phase,grid_variable,crust,select=select))
}
}
vals[grid_y_i,grid_x_i]<-value
cat("point grid_y_i=",grid_y_i,"grid_x_i=",grid_x_i,"finished, ",round(((grid_y_i-1)*grid_x_n+grid_x_i)/(grid_y_n*grid_x_n)*100,2),"% complete\n")
flush.console()
}
}

if(TRUE){
outfile<-paste0(getwd(),"/",working_file,"_",grid_phase,"_",grid_variable,".csv")
write.csv(vals,outfile)
system('shutdown -l')
}

if(return_plot){
bottom<-unlist(lapply((1:grid_x_n),function(grid_x_i){eval(parse(text=grid_temp_def))}))
left<-unlist(lapply((1:grid_y_n),function(grid_y_i){eval(parse(text=grid_press_def))}))
return(
graphics::filled.contour(t(flip_y(vals)),plot.axes={
graphics::axis(1,(0:(length(bottom)-1))/(length(bottom)-1),bottom);
graphics::axis(2,(0:(length(left)-1))/(length(left)-1),left)}
,levels=c(0,1)
)
)
}else{
return(flip_y(vals))
}
}

load_grid<-function(grid_temp_def,grid_press_def){
grid_temp_def<-"600+2*(x_i-1)"
grid_press_def<-"2.5+0.07*(y_i-1)"
outfile<-paste0(gsub("code","Projects/",getwd()),working_file,"/Outputs/",working_file,"_grid.rds")
store_r$crust_r<-readRDS(outfile)
crust<-readRDS(outfile)
x_n<-length(crust[[1]])
y_n<-length(crust)
input_pt<-rep(list(rep(list(c(Pressure=NULL,Temperature=NULL)),x_n)),y_n) 
for(y_i in 1:y_n){
for(x_i in 1:x_n){
input_pt[[y_i]][[x_i]]<-matrix(c(eval(parse(text=grid_press_def)),eval(parse(text=grid_temp_def))),nrow=1,dimnames = list("PT",c("Pressure","Temperature")))
}
}
}
#PT_paths_to_grid(251,180,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",1/0.03884347562821677263094156827127,crust,write_to_file=TRUE,log_off=FALSE)
#P-T paths to Grid
PT_paths_to_grid<-function(grid_x_n,grid_y_n,grid_temp_def,grid_press_def,PT_aspect_ratio=1,crust=crust_out(),write_to_file=FALSE,log_off=FALSE){
crust_grid<-rep(list(rep(list(NULL),grid_x_n)),grid_y_n)
max_d<-matrix(,grid_y_n,grid_x_n)
for(grid_y_i in 1:grid_y_n){
for(grid_x_i in 1:grid_x_n){
#Find bordering 8
first_down<-NULL
first_left<-NULL
first_up<-NULL
first_right<-NULL
#assume x_i increments paths of y_i
#given the line T_c find the closest P_c's to eval parse P
#P_c = P_a + (P_b-P_a)/(T_b-T_a)*(T_c-T_a)
P_cs<-NULL
T_c<-eval(parse(text=grid_temp_def))
for(y_i in 1:y_n){
temp_y_i<-unlist(input_pt[[y_i]])[seq(2,length(unlist(input_pt[[y_i]])),2)]
names(temp_y_i)<-(1:length(temp_y_i))
temp_y_i<-sort(temp_y_i)
pnt_a<-as.numeric(names(which(temp_y_i<T_c)[length(which(temp_y_i<T_c))]))
pnt_b<-as.numeric(names(which(temp_y_i>T_c)[1]))
#error handling
if(is.real(pnt_a,1,x_n)&is.real(pnt_b,1,x_n)){
P_a<-input_pt[[y_i]][[pnt_a]][,"Pressure"]
P_b<-input_pt[[y_i]][[pnt_b]][,"Pressure"]
T_a<-input_pt[[y_i]][[pnt_a]][,"Temperature"]
T_b<-input_pt[[y_i]][[pnt_b]][,"Temperature"]
P_c<-P_a+(P_b-P_a)/(T_b-T_a)*(T_c-T_a)
P_cs<-rbind(P_cs,matrix(c(P_c,y_i,pnt_a,pnt_b),1,4))
}
}
#sort by column 1
if(!is.null(P_cs)){
P_cs<-P_cs[order(P_cs[,1]),]
P_eval<-eval(parse(text=grid_press_def))
first_down<-P_cs[which(P_cs[,1]<P_eval)[length(which(P_cs[,1]<P_eval))],]
first_up<-P_cs[which(P_cs[,1]>P_eval)[1],]
}
#given the line P_c find the closest T_c's to eval parse T
#T_c = T_a + (T_b-T_a)/(P_b-P_a)*(P_c-P_a)
T_cs<-NULL
P_c<-eval(parse(text=grid_press_def))
for(y_i in 1:y_n){
press_y_i<-unlist(input_pt[[y_i]])[seq(1,length(unlist(input_pt[[y_i]])),2)]
names(press_y_i)<-(1:length(press_y_i))
press_y_i<-sort(press_y_i)
#max below
pnt_a<-as.numeric(names(which(press_y_i<=P_c)[length(which(press_y_i<=P_c))]))
#min above
pnt_b<-as.numeric(names(which(press_y_i>P_c)[1]))
#error handling
if(is.real(pnt_a,1,x_n)&is.real(pnt_b,1,x_n)){
P_a<-input_pt[[y_i]][[pnt_a]][,"Pressure"]
P_b<-input_pt[[y_i]][[pnt_b]][,"Pressure"]
T_a<-input_pt[[y_i]][[pnt_a]][,"Temperature"]
T_b<-input_pt[[y_i]][[pnt_b]][,"Temperature"]
T_c<-T_a+(T_b-T_a)/(P_b-P_a)*(P_c-P_a)
T_cs<-rbind(T_cs,matrix(c(T_c,y_i,pnt_a,pnt_b),1,4))
}
}
#sort by column 1
if(!is.null(T_cs)){
T_cs<-T_cs[order(T_cs[,1]),]
T_eval<-eval(parse(text=grid_temp_def))
first_left<-T_cs[which(T_cs[,1]<T_eval)[length(which(T_cs[,1]<T_eval))],]
first_right<-T_cs[which(T_cs[,1]>T_eval)[1],]
}
bordering_8<-NULL
# at least one opposite pair must exist (i.e. you must be within modelled space)
if((is.real(first_left)&is.real(first_right))|(is.real(first_up)&is.real(first_down))){
#compile bordering 8 points
bordering_8<-matrix(,8,4)
if(is.real(first_left)){bordering_8[1,]<-c(first_left[2],first_left[3],input_pt[[first_left[2]]][[first_left[3]]])}
if(is.real(first_left)){bordering_8[2,]<-c(first_left[2],first_left[4],input_pt[[first_left[2]]][[first_left[4]]])}
if(is.real(first_right)){bordering_8[3,]<-c(first_right[2],first_right[3],input_pt[[first_right[2]]][[first_right[3]]])}
if(is.real(first_right)){bordering_8[4,]<-c(first_right[2],first_right[4],input_pt[[first_right[2]]][[first_right[4]]])}
if(is.real(first_up)){bordering_8[5,]<-c(first_up[2],first_up[3],input_pt[[first_up[2]]][[first_up[3]]])}
if(is.real(first_up)){bordering_8[6,]<-c(first_up[2],first_up[4],input_pt[[first_up[2]]][[first_up[4]]])}
if(is.real(first_down)){bordering_8[7,]<-c(first_down[2],first_down[3],input_pt[[first_down[2]]][[first_down[3]]])}
if(is.real(first_down)){bordering_8[8,]<-c(first_down[2],first_down[4],input_pt[[first_down[2]]][[first_down[4]]])}
colnames(bordering_8)<-c("path","point","Pressure","Temperature")
#find closest 2 from 2 different paths given aspect ratio
#apply aspect ratio
asp_P<-matrix(bordering_8[,3]*PT_aspect_ratio,8,1)
colnames(asp_P)<-"aspect P"
bordering_8<-cbind(bordering_8,asp_P)
#calculate distances between c and closest points
pnt_loc<-c(T_eval,P_eval*PT_aspect_ratio)
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
dists<-matrix(,8,1)
colnames(dists)<-"dists"
for(j in 1:8){
dists[j,1]<-euc.dist(pnt_loc,bordering_8[j,4:5])
}
bordering_8<-cbind(bordering_8,dists)
#sort by dists
bordering_8<-bordering_8[order(bordering_8[,6]),]
}
path_a_2pnts<-bordering_8[which(bordering_8[,1]==bordering_8[1,1])[1:2],]
path_b_2pnts<-bordering_8[which(bordering_8[,1]!=bordering_8[1,1])[1:2],]
paths_4pnts<-rbind(path_a_2pnts,path_b_2pnts)
#if unsuccessful to find 2 paths then place blank value
if(any(is.null(paths_4pnts))){
paths_4pnts<-NA
}
if(any(is.na(paths_4pnts))){
max_d[grid_y_i,grid_x_i]<-NA
crust_grid[[grid_y_i]][[grid_x_i]]<-NA
}else{
#normailse proportions based on distances
proportions<-matrix((sum(paths_4pnts[,6])/paths_4pnts[,6])/sum((sum(paths_4pnts[,6])/paths_4pnts[,6])),4,1)
colnames(proportions)<-"proportions"
paths_4pnts<-cbind(paths_4pnts,proportions)
#get list of all phases
all_phases<-NULL
all_columns<-NULL
for(i in 1:4){
all_phases<-union(all_phases,rownames(crust[[paths_4pnts[i,1]]][[paths_4pnts[i,2]]]))
if(nrow(crust[[paths_4pnts[i,1]]][[paths_4pnts[i,2]]])>1){
all_columns<-union(all_columns,colnames(crust[[paths_4pnts[i,1]]][[paths_4pnts[i,2]]]))
}
}
if(is.null(all_columns)){prop_result<-NA}else{
#create proportion result as sum of results*proportions
prop_result<-matrix(0,length(all_phases),length(all_columns))
rownames(prop_result)<-all_phases
colnames(prop_result)<-all_columns
for(i in 1:4){
for(ph in all_phases){
#if phase is present add its proportion to prop_result
chk<-try(crust[[paths_4pnts[i,1]]][[paths_4pnts[i,2]]][ph,],silent=TRUE)
if(class(chk)!="try-error"){
#make sure point has been calculated
if(nrow(crust[[paths_4pnts[i,1]]][[paths_4pnts[i,2]]])>1){
prop_result[ph,]<-prop_result[ph,]+crust[[paths_4pnts[i,1]]][[paths_4pnts[i,2]]][ph,]*proportions[i]
}
}
}
}
}
max_d[grid_y_i,grid_x_i]<-max(path_a_2pnts[,6],path_b_2pnts[,6])
crust_grid[[grid_y_i]][[grid_x_i]]<-prop_result
}
}
cat("point grid_y_i=",grid_y_i,"grid_x_i=",grid_x_i,"finished, ",round(((grid_y_i-1)*grid_x_n+grid_x_i)/(grid_y_n*grid_x_n)*100,2),"% complete\n")
flush.console()
}

if(write_to_file){
#outfile<-paste0(gsub("code","Projects/",getwd()),working_file,"/",working_file,".RData")
#save.image(file=paste0(projects_directory,"/",working_file,"/",working_file,".RData"))
outfile<-paste0(gsub("code","Projects/",getwd()),working_file,"/Outputs/",working_file,"_grid.rds")
saveRDS(crust_grid,outfile)
}
if(log_off){
system('shutdown -l')
}
return(crust_grid)
}


# Contour Grid
if(FALSE){
Grid_variable<-"wt%"
Grid_variable_phase<-"Melt_rs"
crust_in<-crust_out()
grid_x_n<-251
grid_y_n<-251
grid_temp_def<-"600+2*(grid_x_i-1)"
grid_press_def<-"2.5+0.07*(grid_y_i-1)"
x_int<-5
y_int<-10
lower_lim<-NA
upper_lim<-NA
filled<-TRUE
contour_levels<-c(0,0.1,1,10,20,30,40,50,60,70,80,90,100)
select=5
#contour_levels<-c(3.2,3.3,3.4,3.5,3.6,3.7,3.8)
#Crd max 30, min 0, conts at c(0,0.1,5,10,15,20,25,30,35,40,45)
#Ms max 43, min 0, conts at 0,0.1,5,10,15,20,25,30,35,40,45
#melt c(0,0.1,1,10,20,30,40,50,60,70,80,90,100)
#Bio OH c(3.2,3.3,3.4,3.5,3.6,3.7,3.8)
#Bulk H2O seq(0,2.6,0.1)
contour_data("wt%","Gt_rs",151,41,crust_out(),"600+2*(grid_x_i-1)","4+0.1*(grid_y_i-1)",graph="filled",write_to_file=TRUE,file_format=".png")

contour_data("wt%","Bio_rs",251,251,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",filled="yes",write_to_file=TRUE)
contour_data("wt%","Ms_rs",251,251,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=TRUE,file_format=".png",lower_lim=0.001)
contour_data("wt%","Bio_rs",251,251,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",filled="both",write_to_file=FALSE)
contour_data("wt%","Melt_rs",251,251,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",filled=FALSE,write_to_file=FALSE)
contour_data("H2O","Bulk_rs",251,251,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",filled="both",write_to_file=FALSE,contour_levels=seq(0,2.6,0.1))
contour_data("Min_formula","Bio_rs",251,251,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",filled="both",write_to_file=FALSE,contour_levels=c(3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8),select=5)
}

#contour_data("wt%","Ms_rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png")
#contour_data("wt%","Melt_rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps")
contour_all<-function(){
for(Grid_variable_phase in sort(list_all_phases(crust))){
contour_data("wt%",Grid_variable_phase,251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=TRUE,file_format=".png",lower_lim=0.001)
contour_data("wt%",Grid_variable_phase,251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png")
contour_data("wt%",Grid_variable_phase,251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps")
}
contour_data("H2O","Bulk_rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0,2.6,0.2))
contour_data("H2O","Bulk_rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png",contour_levels=seq(0,2.6,0.2))
contour_data("H2O","Bulk_rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",contour_levels=seq(0,2.6,0.2))
contour_data("Min_formula","Bio_rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=TRUE,file_format=".png",contour_levels=c(3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8),select=5)
contour_data("Min_formula","Bio_rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png",contour_levels=c(3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8),select=5)
contour_data("Min_formula","Bio_rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",contour_levels=c(3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8),select=5)
}


if(FALSE){
contour_data("Custom","H2O in Biotite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100")
contour_data("Custom","H2O in Phengite",251,180,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust_in)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust_in)[[2]]/100")
contour_data("Custom","H2O in Cordierite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100")
contour_data("Custom","H2O in hydrous silicates",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100")


contour_data("Read","H2O saturation of Melt",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0,100,10))
contour_data("Read","H2O saturation of Melt",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0,100,10))


contour_data("Custom","H2O in Biotite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="contour",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100")
contour_data("Custom","H2O in Phengite",251,180,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="contour",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust_in)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust_in)[[2]]/100")
contour_data("Custom","H2O in Cordierite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="contour",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100")
contour_data("Custom","H2O in hydrous silicates",251,180,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="contour",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust_in)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust_in)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust_in)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust_in)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust_in)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust_in)[[2]]/100")


contour_data("Custom","H2O in Bulk",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="contour",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"H2O\",\"Bulk_rs\",crust)[[2]]")

contour_data("Custom","Pore fluid lost",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=FALSE,file_format=".png",,contour_levels=seq(14.5,16.3,0.1),custom_def="grid_data(\"mass\",\"H2O_es_cumul\",crust)[[2]]")
max(grid_data("mass","H2O_es_cumul",crust_out())[[2]])


contour_data("Custom","H2O percent Biotite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,5,0.1),custom_def="grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]")
contour_data("Custom","H2O percent in Phengite",251,180,crust_out(),"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,5,0.1),custom_def="grid_data(\"H2O\",\"Ms_rs\",crust_in)[[2]]")
contour_data("Custom","H2O percent in Cordierite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,5,0.1),custom_def="grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]")


 max(grid_data("H2O","Crd_rs",crust)[[2]]*grid_data("wt%","Crd_rs",crust)[[2]]/100)

contour_data("Custom","H2O in Micas",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100")
contour_data("Custom","H2O hydrous silicates",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100")
contour_data("Custom","H2O in Melt",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Liq_rs\",crust)[[2]]*grid_data(\"H2O\",\"Liq_rs\",crust)[[2]]/100")
contour_data("Custom","H2O in Rock",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Liq_rs\",crust)[[2]]*grid_data(\"H2O\",\"Liq_rs\",crust)[[2]]/100")

#idea - graph wt% of silicate as colours and H2O content as intensity, will then see that phengite does not dehydrate it only breaksdown while Bio and Crd dehydrate and then breakdown, when Bio is sole hydrous dehydrates drastically so that it doesn't melt (does it decrease in mode at same time?)
#problem Bulk_rs H2O is not necessarily normalised

(grid_data("wt%","Ms_rs",crust_out())[[2]]*grid_data("H2O","Ms_rs",crust_out())[[2]]/100+grid_data("wt%","Bio_rs",crust_out())[[2]]*grid_data("H2O","Bio_rs",crust_out())[[2]]/100+grid_data("wt%","Crd_rs",crust_out())[[2]]*grid_data("H2O","Crd_rs",crust_out())[[2]]/100+grid_data("wt%","H2O_rs",crust_out())[[2]]+grid_data("wt%","Liq_rs",crust_out())[[2]]*grid_data("H2O","Liq_rs",crust_out())[[2]]/100)/grid_data("H2O","Bulk_rs",crust_out())
graphics::filled.contour(flip_y(t((grid_data("wt%","Ms_rs",crust_out())[[2]]*grid_data("H2O","Ms_rs",crust_out())[[2]]/100+grid_data("wt%","Bio_rs",crust_out())[[2]]*grid_data("H2O","Bio_rs",crust_out())[[2]]/100+grid_data("wt%","Crd_rs",crust_out())[[2]]*grid_data("H2O","Crd_rs",crust_out())[[2]]/100+grid_data("wt%","H2O_rs",crust_out())[[2]]+grid_data("wt%","Liq_rs",crust_out())[[2]]*grid_data("H2O","Liq_rs",crust_out())[[2]]/100)/grid_data("H2O","Bulk_rs",crust_out())[[2]])))


#H2O in silicates+H2O_rs+Liq_rs
graphics::filled.contour(t(flip_y(grid_data("wt%","Ms_rs",crust_out())[[2]]*grid_data("H2O","Ms_rs",crust_out())[[2]]/100)))

#H2O in silicates+H2O_rs+Liq_rs
graphics::filled.contour(t(flip_y(grid_data("wt%","Ms_rs",crust_out())[[2]]*grid_data("H2O","Ms_rs",crust_out())[[2]]/100+grid_data("wt%","Bio_rs",crust_out())[[2]]*grid_data("H2O","Bio_rs",crust_out())[[2]]/100+grid_data("wt%","Crd_rs",crust_out())[[2]]*grid_data("H2O","Crd_rs",crust_out())[[2]]/100+grid_data("wt%","H2O_rs",crust_out())[[2]]+grid_data("wt%","Liq_rs",crust_out())[[2]]*grid_data("H2O","Liq_rs",crust_out())[[2]]/100)))

#H2O prop in silicates+H2O_rs+Liq_rs
graphics::filled.contour(t(flip_y((grid_data("wt%","Ms_rs",crust_out())[[2]]*grid_data("H2O","Ms_rs",crust_out())[[2]]/100+grid_data("wt%","Bio_rs",crust_out())[[2]]*grid_data("H2O","Bio_rs",crust_out())[[2]]/100+grid_data("wt%","Crd_rs",crust_out())[[2]]*grid_data("H2O","Crd_rs",crust_out())[[2]]/100+grid_data("wt%","H2O_rs",crust_out())[[2]]+grid_data("wt%","Liq_rs",crust_out())[[2]]*grid_data("H2O","Liq_rs",crust_out())[[2]]/100)/(grid_data("H2O","Bulk_rs",crust_out())[[2]]/(grid_data("H2O","Bulk_rs",crust_out())[[2]]+grid_data("MGO","Bulk_rs",crust_out())[[2]]+grid_data("AL2O3","Bulk_rs",crust_out())[[2]]+grid_data("SIO2","Bulk_rs",crust_out())[[2]]+grid_data("K2O","Bulk_rs",crust_out())[[2]]+grid_data("CAO","Bulk_rs",crust_out())[[2]]+grid_data("TIO2","Bulk_rs",crust_out())[[2]]+grid_data("FEO","Bulk_rs",crust_out())[[2]]+grid_data("O","Bulk_rs",crust_out())[[2]]+grid_data("NA2O","Bulk_rs",crust_out())[[2]])))))

#H2O in silicates
# (ph_wt%*ph_H2O/100)
contour_data("Custom","H2O hydrous silicates",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100")

#H2O prop in silicates
#(ph_wt%*ph_H2O/100)/(Bulk_H2O/Bulk_comps) - what wt% of the rs is water in biotite
graphics::filled.contour(t(flip_y((grid_data("wt%","Ms_rs",crust_out())[[2]]*grid_data("H2O","Ms_rs",crust_out())[[2]]/100+grid_data("wt%","Bio_rs",crust_out())[[2]]*grid_data("H2O","Bio_rs",crust_out())[[2]]/100+grid_data("wt%","Crd_rs",crust_out())[[2]]*grid_data("H2O","Crd_rs",crust_out())[[2]]/100)/(grid_data("H2O","Bulk_rs",crust_out())[[2]]/(grid_data("H2O","Bulk_rs",crust_out())[[2]]+grid_data("MGO","Bulk_rs",crust_out())[[2]]+grid_data("AL2O3","Bulk_rs",crust_out())[[2]]+grid_data("SIO2","Bulk_rs",crust_out())[[2]]+grid_data("K2O","Bulk_rs",crust_out())[[2]]+grid_data("CAO","Bulk_rs",crust_out())[[2]]+grid_data("TIO2","Bulk_rs",crust_out())[[2]]+grid_data("FEO","Bulk_rs",crust_out())[[2]]+grid_data("O","Bulk_rs",crust_out())[[2]]+grid_data("NA2O","Bulk_rs",crust_out())[[2]])))))
contour_data("Custom","H2O prop in silicates",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=c(0,0.001,1,10,20,30,40,50,60,70,80,90,100),custom_def="(grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100)/(grid_data(\"H2O\",\"Bulk_rs\",crust_out())[[2]]/(grid_data(\"H2O\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"MGO\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"AL2O3\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"SIO2\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"K2O\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"CAO\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"TIO2\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"FEO\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"O\",\"Bulk_rs\",crust_out())[[2]]+grid_data(\"NA2O\",\"Bulk_rs\",crust_out())[[2]]))")
#error on picking up Ms


contour_data("Custom","H2O prop in silicates",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png",,contour_levels=c(0,0.001,1,10,20,30,40,50,60,70,80,90,100),custom_def="(grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100)/(grid_data(\"H2O\",\"Bulk_rs\",crust)[[2]]/(grid_data(\"H2O\",\"Bulk_rs\",crust)[[2]]+grid_data(\"MGO\",\"Bulk_rs\",crust)[[2]]+grid_data(\"AL2O3\",\"Bulk_rs\",crust)[[2]]+grid_data(\"SIO2\",\"Bulk_rs\",crust)[[2]]+grid_data(\"K2O\",\"Bulk_rs\",crust)[[2]]+grid_data(\"CAO\",\"Bulk_rs\",crust)[[2]]+grid_data(\"TIO2\",\"Bulk_rs\",crust)[[2]]+grid_data(\"FEO\",\"Bulk_rs\",crust)[[2]]+grid_data(\"O\",\"Bulk_rs\",crust)[[2]]+grid_data(\"NA2O\",\"Bulk_rs\",crust)[[2]]))")
#H2O prop in melt
graphics::filled.contour(t(flip_y((grid_data("wt%","Liq_rs",crust_out())[[2]]*grid_data("H2O","Liq_rs",crust_out())[[2]]/100)/(grid_data("H2O","Bulk_rs",crust_out())[[2]]/(grid_data("H2O","Bulk_rs",crust_out())[[2]]+grid_data("MGO","Bulk_rs",crust_out())[[2]]+grid_data("AL2O3","Bulk_rs",crust_out())[[2]]+grid_data("SIO2","Bulk_rs",crust_out())[[2]]+grid_data("K2O","Bulk_rs",crust_out())[[2]]+grid_data("CAO","Bulk_rs",crust_out())[[2]]+grid_data("TIO2","Bulk_rs",crust_out())[[2]]+grid_data("FEO","Bulk_rs",crust_out())[[2]]+grid_data("O","Bulk_rs",crust_out())[[2]]+grid_data("NA2O","Bulk_rs",crust_out())[[2]])))))


grid_data("H2O","Bulk_rs",crust)[[2]]*grid_data("H2O","Ms_rs",crust)[[2]]

contour_data("Custom","H2O in rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"H2O_rs\",crust)[[2]]")

contour_data("Custom","H2O in rs",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.2),custom_def="grid_data(\"wt%\",\"Bulk_rs\",crust)[[2]]*grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]*grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100+grid_data(\"wt%\",\"H2O_rs\",crust)[[2]]")

contour_data("Custom","H2O Biotite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(3.0000000000000000001,4,0.1),custom_def="grid_data(\"H2O\",\"Bio_rs\",crust)[[2]]")
contour_data("Custom","H2O Phengite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(3.0000000000000000001,5,0.1),custom_def="grid_data(\"H2O\",\"Ms_rs\",crust)[[2]]")
contour_data("Custom","H2O Cordierite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=FALSE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,2,0.1),custom_def="grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]")
contour_data("Custom","H2O Melt",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=seq(0.0000000000000000001,17,1),custom_def="grid_data(\"H2O\",\"Liq_rs\",crust)[[2]]")

contour_data("Custom","wt Biotite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=TRUE,file_format=".png",,contour_levels=c(0,0.001,1,10,20,30,40,50,60,70,80,90,100),custom_def="grid_data(\"wt%\",\"Bio_rs\",crust)[[2]]")
contour_data("Custom","wt Phengite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=TRUE,file_format=".png",,contour_levels=seq(3.0000000000000000001,5,0.1),custom_def="grid_data(\"wt%\",\"Ms_rs\",crust)[[2]]")
contour_data("Custom","wt Cordierite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2,0.1),custom_def="grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]")
contour_data("Custom","wt Melt",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=TRUE,file_format=".png",,contour_levels=seq(0.0000000000000000001,17,1),custom_def="grid_data(\"wt%\",\"Liq_rs\",crust)[[2]]")


contour_data("Custom","wt Melt",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=c(0.001,1,10,20,30,40,50,60,70,80,90,100),custom_def="grid_data(\"wt%\",\"Melt_rs\",crust)[[2]]")
contour_data("Custom","wt Melt",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="both",write_to_file=TRUE,file_format=".ps",,contour_levels=c(0.001,1,10,20,30,40,50,60,70,80,90,100),custom_def="grid_data(\"mass\",\"Melt_rs\",crust)[[2]]/(100-grid_data(\"mass\",\"H2O_rs\",crust)[[2]])*100")

contour_data("Custom","H2O in Cordierite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=FALSE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.05),custom_def="grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100")

list("H2O in Bio",grid_data("wt%","Bio_rs",crust_out())[[2]]*grid_data("H2O","Bio_rs",crust_out())[[2]]/100)
list("H2O in Ms",grid_data("wt%","Ms_rs",crust_out())[[2]]*grid_data("H2O","Ms_rs",crust_out())[[2]]/100)
list("H2O in Crd",grid_data("wt%","Crd_rs",crust_out())[[2]]*grid_data("H2O","Crd_rs",crust_out())[[2]]/100)
list("H2O in Micas",grid_data("wt%","Ms_rs",crust_out())[[2]]*grid_data("H2O","Ms_rs",crust_out())[[2]]/100+grid_data("wt%","Bio_rs",crust_out())[[2]]*grid_data("H2O","Bio_rs",crust_out())[[2]]/100)
}

list_all_phases<-function(crust=crust_out()){
      all_columns<-NULL
	  all_phases<-NULL
      for(y_i in 1:length(crust)){
        for(x_i in 1:length(crust[[1]])){
		  all_columns<-union(all_columns,colnames(crust[[y_i]][[x_i]]))
          all_phases<-union(all_phases,rownames(crust[[y_i]][[x_i]]))
        }
      }
	  return(all_phases)
}



#what is this even plotting?
lev<-c(3,seq(3.01,4,0.01))
bottom<-seq(1,x_n,length.out=6)
left<-(seq(1,y_n,length.out=6)-1)
col_num<-0:length(lev)/length(lev)
col_num[length(col_num)]<-0.93
col_use<-rev(grDevices::gray(col_num))
#col_use[1]<-gray(0.01)
outfile_path<-paste0(sub("/code","/Projects",getwd()),"/",working_file,"_H.ps")
postscript(file=outfile_path,onefile=TRUE,horizontal=TRUE)
graphics::filled.contour(t(pull_data),levels=lev)
dev.off()

,col=col_use

if(FALSE){
low_lim<-0.1
vals<-fhfa
pres<-matrix(,grid_y_n,grid_x_n)
for(y_i in 1:grid_y_n){
for(x_i in 1:grid_x_n){
if(vals[y_i,x_i]<low_lim){pres[y_i,x_i]<-0}else{pres[y_i,x_i]<-1}
}
}
image(t(flip_y(pres)))

outfile_path<-paste0(getwd(),"/pres.csv")
write.csv(pres,outfile_path,row.names = FALSE)

#Crd max 30, min 0, conts at c(0,0.1,5,10,15,20,25,30,35,40,45)
#Ms max 43, min 0, conts at 0,0.1,5,10,15,20,25,30,35,40,45
#melt c(0,0.1,1,10,20,30,40,50,60,70,80,90,100)
#Bio c(3.2,3.3,3.4,3.5,3.6,3.7,3.8)

grid_x_n<-251
grid_y_n<-251
grid_temp_def<-"600+2*(grid_x_i-1)"
grid_press_def<-"2.5+0.07*(grid_y_i-1)"
file_in<-"PA_650,5_grid_Bulk_rs H2O.csv"
x_int<-5
y_int<-10
bottom<-unlist(lapply((1:grid_x_n),function(grid_x_i){eval(parse(text=grid_temp_def))}))
left<-unlist(lapply((1:grid_y_n),function(grid_y_i){eval(parse(text=grid_press_def))}))
#wet<-as.matrix(read.csv(paste0(getwd(),"/",file_in)))
wet<-as.matrix(read.csv(paste0(projects_directory,"/",working_file,"/Outputs/",file_in)))
#contour(t(flip_y(wet)),levels=c(0,0.1,5,10,15,20,25,30,35,40,45),axes=FALSE,plot.axes={
#contour(t(wet),levels=c(3.2,3.3,3.4,3.5,3.6,3.7,3.8),axes=FALSE,plot.axes={
contour(t(flip_y(wet)),levels=seq(0.5,2.6,0.01),axes=FALSE,plot.axes={
graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
}
)

graphics::filled.contour(t(wet),levels=seq(0.4,2.6,0.2))

low_lim<-0.1
vals<-fhfa
pres<-matrix(,grid_y_n,grid_x_n)
for(y_i in 1:grid_y_n){
for(x_i in 1:grid_x_n){
if(vals[y_i,x_i]<low_lim){pres[y_i,x_i]<-1}
}
}
image(t(flip_y(pres)))


contour(t(flip_y(fhfa)))

fhfa<-as.matrix(read.csv(paste0(getwd(),"/PA_650,5_Bio_rs_Min_formula.csv")))
contour(t(flip_y(fhfa)),levels=c(0,0.1))
graphics::filled.contour(t(flip_y(fhfa)))

}

#Nearest point on path Gridding {npop(grid_variable,grid_phase,grid_start,grid_incr,grid_x_n)}
#Essentially path specific kriging, calculates closest value to each grid increment by using spatial relation and value matrix
npop<-function(grid_variable,grid_phase,grid_start,grid_incr,grid_x_n){
#grid_start<-640
#grid_incr<-2
#grid_x_n<-256	#256
#grid_phase<-"Bt_formula"
#grid_variable<-"12"
grid_mat<-matrix(,y_n,x_n)
grid_val<-NULL
grid_chk<-NULL
all_temps<-unlist(input_pt)[seq(2,length(unlist(input_pt)),2)]
all_pres<-unlist(input_pt)[seq(1,length(unlist(input_pt)),2)]
for(y_i in 1:y_n){
temps<-all_temps[((y_i-1)*grid_x_n+1):(y_i*grid_x_n)]
pres<-all_pres[((y_i-1)*grid_x_n+1):(y_i*grid_x_n)]
for(x_i in 1:x_n){
grid_x_i<-seq(grid_start,by=grid_incr,length.out=grid_x_n)[x_i]
#find if exact temp exists
if(any(temps==grid_x_i)){grid_mat[y_i,x_i]<-get_val(y_i,which(temps==grid_x_i),grid_phase,grid_variable)
grid_chk<-c(grid_chk,grid_x_i)
}else{
#interpolate between first point above and one below
x_b<-temps[which(temps>grid_x_i)[1]]
x_a<-temps[which(temps>grid_x_i)[1]-1]
y_b<-pres[which(temps>grid_x_i)[1]]
y_a<-pres[which(temps>grid_x_i)[1]-1]
if(!any(unlist(lapply(c(x_a,x_b,y_a,y_b),is.na)))){
#get proportion from distance
vec_gradient<-(y_b-y_a)/(x_b-x_a)
grid_y_i<-vec_gradient*(grid_x_i-x_a)+y_a
vec_prop<-sqrt((grid_x_i-x_a)^2+(grid_y_i-y_a)^2)/sqrt((x_b-x_a)^2+(y_b-y_a)^2)
grid_mat[y_i,x_i]<-vec_prop*get_val(y_i,which(temps>grid_x_i)[1],grid_phase,grid_variable)+(1-vec_prop)*get_val(y_i,which(temps>grid_x_i)[1]-1,grid_phase,grid_variable)
grid_chk<-c(grid_chk,vec_prop*x_b+(1-vec_prop)*x_a)
}
}
}
}
#Flip matrix so origin is bottom left
return(flip_y(grid_mat))
}



#Clean ups
#for(y_i in 1:y_n){
#for(x_i in 1:x_n){
#if(is.na(grid_mat[y_i,x_i])){grid_mat[y_i,x_i]<--1}
#}
#}

#for(y_i in 1:y_n){
#for(x_i in 1:x_n){
#if(!is.na(grid_mat[y_i,x_i])){if(grid_mat[y_i,x_i]>0){if(grid_mat[y_i,x_i]<28){grid_mat[y_i,x_i]<-27}}}
#}
#}


#for(y_i in 1:y_n){
#for(x_i in 1:x_n){
#if(!is.na(grid_mat[y_i,x_i])){if(grid_mat[y_i,x_i]<3){grid_mat[y_i,x_i]<-3}}
#}
#}

#contour(t(grid_mat),levels=c(0.01,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90))
#contour(t(grid_mat),levels=c(0.01,seq(2,60,2)))

#lev<-seq(0.5,2,0.1)

#Plotting
if(FALSE){
lev<-c(3,seq(3.01,4,0.01))
bottom<-seq(1,x_n,length.out=6)
left<-(seq(1,y_n,length.out=6)-1)
col_num<-0:length(lev)/length(lev)
col_num[length(col_num)]<-0.93
col_use<-rev(grDevices::gray(col_num))
#col_use[1]<-gray(0.01)
outfile_path<-paste0(sub("/code","/Projects",getwd()),"/",working_file,"_H.ps")
postscript(file=outfile_path,onefile=TRUE,horizontal=TRUE)
graphics::filled.contour(t(grid_mat),levels=lev,col=col_use)
dev.off()

graphics::axis(4,c(0,0.0285714285714286,0.0571428571428571,0.114285714285714,0.228571428571429,0.414285714285714,0.685714285714286,0.942857142857143,1),c("NA",800,400,200,100,50,25,12.5,10))
graphics::axis(2,c(0,0.142857142857143,0.285714285714286,0.428571428571429,0.571428571428572,0.714285714285714,0.857142857142857,1),c(0,10,20,30,40,50,60,70))
graphics::axis(1,c(0,0.01953125,0.1171875,0.21484375,0.3125,0.41015625,0.5078125,0.60546875,1),c(640,650,700,750,800,850,900,950,1))
graphics::filled.contour(t(grid_mat),levels=lev,col=col_use,plot.axes={
graphics::axis(4,c(0,0.0285714285714286,0.0571428571428571,0.114285714285714,0.228571428571429,0.414285714285714,0.685714285714286,0.942857142857143,1),c("NA",800,400,200,100,50,25,12.5,10));
graphics::axis(2,c(0,0.142857142857143,0.285714285714286,0.428571428571429,0.571428571428572,0.714285714285714,0.857142857142857,1),c(0,10,20,30,40,50,60,70));
graphics::axis(1,c(0,0.01953125,0.1171875,0.21484375,0.3125,0.41015625,0.5078125,0.60546875,1),c(640,650,700,750,800,850,900,950,1))
}
)
dev.off()
graphics::filled.contour(t(grid_mat),levels=lev,col=rev(gray(0:length(lev)/length(lev))))
#plot.axes={
#points(lapply(liq_coords2$x,function(i){i/x_n}),lapply(liq_coords2$y,function(i){i/y_n}));
#points(lapply(1:y_n,function(i){i/x_n}),lapply(1:y_n,function(i){i/y_n}));
#		if(!is.null(bottom)){axis(1,(0:(length(bottom)-1))/(length(bottom)-1),bottom)}else{NULL}; 
#		if(!is.null(left)){axis(2,(0:(length(left)-1))/(length(left)-1),left)}else{NULL};
#})
}





#mol% given x_i,y_i and phase
if(FALSE){
mol_percents<-function(y_i,x_i,phase){
#Molar mass in g/mol
Mr_masses<-c("H2O"=18.01528,"MGO"=40.3044,"SIO2"=60.08,"AL2O3"=101.96,"K2O"=94.2,"CAO"=56.0774,"TIO2"=79.866,"FEO"=71.844,"O"=16,"NA2O"=61.9789)
Mr_wts<-crust[[y_i]][[x_i]][phase,comps]
Mr_mols<-vector("numeric")
for(ox in comps){
Mr_mols[ox]<-Mr_wts[ox]/Mr_masses[ox]
}
return(Mr_mols/sum(Mr_mols)*100)
}
}

#Ti independent
#Fe2O3 independet
#O independent
#OH to Fe3+
#balance anion cation
#phengite
#metamorphic facies diagram
#gautier paper AGG - 2016 nicoli
#laurie et al 2013 - moyen and van hunen 2012
#similar to modelled upper surface of slab


#Mineral Formula
Min_formula<-function(y_i,x_i,phase="Bio_rs",oxy_num=24,site_ocup="biotite",crust=crust_out()){
#Given
#		y_i
#		x_i
#		phase (name of phase in crust)
#		oxy_num
#		site_ocup (formula to fill site occupancy)
#Calculate the mineral formula for phase
#DHZ test comp
#Mr_wts<-c("SIO2"=51.63,"AL2O3"=7.39,"FE2O3"=2.5,"FEO"=5.3,"MNO"=0.17,"MGO"=18.09,"CAO"=12.32,"NA2O"=0.61,"H2O"=2.31)
#comps<-names(Mr_wts)
#comps<-c("SIO2","AL2O3","FE2O3","FEO","MNO","MGO","CAO")
#Molar mass in g/mol
#fixtag - add all possibilitites
Mr_masses<-c("H2O"=18.01528,"MGO"=40.3044,"SIO2"=60.08,"AL2O3"=101.96,"K2O"=94.2,"CAO"=56.0774,"TIO2"=79.866,"FEO"=71.844,"FE2O3"=159.7,"O"=16,"O2"=32,"NA2O"=61.9789,"MNO"=70.94)
Mr_oxygen_per_oxide<-c("H2O"=1,"MGO"=1,"SIO2"=2,"AL2O3"=3,"K2O"=1,"CAO"=1,"TIO2"=2,"FEO"=1,"FE2O3"=3,"NA2O"=1,"MNO"=1)
Mr_cation_per_anion<-c("H2O"=2,"MGO"=1,"SIO2"=1/2,"AL2O3"=2/3,"K2O"=2,"CAO"=1,"TIO2"=1/2,"FEO"=1,"FE2O3"=2/3,"NA2O"=2,"MNO"=1)

chk_phs<-try(crust[[y_i]][[x_i]][phase,],silent=TRUE)
#fixtag - find sensible error return
if(class(chk_phs)=="try-error"){return(paste(phase,"not found at y_i=",y_i,"and x_i=",x_i))}
Mr_wts<-crust[[y_i]][[x_i]][phase,comps]
Mr_mols<-vector("numeric")
for(ox in comps){
Mr_mols[ox]<-Mr_wts[ox]/Mr_masses[ox]
}

#Convert any O to Fe2O3
if(any(comps=="O")){
names(Mr_mols)[which(names(Mr_mols)=="O")]<-"FE2O3"
Mr_mols["FEO"]<-Mr_mols["FEO"]-Mr_mols["FE2O3"]*2
}else{
#Convert any O2 to Fe2O3
if(any(comps=="O2")){
Mr_mols["O2"]<-Mr_mols["O2"]*2
names(Mr_mols)[which(names(Mr_mols)=="O2")]<-"FE2O3"
Mr_mols["FEO"]<-Mr_mols["FEO"]-Mr_mols["FE2O3"]*2
}
}

#Calculate atomic proportions of oxygen from each mol (by multiplying by number of oxygen atoms in each oxide)
Mr_oxygens<-vector("numeric")
for(ox in names(Mr_mols)){
Mr_oxygens[ox]<-Mr_mols[ox]*Mr_oxygen_per_oxide[ox]
}

#Recast numbers of anions on basis of number of oxygen atoms by multiplying by oxy_num/sum(Mr_atom_prop))
Mr_anions<-Mr_oxygens*(oxy_num/sum(Mr_oxygens))

#Calculate number of cations in formula  by multiplying by the cation/anion proportion
Mr_cations<-vector("numeric")
for(ox in names(Mr_anions)){
Mr_cations[ox]<-Mr_anions[ox]*Mr_cation_per_anion[ox]
}

#Calculate cation charge by multiplying by oxygens per oxide*2
Mr_cation_charge<-vector("numeric")
for(ox in names(Mr_cations)){
Mr_cation_charge[ox]<-Mr_cations[ox]*Mr_oxygen_per_oxide[ox]*2
}

#Try fill site occupancies
site_ocuppancy<-NULL
if(site_ocup=="muscovite"){
#Fill Biotite formula units
#K2		Al4		[Si6Al2O20]		(OH,F)4
#here we assume tetrehedral sites not filled by Si are occupied by Al and the remaining Al are in the octahedral site
site_ocuppancy<-c("K2"=0,"(Mg,Fe2+)6-4"=0,"(Fe3+,Al,Ti)0-2"=0,"[Si6-5Al2-3O20]"=20,"(OH,F)4"=0)
if(any(names(Mr_cations)=="K2O")){site_ocuppancy[1]<-Mr_cations["K2O"]}
if(any(names(Mr_cations)=="MGO")){site_ocuppancy[2]<-site_ocuppancy[2]+Mr_cations["MGO"]}
if(any(names(Mr_cations)=="FEO")){site_ocuppancy[2]<-site_ocuppancy[2]+Mr_cations["FEO"]}
if(any(names(Mr_cations)=="FE2O3")){site_ocuppancy[3]<-site_ocuppancy[3]+Mr_cations["FE2O3"]}
if(any(names(Mr_cations)=="AL2O3")){site_ocuppancy[3]<-site_ocuppancy[3]+Mr_cations["AL2O3"]-8+Mr_cations["SIO2"]}
if(any(names(Mr_cations)=="TIO2")){site_ocuppancy[3]<-site_ocuppancy[3]+Mr_cations["TIO2"]}
if(any(names(Mr_cations)=="SIO2")){site_ocuppancy[4]<-site_ocuppancy[4]+Mr_cations["SIO2"]}
if(any(names(Mr_cations)=="AL2O3")){site_ocuppancy[4]<-site_ocuppancy[4]+8-Mr_cations["SIO2"]}
if(any(names(Mr_cations)=="H2O")){site_ocuppancy[5]<-site_ocuppancy[5]+Mr_cations["H2O"]}
}
if(site_ocup=="biotite"){
#Fill Biotite formula units
#K2		(Mg,Fe2+)6-4	(Fe3+,Al,Ti)0-2		[Si6-5Al2-3O20]		(OH,F)4
#here we assume tetrehedral sites not filled by Si are occupied by Al and the remaining Al are in the octahedral site
site_ocuppancy<-c("K2"=0,"(Mg,Fe2+)6-4"=0,"(Fe3+,Al,Ti)0-2"=0,"[Si6-5Al2-3O20]"=20,"(OH,F)4"=0)
if(any(names(Mr_cations)=="K2O")){site_ocuppancy[1]<-Mr_cations["K2O"]}
if(any(names(Mr_cations)=="MGO")){site_ocuppancy[2]<-site_ocuppancy[2]+Mr_cations["MGO"]}
if(any(names(Mr_cations)=="FEO")){site_ocuppancy[2]<-site_ocuppancy[2]+Mr_cations["FEO"]}
if(any(names(Mr_cations)=="FE2O3")){site_ocuppancy[3]<-site_ocuppancy[3]+Mr_cations["FE2O3"]}
if(any(names(Mr_cations)=="AL2O3")){site_ocuppancy[3]<-site_ocuppancy[3]+Mr_cations["AL2O3"]-8+Mr_cations["SIO2"]}
if(any(names(Mr_cations)=="TIO2")){site_ocuppancy[3]<-site_ocuppancy[3]+Mr_cations["TIO2"]}
if(any(names(Mr_cations)=="SIO2")){site_ocuppancy[4]<-site_ocuppancy[4]+Mr_cations["SIO2"]}
if(any(names(Mr_cations)=="AL2O3")){site_ocuppancy[4]<-site_ocuppancy[4]+8-Mr_cations["SIO2"]}
if(any(names(Mr_cations)=="H2O")){site_ocuppancy[5]<-site_ocuppancy[5]+Mr_cations["H2O"]}
}

if(site_ocup=="feldspar"){
#Fill feldspar formula units
#	
}

return(c(site_ocuppancy,Mr_cations,"Total Cation Charge"=sum(Mr_cation_charge)))
}


#Plot Ti vs H in PT space
#plot ideal vector from fully hydrated Ti biotite to anhydrous end member, might be simialr to high P-T gradient line
#limit the temperature for Ti H points, check if still span different Ti H ranges
#high P-T gradient only substitutes H2O in Bt for H2O in Ms thus no melting
#inflection in Ti H substituion is due to abscence of Ilm (Ti buffering phase)
# behaviour is due to multiple substiution reactions the dominant of which is Ti for H2O


#3 wt% water = wet solidus
# fix fhfa angle solidus - make new diagram showing this
if(FALSE){
#Plot Bt formula units
Bt_unit<-5
grid_mat<-matrix(,y_n,x_n)
for(y_i in 1:y_n){
for(x_i in 1:x_n){
chk_phs<-try(crust[[y_i]][[x_i]]["Bio(TCC)_rs",],silent=TRUE)
if(class(chk_phs)=="try-error"){grid_mat[y_i,x_i]<-0}else{
grid_mat[y_i,x_i]<-Bt_formula(y_i,x_i)[Bt_unit]
}
}
}
}


#Sean SEAN

#Accepts custom definition of the form grid_data("H2O","Bulk_rs",crust_out())
	# or of the form list("K2O/Na2O",grid_data("K2O","Melt_rs",crust_out())[[2]]/grid_data("Na2O","Melt_rs",crust_out())[[2]])
	# or of the form list("K2O/Na2O",grid_data("K2O","Melt_rs",crust_out())[[2]]/grid_data("Na2O","Melt_rs",crust_out())[[2]])
	#("variable","where it's from",crust_out())
	#list format: (title, content, crust_out())
	#("vol%",
	#grid_data("vol%","Melt_rs",crust_out())	works
	#grid_data("vol%","Gt_rs",crust_out())/(100-grid_data("vol%","Melt_rs",crust_out()))*100	non-numeric argument to binary operator
	#list(grid_data("vol%","Melt_rs",crust_out()),grid_data("vol%","Melt_rs",crust_out())) doesn't work
	#grid_data("vol%","Melt_rs",crust_out())+grid_data("vol%","Gt_rs",crust_out())	non-numeric argument to binary operator
	#list(grid_data("vol%","Melt_rs",crust_out()))	nope
	#c(grid_data("vol%","Melt_rs",crust_out()))	works
	#c(grid_data("vol%","Melt_rs",crust_out()),grid_data("vol%","Gt_rs",crust_out()))	Only plots Melt
	#c(grid_data("vol%","Melt_rs",crust_out())-1)
	
	#list("K2O/Na2O",grid_data("K2O","Melt_rs",crust_out())[[2]]/grid_data("Na2O","Melt_rs",crust_out())[[2]])	NA and Inf output
	#list("K2O/Na2O",grid_data("K2O","Melt_rs",crust_out())/grid_data("Na2O","Melt_rs",crust_out()))	non-numeric argument to binary operator
	#grid_data("K2O","Melt_rs",crust_out())[[2]]/grid_data("Na2O","Melt_rs",crust_out())[[2]]	cannot plot a grid with only 1 row
	#list(grid_data("K2O","Melt_rs",crust_out())[[2]]/grid_data("Na2O","Melt_rs",crust_out())[[2]])	subscript out of bounds
	#list("K2O/Na2O",grid_data("wt%","Melt_rs",crust_out())[[2]]/grid_data("wt%","Melt_rs",crust_out())[[2]])	NA and 1 (1 is correct)
	#list("Gt vol% w/o melt",grid_data("vol%","Melt_rs",crust_out())[[2]]/grid_data("wt%","Melt_rs",crust_out())[[2]])	works
	#list("Gt vol% w/o melt",grid_data("vol%","Gt_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)	 gives output...
	#grid_data("vol%","Gt_rs",crust_out())	works?
	#list("Bulk - melt wt%",100-grid_data("wt%","Melt_rs",crust_out())[[2]]) works!
	# want to get vol% of each major phase normalised to bulk without melt
	#list("Ap vol% normalised to bulk without melt",grid_data("vol%","Ap_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)
	#list("Bio vol% normalised to bulk without melt",grid_data("vol%","Bio_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)
	#list("Kf vol% normalised to bulk without melt",grid_data("vol%","Kf_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)
	#list("q vol% normalised to bulk without melt",grid_data("vol%","q_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)
	#list("Pl vol% normalised to bulk without melt",grid_data("vol%","Pl_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)
	#list("Gt vol% normalised to bulk without melt",grid_data("vol%","Gt_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)
	#list("Crd vol% normalised to bulk without melt",grid_data("vol%","Crd_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)
	#list("Opx vol% normalised to bulk without melt",grid_data("vol%","Opx_rs",crust_out())[[2]]/(100-grid_data("vol%","Melt_rs",crust_out())[[2]])*100)
	# want to get vol% of each major phase normalised to bulk without H2O as a free phase, as this would escape the rock during cooling.
	#phases = and_rs Ap_rs Bio_rs Crd_rs Gt_rs Ilm_rs Kf_rs Melt_rs Mica_rs Opx_rs Pl_rs q_rs ru_rs sill_rs
	#list("and vol% normalised to bulk without free H2O",grid_data("vol%","and_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Ap vol% normalised to bulk without free H2O",grid_data("vol%","Ap_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Bio vol% normalised to bulk without free H2O",grid_data("vol%","Bio_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Crd vol% normalised to bulk without free H2O",grid_data("vol%","Crd_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Gt vol% normalised to bulk without free H2O",grid_data("vol%","Gt_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Ilm vol% normalised to bulk without free H2O",grid_data("vol%","Ilm_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Kf vol% normalised to bulk without free H2O",grid_data("vol%","Kf_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Melt vol% normalised to bulk without free H2O",grid_data("vol%","Melt_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Mica vol% normalised to bulk without free H2O",grid_data("vol%","Mica_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Opx vol% normalised to bulk without free H2O",grid_data("vol%","Opx_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("Pl vol% normalised to bulk without free H2O",grid_data("vol%","Pl_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("q vol% normalised to bulk without free H2O",grid_data("vol%","q_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("ru vol% normalised to bulk without free H2O",grid_data("vol%","ru_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	#list("sill vol% normalised to bulk without free H2O",grid_data("vol%","sill_rs",crust_out())[[2]]/(100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100)
	
	list("La content",grid_data("La","Mnz_rs",crust_out())[[2]]*grid_data("wt%","Mnz_rs",crust_out())[[2]]/grid_data("La","Bulk_rs",crust_out())[[2]])
	
	CAO phases:
	Ap_rs
	Gt_rs
	Kf_rs
	Melt_rs
	Opx_rs
	Pl_1_rs
	Pl_rs
	
	
	list("Ap CAO mass",grid_data("CAO","Ap_rs",crust_out())[[2]]*grid_data("wt%","Ap_rs",crust_out())[[2]]/100)
	list("Gt CAO mass",grid_data("CAO","Gt_rs",crust_out())[[2]]*grid_data("wt%","Gt_rs",crust_out())[[2]]/100)
	list("Kf CAO mass",grid_data("CAO","Kf_rs",crust_out())[[2]]*grid_data("wt%","Kf_rs",crust_out())[[2]]/100)
	list("Melt CAO mass",grid_data("CAO","Melt_rs",crust_out())[[2]]*grid_data("wt%","Melt_rs",crust_out())[[2]]/100)
	list("Opx CAO mass",grid_data("CAO","Opx_rs",crust_out())[[2]]*grid_data("wt%","Opx_rs",crust_out())[[2]]/100)
	list("Pl_1 CAO mass",grid_data("CAO","Pl_1_rs",crust_out())[[2]]*grid_data("wt%","Pl_1_rs",crust_out())[[2]]/100)
	list("Pl CAO mass",grid_data("CAO","Pl_rs",crust_out())[[2]]*grid_data("wt%","Pl_rs",crust_out())[[2]]/100)
	
	Ba Ce Dy Er Eu Gd Hf Ho La Lu Nb Nd Pb Pr Rb Sm Sr Ta Tb Th Tm U V Y Yb Zr
	La Ce Pr Sm Gd Lu Y Sr are compatible in apatite, priority for plots.
	Ap Melt Bio Crd Gt Ilm Kf Opx Pl q sill
if(FALSE){	
	list("La in Ap per weight",grid_data("La","Ap_rs",crust_out())[[2]]/100*grid_data("wt%","Ap_rs",crust_out())[[2]])
	list("La in Melt per weight",grid_data("La","Melt_rs",crust_out())[[2]]/100*grid_data("wt%","Melt_rs",crust_out())[[2]])
	list("La in Bio per weight",grid_data("La","Bio_rs",crust_out())[[2]]/100*grid_data("wt%","Bio_rs",crust_out())[[2]])
	list("La in Crd per weight",grid_data("La","Crd_rs",crust_out())[[2]]/100*grid_data("wt%","Crd_rs",crust_out())[[2]])
	list("La in Gt per weight",grid_data("La","Gt_rs",crust_out())[[2]]/100*grid_data("wt%","Gt_rs",crust_out())[[2]])
	list("La in Ilm per weight",grid_data("La","Ilm_rs",crust_out())[[2]]/100*grid_data("wt%","Ilm_rs",crust_out())[[2]])
	list("La in Kf per weight",grid_data("La","Kf_rs",crust_out())[[2]]/100*grid_data("wt%","Kf_rs",crust_out())[[2]])
	list("La in Opx per weight",grid_data("La","Opx_rs",crust_out())[[2]]/100*grid_data("wt%","Opx_rs",crust_out())[[2]])
	list("La in Pl per weight",grid_data("La","Pl_rs",crust_out())[[2]]/100*grid_data("wt%","Pl_rs",crust_out())[[2]])
	list("La in q per weight",grid_data("La","q_rs",crust_out())[[2]]/100*grid_data("wt%","q_rs",crust_out())[[2]])
	list("La in sill per weight",grid_data("La","sill_rs",crust_out())[[2]]/100*grid_data("wt%","sill_rs",crust_out())[[2]])
	
	list("Ce in Ap per weight",grid_data("Ce","Ap_rs",crust_out())[[2]]/100*grid_data("wt%","Ap_rs",crust_out())[[2]])
	list("Ce in Melt per weight",grid_data("Ce","Melt_rs",crust_out())[[2]]/100*grid_data("wt%","Melt_rs",crust_out())[[2]])
	list("Ce in Bio per weight",grid_data("Ce","Bio_rs",crust_out())[[2]]/100*grid_data("wt%","Bio_rs",crust_out())[[2]])
	list("Ce in Crd per weight",grid_data("Ce","Crd_rs",crust_out())[[2]]/100*grid_data("wt%","Crd_rs",crust_out())[[2]])
	list("Ce in Gt per weight",grid_data("Ce","Gt_rs",crust_out())[[2]]/100*grid_data("wt%","Gt_rs",crust_out())[[2]])
	list("Ce in Ilm per weight",grid_data("Ce","Ilm_rs",crust_out())[[2]]/100*grid_data("wt%","Ilm_rs",crust_out())[[2]])
	list("Ce in Kf per weight",grid_data("Ce","Kf_rs",crust_out())[[2]]/100*grid_data("wt%","Kf_rs",crust_out())[[2]])
	list("Ce in Opx per weight",grid_data("Ce","Opx_rs",crust_out())[[2]]/100*grid_data("wt%","Opx_rs",crust_out())[[2]])
	list("Ce in Pl per weight",grid_data("Ce","Pl_rs",crust_out())[[2]]/100*grid_data("wt%","Pl_rs",crust_out())[[2]])
	list("Ce in q per weight",grid_data("Ce","q_rs",crust_out())[[2]]/100*grid_data("wt%","q_rs",crust_out())[[2]])
	list("Ce in sill per weight",grid_data("Ce","sill_rs",crust_out())[[2]]/100*grid_data("wt%","sill_rs",crust_out())[[2]])
	
	list("Pr in Ap per weight",grid_data("Pr","Ap_rs",crust_out())[[2]]/100*grid_data("wt%","Ap_rs",crust_out())[[2]])
	list("Pr in Melt per weight",grid_data("Pr","Melt_rs",crust_out())[[2]]/100*grid_data("wt%","Melt_rs",crust_out())[[2]])
	list("Pr in Bio per weight",grid_data("Pr","Bio_rs",crust_out())[[2]]/100*grid_data("wt%","Bio_rs",crust_out())[[2]])
	list("Pr in Crd per weight",grid_data("Pr","Crd_rs",crust_out())[[2]]/100*grid_data("wt%","Crd_rs",crust_out())[[2]])
	list("Pr in Gt per weight",grid_data("Pr","Gt_rs",crust_out())[[2]]/100*grid_data("wt%","Gt_rs",crust_out())[[2]])
	list("Pr in Ilm per weight",grid_data("Pr","Ilm_rs",crust_out())[[2]]/100*grid_data("wt%","Ilm_rs",crust_out())[[2]])
	list("Pr in Kf per weight",grid_data("Pr","Kf_rs",crust_out())[[2]]/100*grid_data("wt%","Kf_rs",crust_out())[[2]])
	list("Pr in Opx per weight",grid_data("Pr","Opx_rs",crust_out())[[2]]/100*grid_data("wt%","Opx_rs",crust_out())[[2]])
	list("Pr in Pl per weight",grid_data("Pr","Pl_rs",crust_out())[[2]]/100*grid_data("wt%","Pl_rs",crust_out())[[2]])
	list("Pr in q per weight",grid_data("Pr","q_rs",crust_out())[[2]]/100*grid_data("wt%","q_rs",crust_out())[[2]])
	list("Pr in sill per weight",grid_data("Pr","sill_rs",crust_out())[[2]]/100*grid_data("wt%","sill_rs",crust_out())[[2]])
	
	list("Sm in Ap per weight",grid_data("Sm","Ap_rs",crust_out())[[2]]/100*grid_data("wt%","Ap_rs",crust_out())[[2]])
	list("Sm in Melt per weight",grid_data("Sm","Melt_rs",crust_out())[[2]]/100*grid_data("wt%","Melt_rs",crust_out())[[2]])
	list("Sm in Bio per weight",grid_data("Sm","Bio_rs",crust_out())[[2]]/100*grid_data("wt%","Bio_rs",crust_out())[[2]])
	list("Sm in Crd per weight",grid_data("Sm","Crd_rs",crust_out())[[2]]/100*grid_data("wt%","Crd_rs",crust_out())[[2]])
	list("Sm in Gt per weight",grid_data("Sm","Gt_rs",crust_out())[[2]]/100*grid_data("wt%","Gt_rs",crust_out())[[2]])
	list("Sm in Ilm per weight",grid_data("Sm","Ilm_rs",crust_out())[[2]]/100*grid_data("wt%","Ilm_rs",crust_out())[[2]])
	list("Sm in Kf per weight",grid_data("Sm","Kf_rs",crust_out())[[2]]/100*grid_data("wt%","Kf_rs",crust_out())[[2]])
	list("Sm in Opx per weight",grid_data("Sm","Opx_rs",crust_out())[[2]]/100*grid_data("wt%","Opx_rs",crust_out())[[2]])
	list("Sm in Pl per weight",grid_data("Sm","Pl_rs",crust_out())[[2]]/100*grid_data("wt%","Pl_rs",crust_out())[[2]])
	list("Sm in q per weight",grid_data("Sm","q_rs",crust_out())[[2]]/100*grid_data("wt%","q_rs",crust_out())[[2]])
	list("Sm in sill per weight",grid_data("Sm","sill_rs",crust_out())[[2]]/100*grid_data("wt%","sill_rs",crust_out())[[2]])
	
	list("Gd in Ap per weight",grid_data("Gd","Ap_rs",crust_out())[[2]]/100*grid_data("wt%","Ap_rs",crust_out())[[2]])
	list("Gd in Melt per weight",grid_data("Gd","Melt_rs",crust_out())[[2]]/100*grid_data("wt%","Melt_rs",crust_out())[[2]])
	list("Gd in Bio per weight",grid_data("Gd","Bio_rs",crust_out())[[2]]/100*grid_data("wt%","Bio_rs",crust_out())[[2]])
	list("Gd in Crd per weight",grid_data("Gd","Crd_rs",crust_out())[[2]]/100*grid_data("wt%","Crd_rs",crust_out())[[2]])
	list("Gd in Gt per weight",grid_data("Gd","Gt_rs",crust_out())[[2]]/100*grid_data("wt%","Gt_rs",crust_out())[[2]])
	list("Gd in Ilm per weight",grid_data("Gd","Ilm_rs",crust_out())[[2]]/100*grid_data("wt%","Ilm_rs",crust_out())[[2]])
	list("Gd in Kf per weight",grid_data("Gd","Kf_rs",crust_out())[[2]]/100*grid_data("wt%","Kf_rs",crust_out())[[2]])
	list("Gd in Opx per weight",grid_data("Gd","Opx_rs",crust_out())[[2]]/100*grid_data("wt%","Opx_rs",crust_out())[[2]])
	list("Gd in Pl per weight",grid_data("Gd","Pl_rs",crust_out())[[2]]/100*grid_data("wt%","Pl_rs",crust_out())[[2]])
	list("Gd in q per weight",grid_data("Gd","q_rs",crust_out())[[2]]/100*grid_data("wt%","q_rs",crust_out())[[2]])
	list("Gd in sill per weight",grid_data("Gd","sill_rs",crust_out())[[2]]/100*grid_data("wt%","sill_rs",crust_out())[[2]])
	
	list("Lu in Ap per weight",grid_data("Lu","Ap_rs",crust_out())[[2]]/100*grid_data("wt%","Ap_rs",crust_out())[[2]])
	list("Lu in Melt per weight",grid_data("Lu","Melt_rs",crust_out())[[2]]/100*grid_data("wt%","Melt_rs",crust_out())[[2]])
	list("Lu in Bio per weight",grid_data("Lu","Bio_rs",crust_out())[[2]]/100*grid_data("wt%","Bio_rs",crust_out())[[2]])
	list("Lu in Crd per weight",grid_data("Lu","Crd_rs",crust_out())[[2]]/100*grid_data("wt%","Crd_rs",crust_out())[[2]])
	list("Lu in Gt per weight",grid_data("Lu","Gt_rs",crust_out())[[2]]/100*grid_data("wt%","Gt_rs",crust_out())[[2]])
	list("Lu in Ilm per weight",grid_data("Lu","Ilm_rs",crust_out())[[2]]/100*grid_data("wt%","Ilm_rs",crust_out())[[2]])
	list("Lu in Kf per weight",grid_data("Lu","Kf_rs",crust_out())[[2]]/100*grid_data("wt%","Kf_rs",crust_out())[[2]])
	list("Lu in Opx per weight",grid_data("Lu","Opx_rs",crust_out())[[2]]/100*grid_data("wt%","Opx_rs",crust_out())[[2]])
	list("Lu in Pl per weight",grid_data("Lu","Pl_rs",crust_out())[[2]]/100*grid_data("wt%","Pl_rs",crust_out())[[2]])
	list("Lu in q per weight",grid_data("Lu","q_rs",crust_out())[[2]]/100*grid_data("wt%","q_rs",crust_out())[[2]])
	list("Lu in sill per weight",grid_data("Lu","sill_rs",crust_out())[[2]]/100*grid_data("wt%","sill_rs",crust_out())[[2]])
	
	list("Y in Ap per weight",grid_data("Y","Ap_rs",crust_out())[[2]]/100*grid_data("wt%","Ap_rs",crust_out())[[2]])
	list("Y in Melt per weight",grid_data("Y","Melt_rs",crust_out())[[2]]/100*grid_data("wt%","Melt_rs",crust_out())[[2]])
	list("Y in Bio per weight",grid_data("Y","Bio_rs",crust_out())[[2]]/100*grid_data("wt%","Bio_rs",crust_out())[[2]])
	list("Y in Crd per weight",grid_data("Y","Crd_rs",crust_out())[[2]]/100*grid_data("wt%","Crd_rs",crust_out())[[2]])
	list("Y in Gt per weight",grid_data("Y","Gt_rs",crust_out())[[2]]/100*grid_data("wt%","Gt_rs",crust_out())[[2]])
	list("Y in Ilm per weight",grid_data("Y","Ilm_rs",crust_out())[[2]]/100*grid_data("wt%","Ilm_rs",crust_out())[[2]])
	list("Y in Kf per weight",grid_data("Y","Kf_rs",crust_out())[[2]]/100*grid_data("wt%","Kf_rs",crust_out())[[2]])
	list("Y in Opx per weight",grid_data("Y","Opx_rs",crust_out())[[2]]/100*grid_data("wt%","Opx_rs",crust_out())[[2]])
	list("Y in Pl per weight",grid_data("Y","Pl_rs",crust_out())[[2]]/100*grid_data("wt%","Pl_rs",crust_out())[[2]])
	list("Y in q per weight",grid_data("Y","q_rs",crust_out())[[2]]/100*grid_data("wt%","q_rs",crust_out())[[2]])
	list("Y in sill per weight",grid_data("Y","sill_rs",crust_out())[[2]]/100*grid_data("wt%","sill_rs",crust_out())[[2]])
	
	list("Sr in Ap per weight",grid_data("Sr","Ap_rs",crust_out())[[2]]/100*grid_data("wt%","Ap_rs",crust_out())[[2]])
	list("Sr in Melt per weight",grid_data("Sr","Melt_rs",crust_out())[[2]]/100*grid_data("wt%","Melt_rs",crust_out())[[2]])
	list("Sr in Bio per weight",grid_data("Sr","Bio_rs",crust_out())[[2]]/100*grid_data("wt%","Bio_rs",crust_out())[[2]])
	list("Sr in Crd per weight",grid_data("Sr","Crd_rs",crust_out())[[2]]/100*grid_data("wt%","Crd_rs",crust_out())[[2]])
	list("Sr in Gt per weight",grid_data("Sr","Gt_rs",crust_out())[[2]]/100*grid_data("wt%","Gt_rs",crust_out())[[2]])
	list("Sr in Ilm per weight",grid_data("Sr","Ilm_rs",crust_out())[[2]]/100*grid_data("wt%","Ilm_rs",crust_out())[[2]])
	list("Sr in Kf per weight",grid_data("Sr","Kf_rs",crust_out())[[2]]/100*grid_data("wt%","Kf_rs",crust_out())[[2]])
	list("Sr in Opx per weight",grid_data("Sr","Opx_rs",crust_out())[[2]]/100*grid_data("wt%","Opx_rs",crust_out())[[2]])
	list("Sr in Pl per weight",grid_data("Sr","Pl_rs",crust_out())[[2]]/100*grid_data("wt%","Pl_rs",crust_out())[[2]])
	list("Sr in q per weight",grid_data("Sr","q_rs",crust_out())[[2]]/100*grid_data("wt%","q_rs",crust_out())[[2]])
	list("Sr in sill per weight",grid_data("Sr","sill_rs",crust_out())[[2]]/100*grid_data("wt%","sill_rs",crust_out())[[2]])
}	
	
	#Need (amount or prop of trace)*(wt%)/100/(bulk trace)*100
	#100's cancel out
	#(amount of trace)*(wt%)/(Bulk trace)
	list("Ce in Ap as a function of total Ce",grid_data("Ce","Ap_rs",crust_out())[[2]]*grid_data("wt%","Ap_rs",crust_out())[[2]]/grid_data("Ce","Bulk_rs",crust_out())[[2]])
	
	
	PG6_partition_traces_basic
	PG6_partition_traces_basic_grid_Ce in Ap per weight
	grid_x_n<-96
	grid_y_n<-40
	grid_temp_def<-"600+2*(grid_x_i-1)"
	grid_press_def<-"0.1+0.1*(grid_y_i-1)"
	file_in<-"PG6_partition_traces_basic_grid_Ce in Ap per weight.csv"
	x_int<-5
	y_int<-10
	bottom<-unlist(lapply((1:grid_x_n),function(grid_x_i){eval(parse(text=grid_temp_def))}))
	left<-unlist(lapply((1:grid_y_n),function(grid_y_i){eval(parse(text=grid_press_def))}))
	#wet<-as.matrix(read.csv(paste0(getwd(),"/",file_in)))
	wet<-as.matrix(read.csv(paste0(projects_directory,"/",working_file,"/Outputs/",file_in)))
	#contour(t(flip_y(wet)),levels=c(0,0.1,5,10,15,20,25,30,35,40,45),axes=FALSE,plot.axes={
	#contour(t(wet),levels=c(3.2,3.3,3.4,3.5,3.6,3.7,3.8),axes=FALSE,plot.axes={
	contour(t(flip_y(wet)),levels=seq(0.5,2.6,0.01),axes=FALSE,plot.axes={
	graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
	graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
	}
	)

	graphics::filled.contour(t(wet),levels=seq(0.4,2.6,0.2))
	
	
	#grid plots for trace in phase
	grid_x_n<-96
	grid_y_n<-40
	grid_temp_def<-"600+2*(grid_x_i-1)"
	grid_press_def<-"0.1+0.1*(grid_y_i-1)"
	file_in<-"PG6_partition_traces_basic_grid_Ce in Ap per weight.csv"
	x_int<-5
	y_int<-10
	bottom<-unlist(lapply((1:grid_x_n),function(grid_x_i){eval(parse(text=grid_temp_def))}))
	left<-unlist(lapply((1:grid_y_n),function(grid_y_i){eval(parse(text=grid_press_def))}))
	wet<-as.matrix(read.csv(paste0(projects_directory,"/",working_file,"/Outputs/",file_in)))
	contour(t(flip_y(wet)),levels=seq(0.5,2.6,0.01),axes=FALSE,plot.axes={
	graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
	graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
	}
	)
	graphics::filled.contour(t(wet),levels=c(2,4,6,8,10,12,14))
	dev.off()
	
	
	crust[[1]][[1]][which(rownames(crust[[1]][[1]])=="Pl_rs"),25]
	x = 151, y = 36
	apply(crust,c(36,151),mean)
	
	contour_data("Custom","H2O in Cordierite",251,180,crust,"600+2*(grid_x_i-1)","2.5+0.07*(grid_y_i-1)",graph="hist",write_to_file=FALSE,file_format=".png",,contour_levels=seq(0.0000000000000000001,2.6,0.05),custom_def="grid_data(\"wt%\",\"Crd_rs\",crust)[[2]]*grid_data(\"H2O\",\"Crd_rs\",crust)[[2]]/100")

	contour_data("wt%","Gt_rs",151,41,crust_out(),"600+2*(grid_x_i-1)","4+0.1*(grid_y_i-1)",graph="filled",write_to_file=TRUE,file_format=".png")
	contour_data("wt%","Gt_rs",151,41,crust_out(),"600+2*(grid_x_i-1)","4+0.1*(grid_y_i-1)",graph="filled",write_to_file=TRUE,file_format=".ps")
	contour_data("wt%","Gt_rs",151,41,crust_out(),"600+2*(grid_x_i-1)","4+0.1*(grid_y_i-1)",graph="filled",write_to_file=TRUE,file_format=".png")
	pull_data<-grid_data("wt%","Gt_rs",crust_in=crust_out(),select=select)[[2]]
	graphics::filled.contour(t(flip_y(pull_data)),levels=contour_levels,axes=FALSE)
	
	#weighted addition
	bulk <- c(67.864,0.439,14.312,2.993,0.049,0.868,1.326,3.032,4.914,4,0.205,100)
	names(bulk) <- c("SIO2","TIO2","Al2O3","FEO","MNO","MGO","CAO","NA2O","K2O","H2O","P2O5","mass")
	PlCore <- c(60.015,0.007,24.711,0,0,0,6.955,6.955,0.83,0,0.097,-9.09)
	names(PlCore) <- c("SIO2","TIO2","Al2O3","FEO","MNO","MGO","CAO","NA2O","K2O","H2O","P2O5","mass")
	effBulk <- .wtd.add(rbind(bulk,PlCore))
	
	values<-effBulk
	total<-100
	override<-FALSE
	normalise<-function(values,total=100,override=FALSE){
		if(override == TRUE){
			#does not include mass, normalises to total.
			if(any(match(colnames(values),"mass"))){
				massPos <- which(colnames(values)=="mass")
				values[,1:(massPos-1)] <- values[,1:(massPos-1)]/sum(values[,-massPos])*total
				values[,"mass"]<-total
			}
		} else if(any(match(colnames(values),"mass"))){
			total<-values[,"mass"]
			massPos <- which(colnames(values)=="mass")
			#values<-values[,-massPos]
			#values<-values/sum(values)*total
			values[,1:(massPos-1)] <- values[,1:(massPos-1)]/sum(values[,-massPos])*total
			#values<-c(values[,1:(massPos-1)],total,values[,(massPos+1):length(values)])
		} else {
			values<-values/sum(values)*total
		}
		return(values)
	}
	

contour_data<-function(Grid_variable,Grid_variable_phase="Bulk_rs",grid_x_n=251,grid_y_n=251,crust_in=crust_out(),grid_temp_def="600+2*(grid_x_i-1)",grid_press_def="2.5+0.07*(grid_y_i-1)",x_int=5,y_int=10,lower_lim=NA,upper_lim=NA,graph="filled",contour_levels=c(0,0.001,1,10,20,30,40,50,60,70,80,90,100),select=1,write_to_file=TRUE,file_format=".ps",custom_def=""){
	bottom<-unlist(lapply((1:grid_x_n),function(grid_x_i){eval(parse(text=grid_temp_def))}))
	left<-unlist(lapply((1:grid_y_n),function(grid_y_i){eval(parse(text=grid_press_def))}))
	if(Grid_variable=="Read"){
		pull_data<-as.matrix(read.csv("H2Omelt.csv"))
	}
	if(Grid_variable=="Custom"){
		pull_data<-eval(parse(text=custom_def))
	}
	if(Grid_variable!="Read"&Grid_variable!="Custom"){
		pull_data<-grid_data(Grid_variable,Grid_variable_phase,crust_in,select=select)[[2]]
	}
	if(!is.na(lower_lim)){
		for(y_i in 1:grid_y_n){
			for(x_i in 1:grid_x_n){
				if(!is.na(pull_data[y_i,x_i])){if(pull_data[y_i,x_i]<lower_lim){pull_data[y_i,x_i]<-NA}}
			}
		}
	}
	if(!is.na(upper_lim)){
		for(y_i in 1:grid_y_n){
			for(x_i in 1:grid_x_n){
				if(!is.na(pull_data[y_i,x_i])){if(pull_data[y_i,x_i]>upper_lim){pull_data[y_i,x_i]<-NA}}
			}
		}
	}
	if(all(is.na(pull_data))){cat(paste(Grid_variable,Grid_variable_phase,"\nis null for all values in selected space\n"))}else{
		file_name<-gsub("%","",Grid_variable)
		if(write_to_file){
			outfile<-paste0(projects_directory,"/",working_file,"/Outputs/",working_file,"_",Grid_variable_phase,"_",file_name,"_",graph,file_format)
			if(file_format==".ps"){
				grDevices::postscript(file=outfile,onefile=TRUE,horizontal=TRUE)
			}
			if(file_format==".png"){
				grDevices::png(file=outfile)
			}
		}

		if(graph=="filled"){
			#Colour options = gray.colors,heat.colours,terrain.colours,rainbow,topo.colours		
			graphics::filled.contour(t(flip_y(pull_data)),levels=contour_levels,axes=FALSE,plot.axes={
				graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
				graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
				},color.palette=eval(parse(text="gray.colors"))
			)
		}
		if(graph=="contour"){
			contour(t(flip_y(pull_data)),levels=contour_levels,axes=FALSE,plot.axes={
				graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
				graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
				}
			)
		}
		if(graph=="both"){
			z<-flip_y(t(pull_data))
			x<-1:nrow(z)
			y<-1:ncol(z)
			levels=contour_levels
			draw.contour <- function()
			{
			  contour(x=x, y=y, z=z, add=TRUE, 
					  levels=levels,
					  drawlabels=FALSE,
					  xlim=rev(range(x)),
					  ylim=rev(range(y)))
			}
			plot(NA,xlim=rev(range(x)),
						   ylim=rev(range(y)),xlab="Temperature",ylab="Pressure",
						   frame=FALSE,axes=F,xaxs="i",yaxs="i")
			graphics::.filled.contour(x=x, y=y, z=z,
						   levels=levels,
						   col=gplots::colorpanel(length(levels) + 1, "white", "grey10"))
			draw.contour()
			graphics::axis(1, c(1,round(length(x)/2),length(x)),label=c(bottom[length(x)],bottom[round(length(x)/2)],bottom[1]), tcl=-0.5)
			graphics::axis(2, c(1,round(length(y)/2),length(y)),label=c(left[length(y)],left[round(length(y)/2)],left[1]), tcl=-0.5)
			graphics::box()
		}
		if(graph=="hist"){
			hist(pull_data,main=Grid_variable_phase,xlab=Grid_variable)
		}
		if(write_to_file){
			grDevices::dev.off()
		}
	}
}



contour_data<-function(Grid_variable,Grid_variable_phase,grid_x_n=96,grid_y_n=40,crust_in=crust_out(),grid_temp_def="600+2*(grid_x_i-1)",grid_press_def="0.1+0.1*(grid_y_i-1)",x_int=5,y_int=9,lower_lim=NA,upper_lim=NA,graph="filled",contour_levels=c(0,1,2),select=1,write_to_file=TRUE,file_format=".png",custom_def=""){
	bottom<-unlist(lapply((1:grid_x_n),function(grid_x_i){eval(parse(text=grid_temp_def))}))
	left<-unlist(lapply((1:grid_y_n),function(grid_y_i){eval(parse(text=grid_press_def))}))
	if(Grid_variable=="Read"){
		pull_data<-as.matrix(read.csv("H2Omelt.csv"))
	}
	if(Grid_variable=="Custom"){
		pull_data<-eval(parse(text=custom_def))
	}
	if(Grid_variable!="Read"&Grid_variable!="Custom"){
		pull_data<-grid_data(Grid_variable,Grid_variable_phase,crust_in,select=select)[[2]]
	}
	if(!is.na(lower_lim)){
		for(y_i in 1:grid_y_n){
			for(x_i in 1:grid_x_n){
				if(!is.na(pull_data[y_i,x_i])){if(pull_data[y_i,x_i]<lower_lim){pull_data[y_i,x_i]<-NA}}
			}
		}
	}
	if(!is.na(upper_lim)){
		for(y_i in 1:grid_y_n){
			for(x_i in 1:grid_x_n){
				if(!is.na(pull_data[y_i,x_i])){if(pull_data[y_i,x_i]>upper_lim){pull_data[y_i,x_i]<-NA}}
			}
		}
	}
	if(all(is.na(pull_data))){cat(paste(Grid_variable,Grid_variable_phase,"\nis null for all values in selected space\n"))}else{
		file_name<-gsub("%","",Grid_variable)
		if(write_to_file){
			outfile<-paste0(projects_directory,"/",working_file,"/Outputs/",working_file,"_",Grid_variable_phase,"_",file_name,"_",graph,file_format)
			if(file_format==".ps"){
				grDevices::postscript(file=outfile,onefile=TRUE,horizontal=TRUE)
			}
			if(file_format==".png"){
				grDevices::png(file=outfile)
			}
		}

		if(graph=="filled"){
			#Colour options = gray.colors,heat.colours,terrain.colours,rainbow,topo.colours		
			graphics::filled.contour(t(flip_y(pull_data)),levels=contour_levels,axes=FALSE,plot.axes={
				graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
				graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
				},col=gplots::colorpanel(length(levels) + 1, "white", "red")
			)
		}
		#color.palette=eval(parse(text="heat.colors"))
		if(graph=="contour"){
			contour(t(flip_y(pull_data)),levels=contour_levels,axes=FALSE,plot.axes={
				graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
				graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
				}
			)
		}
		if(graph=="both"){
			z<-flip_y(t(pull_data))
			x<-1:nrow(z)
			y<-1:ncol(z)
			levels=contour_levels
			draw.contour <- function()
			{
			  contour(x=x, y=y, z=z, add=TRUE, 
					  levels=levels,
					  drawlabels=FALSE,
					  xlim=rev(range(x)),
					  ylim=rev(range(y)))
			}
			plot(NA,xlim=rev(range(x)),
						   ylim=rev(range(y)),xlab="Temperature",ylab="Pressure",
						   frame=FALSE,axes=F,xaxs="i",yaxs="i")
			graphics::.filled.contour(x=x, y=y, z=z,
						   levels=levels,
						   col=gplots::colorpanel(length(levels) + 1, "white", "grey10"))
			draw.contour()
			graphics::axis(1, c(1,round(length(x)/2),length(x)),label=c(bottom[length(x)],bottom[round(length(x)/2)],bottom[1]), tcl=-0.5)
			graphics::axis(2, c(1,round(length(y)/2),length(y)),label=c(left[length(y)],left[round(length(y)/2)],left[1]), tcl=-0.5)
			graphics::box()
		}
		if(graph=="hist"){
			hist(pull_data,main=Grid_variable_phase,xlab=Grid_variable)
		}
		if(write_to_file){
			grDevices::dev.off()
		}
	}
}



contour_data("wt%","Gt_rs",96,40,crust_out(),"600+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",graph="filled",write_to_file=FALSE,file_format=".png")

contour_data("Custom","Ce in Ap per weight",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",graph="filled",write_to_file=FALSE,,contour_levels=c(0,2,4,8,10,12,14),x_int=7,upper_lim=850,custom_def="grid_data(\"Ce\",\"Ap_rs\",crust)[[2]]*grid_data(\"wt%\",\"Ap_rs\",crust)[[2]]/100")

#levels go from 0 to 13 in intervals of 1.
contour_data("Custom","Ce in Ap per weight",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",graph="filled",write_to_file=FALSE,,contour_levels=seq(0,13,1),x_int=7,upper_lim=850,custom_def="grid_data(\"Ce\",\"Ap_rs\",crust)[[2]]*grid_data(\"wt%\",\"Ap_rs\",crust)[[2]]/100")

#levels go from 0.001 to 13 in intervals of 1. This is how to exclude 0
contour_data("Custom",,"Ce in Ap per weight",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,graph="filled",write_to_file=FALSE,,contour_levels=seq(2,13,1),custom_def="grid_data(\"Ce\",\"Ap_rs\",crust)[[2]]*grid_data(\"wt%\",\"Ap_rs\",crust)[[2]]/100")

Ap Melt Bio Crd Gt Ilm Kf Opx Pl q sill
#Next, get trace in element as a function of total trace (bulk)
#Need (amount or prop of trace)*(wt%)/100/(bulk trace)*100
	#100's cancel out
	#(amount of trace)*(wt%)/(Bulk trace)
	list("Ce in Ap as a function of total Ce",grid_data("Ce","Ap_rs",crust_out())[[2]]*grid_data("wt%","Ap_rs",crust_out())[[2]]/grid_data("Ce","Bulk_rs",crust_out())[[2]])
	list("Dy in Melt as a function of total Dy",grid_data("Dy","Melt_rs",crust_out())[[2]]*grid_data("wt%","Melt_rs",crust_out())[[2]]/grid_data("Dy","Bulk_rs",crust_out())[[2]])
#Ap
contour_data("Custom",,"Ce in Ap as a function of total Ce",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,
	graph="filled",write_to_file=FALSE,,contour_levels=seq(0.00000000000001,30,5),
	custom_def="grid_data(\"Ce\",\"Ap_rs\",crust)[[2]]*grid_data(\"wt%\",\"Ap_rs\",crust)[[2]]/grid_data(\"Ce\",\"Bulk_rs\",crust)[[2]]")
#Melt
contour_data("Custom",,"Ce in Melt as a function of total Ce",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,
	graph="filled",write_to_file=FALSE,,contour_levels=seq(0.00000000000001,30,5),
	custom_def="grid_data(\"Ce\",\"Melt_rs\",crust)[[2]]*grid_data(\"wt%\",\"Melt_rs\",crust)[[2]]/grid_data(\"Ce\",\"Bulk_rs\",crust)[[2]]")

#histogram
contour_data("Custom",,"Ce in Ap as a function of total Ce",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,graph="hist",write_to_file=FALSE,,contour_levels=seq(0.00000000000001,50,5),custom_def="grid_data(\"Ce\",\"Ap_rs\",crust)[[2]]*grid_data(\"wt%\",\"Ap_rs\",crust)[[2]]/grid_data(\"Ce\",\"Bulk_rs\",crust)[[2]]")

#contour
contour_data("Custom",,"Ce in Ap as a function of total Ce",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,graph="contour",write_to_file=FALSE,,contour_levels=seq(0.00000000000001,50,5),custom_def="grid_data(\"Ce\",\"Ap_rs\",crust)[[2]]*grid_data(\"wt%\",\"Ap_rs\",crust)[[2]]/grid_data(\"Ce\",\"Bulk_rs\",crust)[[2]]")

#both
contour_data("Custom",,"Ce in Ap as a function of total Ce",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,graph="both",write_to_file=FALSE,,contour_levels=seq(0.00000000000001,50,5),custom_def="grid_data(\"Ce\",\"Ap_rs\",crust)[[2]]*grid_data(\"wt%\",\"Ap_rs\",crust)[[2]]/grid_data(\"Ce\",\"Bulk_rs\",crust)[[2]]")

###########################################################
#######Saving (2D?) matrix to file, reading in again#######
#Useful to compare data from 2 different runs
library(MASS)
write.matrix(x,'file.txt',sep = "\t")
#reading
x2 <- as.matrix(read.table("file.txt", as.is = TRUE))
# make mat2 a true matrix
colnames(x2) <- NULL
mat2 <- unname(x2)
#all.equal(x, x2)
# [1] TRUE
###########################################################




#Contains more information on customising plots
contour_data<-function(Grid_variable,Grid_variable_phase="Bulk_rs",plot_title="Title",grid_x_n=251,grid_y_n=251,crust_in=crust_out(),grid_temp_def="600+2*(grid_x_i-1)",grid_press_def="2.5+0.07*(grid_y_i-1)",x_int=5,y_int=10,lower_lim=NA,upper_lim=NA,graph="filled",contour_levels=c(0,0.001,1,10,20,30,40,50,60,70,80,90,100),select=1,write_to_file=TRUE,file_format=".ps",custom_def="",colorName="YlOrRd",xlab = "Temperature (C)",ylab = "Pressure (kbar)"){
	bottom<-unlist(lapply((1:grid_x_n),function(grid_x_i){eval(parse(text=grid_temp_def))}))
	left<-unlist(lapply((1:grid_y_n),function(grid_y_i){eval(parse(text=grid_press_def))}))
	#values for A/CNK
	mwAl2O3<-101.96128; mwCaO<-56.0774; mwNa2O<-61.97894; mwK2O<-94.19600
	if(Grid_variable=="Read"){
		pull_data<-as.matrix(read.csv("H2Omelt.csv"))
	}
	if(Grid_variable=="Custom"){
		pull_data<-eval(parse(text=custom_def))
	}
	if(Grid_variable!="Read"&Grid_variable!="Custom"){
		pull_data<-grid_data(Grid_variable,Grid_variable_phase,crust_in,select=select)[[2]]
	}
	if(!is.na(lower_lim)){
		for(y_i in 1:grid_y_n){
			for(x_i in 1:grid_x_n){
				if(!is.na(pull_data[y_i,x_i])){if(pull_data[y_i,x_i]<lower_lim){pull_data[y_i,x_i]<-NA}}
			}
		}
	}
	if(!is.na(upper_lim)){
		for(y_i in 1:grid_y_n){
			for(x_i in 1:grid_x_n){
				if(!is.na(pull_data[y_i,x_i])){if(pull_data[y_i,x_i]>upper_lim){pull_data[y_i,x_i]<-NA}}
			}
		}
	}
	if(all(is.na(pull_data))){cat(paste(Grid_variable,Grid_variable_phase,"\nis null for all values in selected space\n"))}else{
		file_name<-gsub("%","",Grid_variable)
		#error may arise when trying to use % in title of saved file. replace with %% instead.
		if(length(grep("%",plot_title))!=0){plot_title<-paste0(strsplit(plot_title,"%")[[1]][1],"%%",strsplit(plot_title,"%")[[1]][2])}else{plot_title1<-plot_title}
		#plot_title<-paste0(strsplit(plot_title,"%")[[1]][1],"%%",strsplit(plot_title,"%")[[1]][2])
		if(length(grep("%",file_name))!=0){file_name<-paste0(strsplit(Grid_variable,"%")[[1]][1],"%%",strsplit(Grid_variable,"%")[[1]][2])}
		#file_name<-paste0(strsplit(Grid_variable,"%")[[1]][1],"%%",strsplit(Grid_variable,"%")[[1]][2])
		if(write_to_file){
			#outfile<-paste0(projects_directory,"/",working_file,"/Outputs/",working_file,"_",Grid_variable_phase,"_",file_name,"_",graph,file_format)
			outfile<-paste0(projects_directory,"/",working_file,"/Outputs/",working_file,"_",file_name,"_",graph,"_",plot_title,file_format)
			cat(outfile)
			if(file_format==".ps"){
				grDevices::postscript(file=outfile,onefile=TRUE,horizontal=TRUE)
			}
			if(file_format==".png"){
				grDevices::png(file=outfile)
			}
		}
		#custom alteration of values (difference between output of 2 different runs
		#pull_data <- pull_data - PG6_ApSat_Bea_Pl_wt
		#pull_data <- pull_data - PG6_ApSat_Bea_Pl_CAO
		# pull_data <- pull_data - PG6_ApSat_Bea_Melt_wt
		#pull_data <- pull_data - PG6_ApSat_Bea_Melt_CAO
		######
		if(graph=="filled"){
			#https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/filled.contour
			#Colour options = gray.colors,heat.colours,terrain.colours,rainbow,topo.colours
			#https://developer.r-project.org/Blog/public/2019/04/01/hcl-based-color-palettes-in-grdevices/
			#https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html
			#or lookup hcl.colors function
			#col_brew<-RColorBrewer::brewer.pal(length(contour_levels)-1,"YlOrRd")
			require("grDevices")
			graphics::filled.contour(t(flip_y(pull_data)),
				levels = contour_levels,
				axes = FALSE,
				#change axis and graph titles here
				plot.title = graphics::title(main = plot_title,
					xlab,
					ylab),
				plot.axes = {
					graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
					graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
				},
				#color.palette=eval(parse(text="heat.colors"))
				#color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE),
				#col_brew<-RColorBrewer::brewer.pal(n=col_levels-1,name=input$Grid_RColorBrewer_colours)
				#col_pal<-colorRampPalette(col_brew)(col_levels)
				#Should add in a reverse palette option:
				col = colorRampPalette(brewer.pal(length(contour_levels)-1,colorName))(length(contour_levels)),
				# col = rev(colorRampPalette(brewer.pal(length(contour_levels)-1,colorName))(length(contour_levels))),
				#col = colorRampPalette(col_brew)(length(contour_levels)),
				key.axes = graphics::axis(4, round(contour_levels,1))
			)
		}
		if(graph=="contour"){
			contour(t(flip_y(pull_data)),levels=contour_levels,axes=FALSE,
				plot.title = graphics::title(main = plot_title,
					xlab, 
					ylab),
				plot.axes={
				graphics::axis(1,((0:(length(bottom)-1))/(length(bottom)-1))[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])],bottom[c(1,seq(0,length(bottom),length(bottom)/x_int)[-1])]);
				graphics::axis(2,((0:(length(left)-1))/(length(left)-1))[c(1,seq(0,length(left),length(left)/y_int)[-1])],left[c(1,seq(0,length(left),length(left)/y_int)[-1])])
				}
			)
		}
		if(graph=="both"){
			z<-flip_y(t(pull_data))
			x<-1:nrow(z)
			y<-1:ncol(z)
			levels=contour_levels
			draw.contour <- function()
			{
			  contour(x=x, y=y, z=z, add=TRUE, 
					  levels=levels,
					  drawlabels=FALSE,
					  xlim=rev(range(x)),
					  ylim=rev(range(y)))
			}
			plot(NA,xlim=rev(range(x)),
						   ylim=rev(range(y)),xlab="Temperature",ylab="Pressure",
						   frame=FALSE,axes=F,xaxs="i",yaxs="i")
			graphics::.filled.contour(x=x, y=y, z=z,
						   levels=levels,
						   col=gplots::colorpanel(length(levels) + 1, "white", "grey10"))
			draw.contour()
			graphics::axis(1, c(1,round(length(x)/2),length(x)),label=c(bottom[length(x)],bottom[round(length(x)/2)],bottom[1]), tcl=-0.5)
			graphics::axis(2, c(1,round(length(y)/2),length(y)),label=c(left[length(y)],left[round(length(y)/2)],left[1]), tcl=-0.5)
			graphics::box()
		}
		if(graph=="hist"){
			hist(pull_data,main=plot_title,xlab=Grid_variable)
		}
		if(write_to_file){
			grDevices::dev.off()
		}
	}
}
list_all_phases<-function(crust=crust_out()){
      all_columns<-NULL
	  all_phases<-NULL
      for(y_i in 1:length(crust)){
        for(x_i in 1:length(crust[[1]])){
		  all_columns<-union(all_columns,colnames(crust[[y_i]][[x_i]]))
          all_phases<-union(all_phases,rownames(crust[[y_i]][[x_i]]))
        }
      }
	  return(all_phases)
}

if(MGO in biotite){
Grid_variable_phase <- "Bio_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Bio_rs",
	plot_title=paste0("wt%_",Grid_variable_phase),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(8.00000000000001,10,0.2),
	select=1,
	write_to_file=FALSE,file_format=".png",
	custom_def=paste0("grid_data(\"MGO\",\"",Grid_variable_phase="Bio_rs","\",crust)[[2]]"),
	colorName="YlOrRd",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")
}

Bio_rs
Crd_rs
Gt_rs
Melt_rs
Opx_rs
Pl_rs
list("Pl_rs wt% change",grid_data("wt%","Pl_rs",crust_out())[[2]]-grid_data("wt%","Pl_rs",crust_out())[[2]][,27])
list("Pl_rs wt% change",function(){if((grid_data("wt%","Pl_rs",crust_out())[[2]]-grid_data("wt%","Pl_rs",crust_out())[[2]][,27]) <= -5){return(0)}else{return(paste0("grid_data(\"wt%\",\"Pl_rs\",crust_out())[[2]]-grid_data(\"wt%\",\"Pl_rs\",crust_out())[[2]][,27]"))}})
list("Pl_rs wt% change",if(grid_data("wt%","Pl_rs",crust_out())[[2]]-grid_data("wt%","Pl_rs",crust_out())[[2]][,27] <= -5){grid_data("wt%","Pl_rs",crust_out())[[2]]=0}else{grid_data("wt%","Pl_rs",crust_out())[[2]]-grid_data("wt%","Pl_rs",crust_out())[[2]][,27]})
if(grid_data("wt%","Pl_rs",crust_out())[[2]]-grid_data("wt%","Pl_rs",crust_out())[[2]][,27] <= -5){grid_data("wt%","Pl_rs",crust_out())[[2]]=0}else{grid_data("wt%","Pl_rs",crust_out())[[2]]-grid_data("wt%","Pl_rs",crust_out())[[2]][,27]}
function(){if(grid_data("wt%","Pl_rs",crust_out())[[2]]-grid_data("wt%","Pl_rs",crust_out())[[2]][,27] <= -5){return(NA)}else{return(grid_data("wt%","Pl_rs",crust_out())[[2]]-grid_data("wt%","Pl_rs",crust_out())[[2]][,27])}}

list("Ap_rs P2O5 portion",grid_data("P2O5","Ap_rs",crust_out())[[2]]*grid_data("mass","Ap_rs",crust_out())[[2]]/grid_data("P2O5","Bulk_rs",crust_out())[[2]])
list("Mnz_rs P2O5 portion",grid_data("P2O5","Mnz_rs",crust_out())[[2]]*grid_data("mass","Mnz_rs",crust_out())[[2]]/grid_data("P2O5","Bulk_rs",crust_out())[[2]])
list("Melt_rs P2O5 portion",grid_data("P2O5","Melt_rs",crust_out())[[2]]*grid_data("mass","Melt_rs",crust_out())[[2]]/grid_data("P2O5","Bulk_rs",crust_out())[[2]])

list("total P2O5 portion",((grid_data("P2O5","Melt_rs",crust_out())[[2]]*grid_data("mass","Melt_rs",crust_out())[[2]])+(grid_data("P2O5","Mnz_rs",crust_out())[[2]]*grid_data("mass","Mnz_rs",crust_out())[[2]])+(grid_data("P2O5","Ap_rs",crust_out())[[2]]*grid_data("mass","Ap_rs",crust_out())[[2]]))/100)
list("total P2O5 portion",((grid_data("P2O5","Melt_rs",crust_out())[[2]]*grid_data("mass","Melt_rs",crust_out())[[2]])+(grid_data("P2O5","Ap_rs",crust_out())[[2]]*grid_data("mass","Ap_rs",crust_out())[[2]]))/100)
list("total CAO portion",((grid_data("CAO","Melt_rs",crust_out())[[2]]*grid_data("mass","Melt_rs",crust_out())[[2]])+(grid_data("CAO","Mnz_rs",crust_out())[[2]]*grid_data("mass","Mnz_rs",crust_out())[[2]])+(grid_data("CAO","Ap_rs",crust_out())[[2]]*grid_data("mass","Ap_rs",crust_out())[[2]]))/100)

#For PG6 T_CAO run
if(PG6 T_CAO){
#Trying to show change in melt wt%, relative to having full CAO in bulk (not accounting for Ca-in-Ap)
Grid_variable_phase<-"Melt_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase,
	plot_title=paste0("wt%_",Grid_variable_phase),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=27,grid_y_n=151,
	crust_in=crust_out(),
	grid_temp_def="(grid_x_i-1)*0.010385+1.056",grid_press_def="600+(grid_y_i-1)*2",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,5,0.5),
	select=1,
	write_to_file=FALSE,file_format=".png",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]-grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]][,27]"),
	colorName="GnBu",
	xlab="CaO wt% Bulk_rs",
	ylab="Temperature (C)")
#Melt wt%
contour_data(Grid_variable="Custom",
	Grid_variable_phase,
	plot_title=paste0("wt%_",Grid_variable_phase),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=27,grid_y_n=151,
	crust_in=crust_out(),
	grid_temp_def="(grid_x_i-1)*0.010385+1.056",grid_press_def="600+(grid_y_i-1)*2",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,100,10),
	select=1,
	write_to_file=TRUE,file_format=".png",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="GnBu",
	xlab="CaO wt% Bulk_rs",
	ylab="Temperature (C)")
#Pl wt%
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("wt%_",Grid_variable_phase="Pl_rs"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=27,grid_y_n=101,
	crust_in=crust_out(),
	grid_temp_def="(grid_x_i-1)*0.009118+1.152927",grid_press_def="650+(grid_y_i-1)*2",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,30,1),
	select=1,
	write_to_file=TRUE,file_format=".ps",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase="Pl_rs","\",crust)[[2]]"),
	colorName="PuBu",
	xlab="CaO wt% Bulk_rs",
	ylab="Temperature (C)")
#Pl wt% change
Grid_variable_phase<-"Pl_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase,
	plot_title=paste0("wt% change in ",Grid_variable_phase),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=27,grid_y_n=101,
	crust_in=crust_out(),
	grid_temp_def="(grid_x_i-1)*0.009118+1.152927",grid_press_def="650+(grid_y_i-1)*2",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(-3.00000000000001,0.5,0.5),
	select=1,
	write_to_file=TRUE,file_format=".ps",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]-grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]][,27]"),
	colorName="YlOrRd",
	xlab="CaO wt% Bulk_rs",
	ylab="Temperature (C)")
#Pl CAO wt%
Grid_variable_phase="Pl_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("",Grid_variable_phase," CAO wt.%"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=27,grid_y_n=101,
	crust_in=crust_out(),
	grid_temp_def="(grid_x_i-1)*0.009118+1.152927",grid_press_def="650+(grid_y_i-1)*2",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,15,1),
	select=1,
	write_to_file=FALSE,file_format=".ps",
	custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="PuBu",
	xlab="CaO wt% Bulk_rs",
	ylab="Temperature (C)")
#Pl CAO wt% change
Grid_variable_phase="Pl_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("CAO wt.% change in ",Grid_variable_phase),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=27,grid_y_n=101,
	crust_in=crust_out(),
	grid_temp_def="(grid_x_i-1)*0.009118+1.152927",grid_press_def="650+(grid_y_i-1)*2",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(-3.00000000000001,0.5,0.5),
	select=1,
	write_to_file=FALSE,file_format=".ps",
	custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]-grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]][,27]"),
	colorName="YlOrRd",
	xlab="CaO wt% Bulk_rs",
	ylab="Temperature (C)")
#Pl CAO mass
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("",Grid_variable_phase," CAO mass"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=27,grid_y_n=151,
	crust_in=crust_out(),
	grid_temp_def="(grid_x_i-1)*0.010385+1.056",grid_press_def="600+(grid_y_i-1)*2",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,2,0.2),
	select=1,
	write_to_file=TRUE,file_format=".ps",
	custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/100"),
	colorName="PuBu",
	xlab="CaO wt% Bulk_rs",
	ylab="Temperature (C)")
#Pl CAO mass change
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("",Grid_variable_phase," CAO mass change"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=27,grid_y_n=151,
	crust_in=crust_out(),
	grid_temp_def="(grid_x_i-1)*0.010385+1.056",grid_press_def="600+(grid_y_i-1)*2",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(-0.50000000000001,0.1,0.05),
	select=1,
	write_to_file=TRUE,file_format=".png",
	custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/100-grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]][,27]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]][,27]/100"),
	colorName="PuBu",
	xlab="CaO wt% Bulk_rs",
	ylab="Temperature (C)")
#For each phase, show change in mode relative to full CAO in bulk
for(Grid_variable_phase in sort(list_all_phases(crust))){
	contour_data(Grid_variable="Custom",
		Grid_variable_phase,
		#plot_title=paste0("CAO in ",Grid_variable_phase),
		plot_title=paste0("CAO change in ",Grid_variable_phase,""),
		grid_x_n=27,grid_y_n=151,
		crust_in=crust_out(),
		grid_temp_def="(grid_x_i-1)*0.010385+1.056",grid_press_def="600+(grid_y_i-1)*2",
		x_int=6,y_int=7,
		lower_lim=NA,upper_lim=NA,
		graph="filled",
		contour_levels=seq(-3.00000000000001,3,0.5),
		select=1,
		write_to_file=TRUE,file_format=".ps",
		custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]-grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]][,27]"),
		colorName="RdBu",
		xlab="CaO wt% Bulk_rs",
		ylab="Temperature (C)")
}
}

eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]")))
eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]-grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]][,27]")))

#Comparing PG6 with apatite and without apatite
if(PG6 Ap comparing with and without apatite){
if(Pl){
#Pl wt%
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("wt%_",Grid_variable_phase="Pl_rs"," range 5"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,30,5),
	select=1,
	write_to_file=TRUE,file_format=".png",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase="Pl_rs","\",crust)[[2]]"),
	colorName="PuBu",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")

PG6_ApSat_Bea_Pl_wt <- eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase="Pl_rs","\",crust)[[2]]")))
#with Ap saving as PG6_ApSat_Bea_Pl_wt
library(MASS)
write.matrix(PG6_ApSat_Bea_Pl_wt,'PG6_ApSat_Bea_Pl_wt.txt',sep = "\t")
#reading in console without Ap
PG6_ApSat_Bea_Pl_wt <- as.matrix(read.table("PG6_ApSat_Bea_Pl_wt.txt", as.is = TRUE))
colnames(PG6_ApSat_Bea_Pl_wt) <- NULL
#might be unnecessary?:
mat2 <- unname(PG6_ApSat_Bea_Pl_wt)
all.equal(PG6_ApSat_Bea_Pl_wt, eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase="Pl_rs","\",crust)[[2]]"))))

#Pl wt% difference (noAp - Ap)
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("wt% ",Grid_variable_phase="Pl_rs"," difference 0-6"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,6,0.5),
	select=1,
	write_to_file=FALSE,file_format=".ps",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase="Pl_rs","\",crust)[[2]]"),
	colorName="PuBu",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")

#Pl CAO wt%
Grid_variable_phase <- "Pl_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("",Grid_variable_phase," CaO wt%"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,16,2),
	select=1,
	write_to_file=TRUE,file_format=".png",
	custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="PuBu",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")

PG6_ApSat_Bea_Pl_CAO <- eval(parse(text=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]")))
#with Ap saving as PG6_ApSat_Bea_Pl_CAO
library(MASS)
write.matrix(PG6_ApSat_Bea_Pl_CAO,'PG6_ApSat_Bea_Pl_CAO.txt',sep = "\t")
#reading in console without Ap
PG6_ApSat_Bea_Pl_CAO <- as.matrix(read.table("PG6_ApSat_Bea_Pl_CAO.txt", as.is = TRUE))
colnames(PG6_ApSat_Bea_Pl_CAO) <- NULL
#might be unnecessary?:
mat2 <- unname(PG6_ApSat_Bea_Pl_CAO)
all.equal(PG6_ApSat_Bea_Pl_CAO, eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase="Pl_rs","\",crust)[[2]]"))))

#Pl CAO difference (noAp - Ap)
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Pl_rs",
	plot_title=paste0("CAO ",Grid_variable_phase," difference -0.4 to 1.8"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(-0.40000000000001,1.8,0.2),
	select=1,
	write_to_file=TRUE,file_format=".png",
	custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="PuBu",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")
}
if(Melt){
#Melt wt%
Grid_variable_phase <- "Melt_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Melt_rs",
	plot_title=paste0("wt%_",Grid_variable_phase),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,100,10),
	select=1,
	write_to_file=TRUE,file_format=".ps",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="PuBu",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")

PG6_ApSat_Bea_Melt_wt <- eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]")))
#with Ap saving as PG6_ApSat_Bea_Melt_wt
library(MASS)
write.matrix(PG6_ApSat_Bea_Melt_wt,'PG6_ApSat_Bea_Melt_wt.txt',sep = "\t")
#reading in console without Ap
PG6_ApSat_Bea_Melt_wt <- as.matrix(read.table("PG6_ApSat_Bea_Melt_wt.txt", as.is = TRUE))
colnames(PG6_ApSat_Bea_Melt_wt) <- NULL
#might be unnecessary?:
PG6_ApSat_Bea_Melt_wt <- unname(PG6_ApSat_Bea_Melt_wt)
all.equal(PG6_ApSat_Bea_Melt_wt, eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"))))

#Melt wt.% difference
Grid_variable_phase <- "Melt_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Melt_rs",
	plot_title=paste0("",Grid_variable_phase," wt% difference -5 to 0"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(-5.00000000000001,0,0.5),
	select=1,
	write_to_file=TRUE,file_format=".ps",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="YlOrRd",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")

#Melt CAO
Grid_variable_phase <- "Melt_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Melt_rs",
	plot_title=paste0("",Grid_variable_phase," CAO wt%"),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,1.2,0.1),
	select=1,
	write_to_file=TRUE,file_format=".png",
	custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="PuBu",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")


PG6_ApSat_Bea_Melt_CAO <- eval(parse(text=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]")))
#with Ap saving as PG6_ApSat_Bea_Melt_CAO
library(MASS)
write.matrix(PG6_ApSat_Bea_Melt_CAO,'PG6_ApSat_Bea_Melt_CAO.txt',sep = "\t")
#reading in console without Ap
PG6_ApSat_Bea_Melt_CAO <- as.matrix(read.table("PG6_ApSat_Bea_Melt_CAO.txt", as.is = TRUE))
colnames(PG6_ApSat_Bea_Melt_CAO) <- NULL
#might be unnecessary?:
PG6_ApSat_Bea_Melt_CAO <- unname(PG6_ApSat_Bea_Melt_CAO)
all.equal(PG6_ApSat_Bea_Melt_CAO, eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"))))

#Melt CaO wt% difference. (No-Ap minus Ap)
Grid_variable_phase <- "Melt_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Melt_rs",
	plot_title=paste0("",Grid_variable_phase," CAO difference -0.2 to 0.2"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(-0.2000000000001,0.2,0.02),
	select=1,
	write_to_file=TRUE,file_format=".png",
	custom_def=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="RdYlGn",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")

eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"))) - PG6_ApSat_Bea_Melt_CAO

eval(parse(text=paste0("grid_data(\"Pressure\",\"Bulk_rs\",crust)[[2]]")))
write.matrix(eval(parse(text=paste0("grid_data(\"Pressure\",\"Bulk_rs\",crust)[[2]]"))),'PG6_noAp_Pressure.txt',sep = "\t")
write.matrix(eval(parse(text=paste0("grid_data(\"Temperature\",\"Bulk_rs\",crust)[[2]]"))),'PG6_noAp_Pressure.txt',sep = "\t")
#PG6_noAp
#Melt CaO wt%
Grid_variable_phase<-"Melt_rs"
PG6_noAp_Melt_CAO <- eval(parse(text=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]")))
library(MASS)
write.matrix(PG6_noAp_Melt_CAO,'PG6_noAp_Melt_CAO.txt',sep = "\t")

#Melt CaO diff
Grid_variable_phase<-"Melt_rs"
PG6_noAp_Melt_CAO_diff <- eval(parse(text=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]"))) - PG6_ApSat_Bea_Melt_CAO
library(MASS)
write.matrix(PG6_noAp_Melt_CAO_diff,'PG6_noAp_Melt_CAO_diff.txt',sep = "\t")

#Pl wt%
Grid_variable_phase<-"Pl_rs"
PG6_noAp_Pl_wt <- eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]")))
library(MASS)
write.matrix(PG6_noAp_Pl_wt,'PG6_noAp_Pl_wt.txt',sep = "\t")

#Melt wt%
Grid_variable_phase<-"Melt_rs"
PG6_noAp_Melt_wt <- eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]")))
library(MASS)
write.matrix(PG6_noAp_Melt_wt,'PG6_noAp_Melt_wt.txt',sep = "\t")

#Pl wt% diff
Grid_variable_phase<-"Pl_rs"
PG6_noAp_Pl_wt_diff <- eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"))) - PG6_ApSat_Bea_Pl_wt
library(MASS)
write.matrix(PG6_noAp_Pl_wt_diff,'PG6_noAp_Pl_wt_diff.txt',sep = "\t")

#Pl CaO wt%
Grid_variable_phase<-"Pl_rs"
PG6_noAp_Pl_CAO <- eval(parse(text=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]")))
library(MASS)
write.matrix(PG6_noAp_Pl_CAO,'PG6_noAp_Pl_CAO.txt',sep = "\t")

#Melt wt% diff
Grid_variable_phase<-"Melt_rs"
PG6_noAp_Melt_wt_diff <- eval(parse(text=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"))) - PG6_ApSat_Bea_Melt_wt
library(MASS)
write.matrix(PG6_noAp_Melt_wt_diff,'PG6_noAp_Melt_wt_diff.txt',sep = "\t")

#Pl CAO diff
Grid_variable_phase<-"Pl_rs"
PG6_noAp_Pl_CAO_diff <- eval(parse(text=paste0("grid_data(\"CAO\",\"",Grid_variable_phase,"\",crust)[[2]]"))) - PG6_ApSat_Bea_Pl_CAO
library(MASS)
write.matrix(PG6_noAp_Pl_CAO_diff,'PG6_noAp_Pl_CAO_diff.txt',sep = "\t")

}
}
if(Ap){
#Ap wt%
Grid_variable_phase <- "Ap_rs"
contour_data(Grid_variable="Custom",
	Grid_variable_phase="Ap_rs",
	plot_title=paste0("wt%_",Grid_variable_phase),
	#plot_title=paste0("Change in wt% Melt, relative to full CaO in bulk composition"),
	grid_x_n=151,grid_y_n=36,
	crust_in=crust_out(),
	grid_temp_def="600+(grid_x_i-1)*2",grid_press_def="0.5+(grid_y_i-1)*0.1",
	x_int=6,y_int=7,
	lower_lim=NA,upper_lim=NA,
	graph="filled",
	contour_levels=seq(0.00000000000001,0.5,0.05),
	select=1,
	write_to_file=TRUE,file_format=".png",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]"),
	colorName="PuBu",
	xlab="Temperature (C)",
	ylab="Pressure (kbar)")
}


Grid_variable_phase<-"Ap_rs"
contour_data("Custom",Grid_variable_phase,paste0("wt%_",Grid_variable_phase),151,36,crust_out(),
	"(grid_x_i-1)*2+600","0.5+0.1*(grid_y_i-1)",x_int=6,y_int=6,graph="filled",write_to_file=TRUE,,
	contour_levels=seq(0.00000000000001,0.6,0.05),file_format=".ps",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/(100-grid_data(\"wt%\",\"","H2O_rs","\",crust)[[2]])*100"))
contour_data("Custom",Grid_variable_phase,paste0("wt%_",Grid_variable_phase),151,36,crust_out(),
	"(grid_x_i-1)*2+600","0.5+0.1*(grid_y_i-1)",x_int=6,y_int=6,graph="filled",write_to_file=TRUE,,
	contour_levels=seq(0.00000000000001,0.6,0.05),file_format=".png",
	custom_def=paste0("grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/(100-grid_data(\"wt%\",\"","H2O_rs","\",crust)[[2]])*100"))
contour_data("Custom",Grid_variable_phase,paste0("vol%_",Grid_variable_phase),151,36,crust_out(),
	"(grid_x_i-1)*2+600","0.5+0.1*(grid_y_i-1)",x_int=6,y_int=6,graph="filled",write_to_file=TRUE,,
	contour_levels=seq(0.00000000000001,0.6,0.05),file_format=".ps",
	custom_def=paste0("grid_data(\"vol%\",\"",Grid_variable_phase,"\",crust)[[2]]/(100-grid_data(\"vol%\",\"","H2O_rs","\",crust)[[2]])*100"))
contour_data("Custom",Grid_variable_phase,paste0("vol%_",Grid_variable_phase),151,36,crust_out(),
	"(grid_x_i-1)*2+600","0.5+0.1*(grid_y_i-1)",x_int=6,y_int=6,graph="filled",write_to_file=TRUE,,
	contour_levels=seq(0.00000000000001,0.6,0.05),file_format=".png",
	custom_def=paste0("grid_data(\"vol%\",\"",Grid_variable_phase,"\",crust)[[2]]/(100-grid_data(\"vol%\",\"","H2O_rs","\",crust)[[2]])*100"))


#plots for Ce
for(Grid_variable_phase in sort(list_all_phases(crust))){
	contour_data("Custom",,paste0("Ce in ",Grid_variable_phase," as a function of total Ce"),96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,
	graph="filled",write_to_file=TRUE,,contour_levels=seq(0.00000000000001,50,5),file_format=".png",
	custom_def=paste0("grid_data(\"Ce\",\"",Grid_variable_phase,"\",crust)[[2]]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/grid_data(\"Ce\",\"Bulk_rs\",crust)[[2]]"))
	contour_data("Custom",,"Ce in Melt as a function of total Ce",96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,
	graph="filled",write_to_file=TRUE,,contour_levels=seq(0.00000000000001,100,5),file_format=".png",
	custom_def="grid_data(\"Ce\",\"Melt_rs\",crust)[[2]]*grid_data(\"wt%\",\"Melt_rs\",crust)[[2]]/grid_data(\"Ce\",\"Bulk_rs\",crust)[[2]]")
}
#plots for Eu
for(Grid_variable_phase in sort(list_all_phases(crust))){
	contour_data("Custom",,paste0("Eu in ",Grid_variable_phase," as a function of total Eu"),96,36,crust_out(),"660+2*(grid_x_i-1)","0.5+0.1*(grid_y_i-1)",7,10,
	graph="filled",write_to_file=TRUE,,contour_levels=seq(0.00000000000001,16,1),file_format=".png",
	custom_def=paste0("grid_data(\"Eu\",\"",Grid_variable_phase,"\",crust)[[2]]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/grid_data(\"Eu\",\"Bulk_rs\",crust)[[2]]"))
	# contour_data("Custom",,paste0("Eu in ",Grid_variable_phase," as a function of total Eu"),96,36,crust_out(),"660+2*(grid_x_i-1)","0.5+0.1*(grid_y_i-1)",7,10,
	# graph="filled",write_to_file=TRUE,,contour_levels=seq(0.00000000000001,100,5),file_format=".png",
	# custom_def=paste0("grid_data(\"Eu\",\"",Grid_variable_phase,"\",crust)[[2]]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/grid_data(\"Eu\",\"Bulk_rs\",crust)[[2]]"))
}
Ba Ce Dy Er Eu Gd Hf Ho La Lu Nb Nd Pb Pr Rb Sm Sr Ta Tb Th Tm U V Y Yb Zr
La Ce Pr Sm Gd Lu Y Sr
traces <- c("Ba","Ce","Dy","Er","Eu","Gd","Hf","Ho","La","Lu","Nb","Nd","Pb","Pr","Rb","Sm","Sr","Ta","Tb","Th","Tm","U","V","Y","Yb","Zr")
#histogram plots for all traces of all phases
for(trc in traces){
for(Grid_variable_phase in sort(list_all_phases(crust))){
	contour_data("Custom",,paste0(trc," in ",Grid_variable_phase," as a function of total ",trc),96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,
	graph="hist",write_to_file=TRUE,,contour_levels=seq(0.00000000000001,20,1),file_format=".png",
	custom_def=paste0("grid_data(\"",trc,"\",\"",Grid_variable_phase,"\",crust)[[2]]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/grid_data(\"",trc,"\",\"Bulk_rs\",crust)[[2]]"))
	contour_data("Custom",,paste0(trc," in Melt_rs as a function of total ",trc),96,40,crust_out(),"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,
	graph="hist",write_to_file=TRUE,,contour_levels=seq(0.00000000000001,20,1),file_format=".png",
	custom_def=paste0("grid_data(\"",trc,"\",\"Melt_rs\",crust)[[2]]*grid_data(\"wt%\",\"Melt_rs\",crust)[[2]]/grid_data(\"",trc,"\",\"Bulk_rs\",crust)[[2]]"))
}
}
#plots for Zr
for(Grid_variable_phase in sort(list_all_phases(crust))){
	contour_data("Custom",,paste0("Zr in ",Grid_variable_phase," as a function of total Zr"),96,40,crust_out(),
		"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,graph="filled",write_to_file=TRUE,,
		contour_levels=seq(0.00000000000001,14,2),file_format=".ps",
		custom_def=paste0("grid_data(\"Zr\",\"",Grid_variable_phase,"\",crust)[[2]]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/grid_data(\"Zr\",\"Bulk_rs\",crust)[[2]]"))
	contour_data("Custom",,paste0("Zr in ",Grid_variable_phase," as a function of total Zr"),96,40,crust_out(),
		"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,graph="filled",write_to_file=TRUE,,
		contour_levels=seq(0.00000000000001,14,2),file_format=".png",
		custom_def=paste0("grid_data(\"Zr\",\"",Grid_variable_phase,"\",crust)[[2]]*grid_data(\"wt%\",\"",Grid_variable_phase,"\",crust)[[2]]/grid_data(\"Zr\",\"Bulk_rs\",crust)[[2]]"))
}

traces <- c("Ba","Ce","Dy","Er","Eu","Gd","Hf","Ho","La","Lu","Nb","Nd","Pb","Pr","Rb","Sm","Sr","Ta","Tb","Th","Tm","U","V","Y","Yb","Zr")
#display melt up to 100% filled plots with PuBu colour
for(trc in traces){
	contour_data("Custom",,paste0(trc," in Melt_rs as a function of total ",trc),96,40,crust_out(),
		"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,graph="filled",write_to_file=TRUE,,
		contour_levels=seq(0.00000000000001,100,5),file_format=".ps",
		custom_def=paste0("grid_data(\"",trc,"\",\"Melt_rs\",crust)[[2]]*grid_data(\"wt%\",\"Melt_rs\",crust)[[2]]/grid_data(\"",trc,"\",\"Bulk_rs\",crust)[[2]]"),colorName="PuBu")
}

#display contour for trace in apatite as a function of total trace
for(trc in traces){
	contour_data("Custom",,paste0(trc," in Ap_rs as a function of total ",trc),96,40,crust_out(),
		"660+2*(grid_x_i-1)","0.1+0.1*(grid_y_i-1)",7,10,graph="contour",write_to_file=TRUE,,
		contour_levels=seq(0.00000000000001,100,5),file_format=".ps",
		custom_def=paste0("grid_data(\"",trc,"\",\"Ap_rs\",crust)[[2]]*grid_data(\"wt%\",\"Ap_rs\",crust)[[2]]/grid_data(\"",trc,"\",\"Bulk_rs\",crust)[[2]]"),colorName="PuBu")
}

#plot vol% for each phase
#axes should be y: Temperature (C); x: K2O (Bulk_rs)
# must normalise to without H2O as a phase (100-grid_data("vol%","H2O_rs",crust_out())[[2]])*100
#Grid_variable_phase is the phase
for(Grid_variable_phase in sort(list_all_phases(crust))){
	contour_data("Custom",Grid_variable_phase,paste0("vol%_",Grid_variable_phase),40,151,crust_out(),
		"(grid_x_i-1)*0.12344+0.1","600+2*(grid_y_i-1)",,,graph="filled",write_to_file=TRUE,,
		contour_levels=seq(0.00000000000001,40,4),file_format=".ps",
		custom_def=paste0("grid_data(\"vol%\",\"",Grid_variable_phase,"\",crust)[[2]]/(100-grid_data(\"vol%\",\"H20_rs\",crust)[[2]])*100"))
	contour_data("Custom",Grid_variable_phase,paste0("vol%_",Grid_variable_phase),40,151,crust_out(),
		"(grid_x_i-1)*0.12344+0.1","600+2*(grid_y_i-1)",,,graph="filled",write_to_file=TRUE,,
		contour_levels=seq(0.00000000000001,40,4),file_format=".png",
		custom_def=paste0("grid_data(\"vol%\",\"",Grid_variable_phase,"\",crust)[[2]]/(100-grid_data(\"vol%\",\"H20_rs\",crust)[[2]])*100"))
}



#plot A/CNK
mwAl2O3<-101.96128; mwCaO<-56.0774; mwNa2O<-61.97894; mwK2O<-94.19600
contour_data("Custom",,"ACNK",40,151,crust_out(),
	"(grid_x_i-1)*0.12344+0.1","600+2*(grid_y_i-1)",,,graph="filled",write_to_file=TRUE,,
	contour_levels=seq(1,2,0.1),file_format=".png",
	custom_def=paste0("(grid_data(\"AL2O3\",\"Bulk_rs\",crust)[[2]]/mwAl2O3)/((grid_data(\"CAO\",\"Bulk_rs\",crust)[[2]]/mwCaO)+(grid_data(\"NA2O\",\"Bulk_rs\",crust)[[2]]/mwNa2O)+(grid_data(\"K2O\",\"Bulk_rs\",crust)[[2]]/mwK2O))"))
contour_data("Custom",,"ACNK",40,151,crust_out(),
	"(grid_x_i-1)*0.12344+0.1","600+2*(grid_y_i-1)",,,graph="filled",write_to_file=TRUE,,
	contour_levels=seq(1,2,0.1),file_format=".ps",
	custom_def=paste0("(grid_data(\"AL2O3\",\"Bulk_rs\",crust)[[2]]/mwAl2O3)/((grid_data(\"CAO\",\"Bulk_rs\",crust)[[2]]/mwCaO)+(grid_data(\"NA2O\",\"Bulk_rs\",crust)[[2]]/mwNa2O)+(grid_data(\"K2O\",\"Bulk_rs\",crust)[[2]]/mwK2O))"))


#plot SiO2 in Melt
contour_data("Custom",,"SiO2 in Melt",40,151,crust_out(),
	"(grid_x_i-1)*0.12344+0.1","600+2*(grid_y_i-1)",,,graph="filled",write_to_file=TRUE,,
	contour_levels=seq(69,75,0.5),file_format=".ps",
	custom_def=paste0("grid_data(\"SIO2\",\"Melt_rs\",crust)[[2]]"))
	
sum(grid_data(grid_data(SIO2,Melt_rs,crust)[[2]]),

# sum of LREE in melt
	list("sum REE in melt",grid_data("La","Melt_rs",crust_out())[[2]]+grid_data("Ce","Melt_rs",crust_out())[[2]]+grid_data("Pr","Melt_rs",crust_out())[[2]]+grid_data("Nd","Melt_rs",crust_out())[[2]]+grid_data("Sm","Melt_rs",crust_out())[[2]]+grid_data("Gd","Melt_rs",crust_out())[[2]])
	c("La","Ce","Pr","Nd","Sm","Gd")

#correctMnzSat
function (kd, c0, pm, cmins, min.props, melt.arg = list(), dont = character(0)) 
{
	#kd = kd.ppx
	#c0 = c0
	#pm = melt wt% (proportion melt?)
	#min.props = mineral proportions
	#melt.arg
		#TT=temp+273.15,
		# mjrs=calc_phases["Melt",major_elements_without_H2O],
		# trc=bpm$cL,
		# H2O=calc_phases["Melt","H2O"])
	#dont?
    c0 <- .sanitize(c0)
	#Sanitize objects that should be one (numeric) line typically mineral proportions or C0. Make sure it's a vector. Remove NA's. Normalize to specified value.
    m.pr <- .sanitize(min.props, normalize.TO = 1)
	#m.pr = sanitized min.props
    lree.nm <- c("Th", "La", "Ce", "Pr", 
        "Nd", "Sm", "Gd")
    mw <- c(232.04, 138.91, 140.12, 140.91, 144.24, 150.3, 157.2, 
        30.974, 15.999)
    names(mw) <- c(lree.nm, "P", "O")
    if (any(toupper(melt.arg$mjrs) == melt.arg$mjrs)) {
        melt.arg$mjrs <- .TrueNames(melt.arg$mjrs)
    } #not sure what this does. .TrueNames is unaccessable
    W.R <- c(-520, -260, -116, 31, 177, 460, 750)
	#W.R = ? some kind of weighting? why are there negative values? W.R usually stands for whole rock
    names(W.R) <- lree.nm
    MW.REE <- mw[lree.nm]
	#MW.REE = molecular weights of LREE
    REE <- melt.arg$trc[lree.nm]
	#REE = trc element conc of liquid (melt)
    M <- rep(0, length(lree.nm))
	#M vector of zeroes, length of LREEs
    Xmz <- 0.99
    milcats <- millications(melt.arg$mjrs)
	#millications = [no. cations in oxide formula] * [Concentration (wt%)] / [M.W.] * 1000
    S <- REE[lree.nm]/MW.REE[lree.nm]
	#S = moles REE of melt. concentration / molecular wt = moles
    ST <- sum(S)
	#ST = total moles in melt
    ee <- mzSaturationWithTh(cats = milcats, REE = REE, H2O = melt.arg$H2O, 
        Xmz = Xmz, Temp = melt.arg$TT - 273.15)
	#ee stores Dmz, Tmz.sat.C, Sum.A.REE.sat
		#D value calculated by (Na + K + 2Ca)/Al * 1/(Al + Si)
		#Temperature of Mz saturation.
		#Sum of REE saturation in melt.
    AiT <- ee[1, "Sum.A.REE.sat"]
    if (ST <= AiT) {
		#if total moles traces in melt from bpm function is <= than calculated from mzSaturationWithTh
        cm <- cmins
        cL <- melt.arg$trc
    }
    else {
        gamma <- exp(W.R/(melt.arg$TT))
		#what is gamma? W.R are those -ve & +ve values for LREE. Divided by temperature?
        .compute <- function(MiT) {
            A <- gamma * AiT * S/(MiT + gamma * AiT)
            ee <- (sum(A) - AiT)^2
            return(ee)
        }
        m.guess <- 0
        out <- optim(m.guess, .compute, method = "Brent", 
            lower = 0, upper = 1)
        R2 <- out$value
        MiT <- out$par
        A <- gamma * AiT * S/(MiT + gamma * AiT)
		#AiT is sum of LREE Sat calculated in mzSaturationWithTh (using Montel 93)
		#S is moles of LREE in melt from BatchPM function.
		#A must be a function of the REE
        M <- S - A
		#S is moles LREE in melt (should contain data for each LREE). Doesn't include P & O.
		#A could be the saturation for each LREE or could be the LREE allocated to Mnz.
			#I lean towards saturation of LREE in melt. Then the excess is S - A, which is LREE for Mnz.
        mf.mnz <- sum((MW.REE + mw["P"] + 4 * mw["O"]) * 
            M)/1e+06
		#mf.mnz = mass fraction mnz? units are wt% i think
		#mf.mnz = total molecular weight
        mnz.prop <- mf.mnz * pm/100
		#mnz.prop = mass fraction mnz * wt% melt/100
        lree.c <- A * MW.REE
        x <- M/sum(M)
        MW.mnz <- sum(x * MW.REE) + mw["P"] + 4 * mw["O"]
        cmnz <- x * MW.REE/MW.mnz * 1e+06
        kd["Mnz", lree.nm] <- cmnz/lree.c
        m.pr["Mnz"] <- mnz.prop
        m.pr <- m.pr/sum(m.pr)
        ee <- BatchPM(kd = kd, c0 = c0, pm = pm, min.props = m.pr)
        cL <- ee$cL
        cm <- ee$cmins
        cL[dont] <- melt.arg$trc[dont]
        all.but.mnz <- setdiff(rownames(cm), "Mnz")
        cm[all.but.mnz, dont] <- cmins[, dont]
        cL[lree.nm] <- lree.c
        cm["Mnz", lree.nm] <- cmnz
    }
    cS <- .sanitize(m.pr %*% cm)
    DD <- cS/cL
    dont <- c(dont, lree.nm)
    invisible(list(cL = .sanitize(cL), cmins = cm, kd = kd, DD = DD, 
        min.props = m.pr, dont = dont))
}
#mzSaturationWithTh
function (cats = milli, REE = filterOut(WR, c("Th", "La", 
    "Ce", "Pr", "Nd", "Sm", "Gd"), 
    1), H2O = 3, Xmz = 0, Temp = 750) 
{
    on.exit(options(show.error.messages = TRUE))
    if (Temp == 0) {
        Temp <- winDialogString("Temperature (degrees C)", 
            "750")
        if (is.null(Temp)) {
            cat("Cancelled.\n")
            options(show.error.messages = FALSE)
            stop()
        }
    }
    Temp <- as.numeric(Temp) + 273.15
    if (is.na(H2O)) 
        H2O <- 1e-07
    if (all(H2O == 0)) {
        H2O <- winDialogString("Water contents in the melt (%)", 
            "3")
        if (is.null(H2O)) {
            cat("Cancelled.\n")
            options(show.error.messages = FALSE)
            stop()
        }
        H2O <- as.numeric(H2O)
    }
    if (Xmz == 0) {
        Xmz <- winDialogString("X mz", "0.99")
        if (is.null(Xmz)) {
            cat("Cancelled.\n")
            options(show.error.messages = FALSE)
            stop()
        }
    }
    Xmz <- as.numeric(Xmz)
    if (nrow(cats) > 1) 
        cats <- cats[rownames(REE), , drop = F]
    lree.nm <- c("Th", "La", "Ce", "Pr", 
        "Nd", "Sm", "Gd")
    mw <- c(232.04, 138.91, 140.12, 140.91, 144.24, 150.3, 157.2, 
        30.974, 15.999)
    names(mw) <- c(lree.nm, "P", "O")
    MW.REE <- mw[lree.nm]
    ree <- t(t(REE)/MW.REE)
	#ree is LREE mols in melt
    ree <- apply(ree, 1, sum, na.rm = TRUE)
	#ree is the sum of LREE in melt
    reex <- ree/Xmz
	#reex could be a stoichiometry of total 'ree' / 'Xmz' where Xmz is set to 0.99. 
	#ree is the total ree mols in melt. so reex is not a stoichiometry of a single phase.
    sums <- apply(cats, 1, sum, na.rm = TRUE)
	#sums is total majors in melt? in millications
    ee <- sapply(1:nrow(cats), function(i) {
        z <- cats[i, ]/sums[i]
        return(z)
    })
	#ee seems to be the proportion of majors in melt.
    x <- t(ee)
	#calculation of Mnz saturation of Montel (1993):
    D.M <- (x[, "Na2O"] + x[, "K2O"] + 2 * x[, "CaO"])/x[, 
        "Al2O3"] * 1/(x[, "Al2O3"] + x[, "SiO2"]) * 
        1
    T.calc <- 13318/(9.5 + 2.34 * D.M + 0.3879 * sqrt(H2O) - 
        log(reex)) - 273.15
    sum.ree.sat <- exp(9.5 + 2.34 * D.M + 0.3879 * sqrt(H2O) - 
        13318/Temp)
    y <- cbind(D.M, round(T.calc, 3), sum.ree.sat)
    colnames(y) <- c("Dmz", "Tmz.sat.C", "Sum.A.REE.sat")
    assign("results", y, .GlobalEnv)
    invisible(y)
}

mzSaturation<-function(cats=milli,REE=filterOut(WR,c("La","Ce","Pr","Nd","Sm","Gd"),1),H2O=3,Xmz=1,Temp=750){ 
    on.exit(options("show.error.messages"=TRUE))
    REE[REE==0]<-NA
    REE<-filterOut(REE,c("La","Ce","Pr","Nd","Sm","Gd"),1)
    browser()
    if(H2O==0){
        H2O<-winDialogString("Water contents in the melt (%)","3")
        if(is.null(H2O)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
        H2O<-as.numeric(H2O)
    }
    
    if(Xmz==0){
        Xmz<-winDialogString("Molar proportion of the LREE monazite (Xmz = 0-1)","0.83")
        if(is.null(Xmz)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}    
    }
    Xmz<-as.numeric(Xmz)
    cats<-cats[rownames(REE),]
    
    oxides<-c("SiO2","TiO2","Al2O3","FeOt","MnO","MgO","CaO","Na2O","K2O","P2O5")
    oxides<-oxides[oxides%in%colnames(cats)]
    cats<-cats[,oxides,drop=FALSE]
    
    sums<-apply(cats,1,sum,na.rm=TRUE)
    ee<-sapply(1:nrow(cats),function(i){
        z<-cats[i,]/sums[i]
        return(z)
    })
    cats<-t(ee)
    
    # Montel 1993
    MW.REE<-c(138.9055,140.12,140.9077,144.24,151.4,154.25)
    names(MW.REE)<-c("La","Ce","Pr","Nd","Sm","Gd")
    
    ree<-t(t(REE)/MW.REE)
    ree<-apply(ree,1,sum,na.rm=TRUE)
    reex<-ree/Xmz
    
	Temp <- as.numeric(Temp) + 273.15
	
    D<-(cats[,"Na2O"]+cats[,"K2O"]+2*cats[,"CaO"])/cats[,"Al2O3"]*1/(cats[,"Al2O3"]+cats[,"SiO2"]) 
    
    T.calc<-13318/(9.5+2.34*D+0.3879*sqrt(H2O)-log(reex))-273.15
    sum.ree.sat <- exp(9.5 + 2.34 * D + 0.3879 * sqrt(H2O) - 
        13318/Temp)
    # Kelsey et al. 2008 GIVES FUNNY RESULTS FOR ANYTHING BUT LEUCOGRANITES
    #FM<-(cats[,"Na2O"]+cats[,"K2O"]+2*(cats[,"CaO"])+cats[,"FeOt"]+cats[,"MgO"])/(cats[,"Al2O3"]*cats[,"SiO2"])
    #ree<-apply(REE,1,sum,na.rm=TRUE)
    
    #D.LREE.K<-566794/ree
    #TMnz.satK.C<--310/(log(D.LREE.K)+1.324*FM-7.5852)-273.15
    
    #y<-cbind(D,round(T.calc,3),FM,round(TMnz.satK.C,3))
    #colnames(y)<-c("Dmz","Tmz.sat.C","FM","TMnz.sat.C (Kelsey)")
    
    y<-cbind(D,round(T.calc,3),sum.ree.sat)
    colnames(y)<-c("Dmz","Tmz.sat.C","sum.ree.sat")
    
    if(!getOption("gcd.shut.up"))print(y)
    
    y<-formatResults(y)
    assign("results",y,.GlobalEnv)
    invisible(y)
}

mz<-correctMnzSat(kd=kd.ppx,
	c0=c0[1:length(c0)-1],
	pm=calc_phases["Melt","wt%"],
	min.props=min.props,
	melt.arg=list(TT=temp+273.15,
		mjrs=calc_phases["Melt",major_elements_without_H2O],
		trc=bpm$cL,
		H2O=calc_phases["Melt","H2O"]),
	cmins=bpm$cmins,dont=character(0))

	crust[[20]][[20]]["Melt_rs",]
mz<-correctMnzSat(kd=kd.ppx,
	c0=c0[1:length(c0)-1],
	pm=calc_phases["Melt","wt%"],
	min.props=min.props,
	melt.arg=list(TT=temp+273.15,
		mjrs=calc_phases["Melt",major_elements_without_H2O],
		trc=bpm$cL,
		H2O=calc_phases["Melt","H2O"]),
	cmins=bpm$cmins,dont=character(0))
###################
#Check accuracy of melt P saturation
P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt",major_elements],
			temp+273.15)
			
getAcc <- function(crust2 = crust,y = y_i, x = x_i){
	Acc_grid <- matrix(0, nrow = y, ncol = x)
	for( i in 1:y){
		for( j in 1:x){
			if(any(rownames(crust2[[i]][[j]]) == "Melt_rs")){
				# browser()
				calc_phases <- crust2[[i]][[j]]
				#store the numbers as {j;i},
				temp <- 660+(j-1)*2
				P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt_rs",major_elements],
					temp+273.15)
				P_Melt <- P_Sat*calc_phases["Melt_rs","wt%"]/100
				Acc_grid[i,j] <- as.numeric(abs((P_Melt-calc_phases["Melt_rs","P2O5"]*calc_phases["Melt_rs","wt%"]/100)/P_Melt*100))
				# browser()
			}
		}
	}
	return(Acc_grid)
}
getT <- function(crust2 = crust,y = y_i, x = x_i){
	Acc_grid <- matrix(0, nrow = y, ncol = x)
	for( i in 1:y){
		for( j in 1:x){
			if(any(rownames(crust2[[i]][[j]]) == "Melt_rs")){
				# browser()
				# calc_phases <- crust2[[i]][[j]]
				#store the numbers as {j;i},
				# P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt_rs",major_elements],
					# temp+273.15)
				# P_Melt <- P_Sat*calc_phases["Melt_rs","wt%"]/100
				# Acc_grid[i,j] <- as.numeric(abs((P_Melt-calc_phases["Melt_rs","P2O5"]*calc_phases["Melt_rs","wt%"]/100)/P_Melt*100))
				temp <- 660+(j-1)*2
				Acc_grid[i,j] <- temp
				# browser()
			}
		}
	}
	return(Acc_grid)
}
meltT <- getT(crust,y_i,x_i)
write(t(flip_y(meltT)), file = paste0("C:/Rcrust/Projects/",working_file,"/Outputs/meltT.txt"),
		ncolumns = x_i,
		append = FALSE)
getMeltP <- function(crust2 = crust,y = y_i, x = x_i){
	Acc_grid <- matrix(0, nrow = y, ncol = x)
	for( i in 1:y){
		for( j in 1:x){
			if(any(rownames(crust2[[i]][[j]]) == "Melt_rs")){
				# browser()
				calc_phases <- crust2[[i]][[j]]
				#store the numbers as {j;i},
				# P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt_rs",major_elements],
					# temp+273.15)
				# P_Melt <- P_Sat*calc_phases["Melt_rs","wt%"]/100
				# Acc_grid[i,j] <- as.numeric(abs((P_Melt-calc_phases["Melt_rs","P2O5"]*calc_phases["Melt_rs","wt%"]/100)/P_Melt*100))
				Acc_grid[i,j] <- as.numeric(calc_phases["Melt_rs","P2O5"])
				# browser()
			}
		}
	}
	return(Acc_grid)
}
meltP <- getMeltP(crust,y_i,x_i)
write(t(flip_y(meltP)), file = paste0("C:/Rcrust/Projects/",working_file,"/Outputs/meltP.txt"),
		ncolumns = x_i,
		append = FALSE)
getMeltSat <- function(crust2 = crust,y = y_i, x = x_i){
	Acc_grid <- matrix(0, nrow = y, ncol = x)
	for( i in 1:y){
		for( j in 1:x){
			if(any(rownames(crust2[[i]][[j]]) == "Melt_rs")){
				# browser()
				calc_phases <- crust2[[i]][[j]]
				#store the numbers as {j;i},
				temp <- 660+(j-1)*2
				P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt_rs",major_elements],
					temp+273.15)
				# P_Melt <- P_Sat*calc_phases["Melt_rs","wt%"]/100
				# Acc_grid[i,j] <- as.numeric(abs((P_Melt-calc_phases["Melt_rs","P2O5"]*calc_phases["Melt_rs","wt%"]/100)/P_Melt*100))
				Acc_grid[i,j] <- as.numeric(P_Sat)
				# browser()
			}
		}
	}
	return(Acc_grid)
}
meltSat <- getMeltSat(crust,y_i,x_i)
write(t(flip_y(meltSat)), file = paste0("C:/Rcrust/Projects/",working_file,"/Outputs/meltSat.txt"),
		ncolumns = x_i,
		append = FALSE)
getAcc <- function(crust2 = crust,y = y_i, x = x_i){
	Acc_grid <- matrix(0, nrow = y, ncol = x)
	for( i in 1:y){
		for( j in 1:x){
			if(any(rownames(crust2[[i]][[j]]) == "Melt_rs")){
				# browser()
				calc_phases <- crust2[[i]][[j]]
				#store the numbers as {j;i},
				temp <- 660+(j-1)*2
				P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt_rs",major_elements],
					temp+273.15)
				Acc_grid[i,j] <- as.numeric(100-abs((P_Sat-calc_phases["Melt_rs","P2O5"])/P_Sat)*100)
				# browser()
			}
		}
	}
	return(Acc_grid)
}
Acc <- getAcc(crust,y_i,x_i)
write(t(flip_y(Acc)), file = paste0("C:/Rcrust/Projects/",working_file,"/Outputs/Acc.txt"),
		ncolumns = x_i,
		append = FALSE)
write("TEST", file = paste0("C:/Rcrust/Projects/",working_file,"/Outputs/TestOutput.txt"),
		ncolumns = x_i,
		append = TRUE)
getMelt <- function(crust2 = crust,y = y_i, x = x_i){
	Acc_grid <- matrix(0, nrow = y, ncol = x)
	for( i in 1:y){
		for( j in 1:x){
			if(any(rownames(crust2[[i]][[j]]) == "Melt_rs")){
				# browser()
				# calc_phases <- crust2[[i]][[j]]
				#store the numbers as {j;i},
				# P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt_rs",major_elements],
					# temp+273.15)
				Acc_grid[i,j] <- as.numeric(crust2[[i]][[j]]["Melt_rs","mass"])
				# browser()
			}
		}
	}
	return(Acc_grid)
}
theMelt <- getMelt(crust,y_i,x_i)
write(t(flip_y(theMelt)), file = paste0("C:/Rcrust/Projects/",working_file,"/Outputs/theMelt.txt"),
		ncolumns = x_i,
		append = FALSE)
#Masses and component proportions check:
	#Multiplies each component weight by the phase mass and sums them. (including P2O5). This gives the total mass of each phase.
	#discrepancies would arise from incorrect or missing component weights
	br <- which(rownames(calc_phases)=="Bulk_rs")
	temp_sum <- c()
	for(i in 1:(br-1)){
		temp_sum <- c(temp_sum,sum(calc_phases[i,c(major_elements,"P2O5")] * calc_phases[i,"mass"] / 100))
	}
	names(temp_sum) <- rownames(calc_phases[1:br-1,])
	print(temp_sum)
	print(sum(temp_sum,calc_phases["Ap","mass"]*0.05,calc_phases["Mnz","mass"]*0.71))
	print(sum(temp_sum,calc_phases["Ap_rs","mass"]*0.05,calc_phases["Mnz_rs","mass"]*0.71))
	#For each component, the component weight of a phase is multiplied by the mass of the phase and summed. This returns the total mass of component.
	br <- which(rownames(calc_phases)=="Bulk_rs")
	temp_sum2 <- c()
	for(i in 6:(5+length(c(major_elements,"P2O5")))){
		temp_sum2 <- c(temp_sum2,sum(calc_phases[1:br-1,i] * calc_phases[1:br-1,"mass"] / 100))
	}
	names(temp_sum2) <- names(calc_phases["Bulk_rs",c(major_elements,"P2O5")])
	# sum(temp_sum2,Ap[,"mass"]*0.05,Mnz[,"mass"]*0.71,melt[,"mass"]*(1-melt[,"P2O5"]/100))
	print(temp_sum2)
	print(sum(temp_sum2,Ap[,"mass"]*0.05,Mnz[,"mass"]*0.71))
	print(sum(temp_sum2,calc_phases["Ap_rs","mass"]*0.05,calc_phases["Mnz_rs","mass"]*0.71))
#masses balance out if the empty mass of Ap, Mnz, melt is accounted for.

(calc_phases["Ap_rs","mass"]*calc_phases["Ap_rs","P2O5"]/100)+
sum((calc_phases["Ap_rs","mass"]*calc_phases["Ap_rs","P2O5"]/100),
	(calc_phases["Ap_rs","mass"]*calc_phases["Ap_rs","P2O5"]/100),
	(calc_phases["Ap_rs","mass"]*calc_phases["Ap_rs","P2O5"]/100))

#Find first points above solidus
getFirstPoints <- function(crust = crust, y_i = y_i, x_i = x_i){
	firstPoints <- ""
	for( i in 1:y_i){
		for( j in 1:x_i){
			if(any(rownames(crust[[i]][[j]]) == "Melt_rs")){
				#store the numbers as {j;i},
				firstPoints <- paste0(firstPoints,paste0("{",j,";",i,"}"),",")
				break
			}
		}
	}
	return(firstPoints)
}
firstPoints <- getFirstPoints(crust,y_i,x_i)
write(firstPoints, file = paste0("C:/Rcrust/Projects/",working_file,"/Outputs/firstPoints.txt"),
		ncolumns = 1,
		append = FALSE)

#function to get points with specific PT for pasting in Rcrust
# temps <- c(700,750,800,850)
# temps <- c(seq(680,750,5))
temps <- seq(680,850,20)
pressures <- c(1.5,3,4.5,6)
# pressures <- c(3)
#solidusPts will include the first points above solidus for each pressure supplied
#if temperatures not supplied, will only return each point above solidus for each pressure
#if pressures not supplied, will not return anything
getPoints <- function(theCrust=crust,y=y_n,x=x_n,T=NULL,P=NULL,reqMelt=TRUE,solidusPts=TRUE){
# browser()
	if(is.null(P)){
		P <- strsplit(winDialogString("Pressures (kbar)",""),",")[[1]]
		if(is.null(P)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
	}
	if(is.null(T)){
		T <- strsplit(winDialogString("Temperatures (C)",""),",")[[1]]
		if(is.null(P)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
	}
	thePoints <- ""
	if(solidusPts){
		for(y_i in 1:y){
			if(!is.na(any(match(input_pt[[y_i]][[x_i]][,"Pressure"],P)))){
				for(x_i in 1:x){
					if(any(rownames(theCrust[[y_i]][[x_i]]) == "Melt_rs")){
						thePoints <- paste0(thePoints,paste0("{",x_i,";",y_i,"}"),",")
						break
					}
				}
			}
		}
	}
	for(y_i in 1:y){
		for(x_i in 1:x){
			# theCrust[[y_i]][[x_i]]
			# browser()
			if(!is.na(any(match(input_pt[[y_i]][[x_i]][,"Pressure"],P))) && !is.na(any(match(input_pt[[y_i]][[x_i]][,"Temperature"],T)))){
				if(reqMelt){
					if(any(rownames(theCrust[[y_i]][[x_i]]) == "Melt_rs")){
						thePoints <- paste0(thePoints,paste0("{",x_i,";",y_i,"}"),",")
					}
				} else {
					thePoints <- paste0(thePoints,paste0("{",x_i,";",y_i,"}"),",")
				}
			}
		}
	}
	return(thePoints)
}
getPoints(T=temps, P= pressures)

#Adapted from GCDkit, Janousek 2006. http://www.gcdkit.org/
millicats <- function (x = c0, print = TRUE, save = FALSE, precision = 3) {
	# browser()
    x <- as.matrix(x)
    if (ncol(x) == 1) {
        x <- t(x)
    }
    oxides <- c("SiO2", "TiO2", "Al2O3", "Fe2O3", 
        "FeO", "FeOt", "MnO", "MgO", 
        "CaO", "Na2O", "K2O", "H2O.PLUS", 
        "CO2", "P2O5", "F", "S")
    # oxides <- oxides[oxides %in% colnames(x)]
	oxides <- toupper(oxides[toupper(oxides) %in% toupper(colnames(x))])
    ee <- sapply(1:nrow(x), function(f) {
		# browser()
		#MW and x.atoms is dependent on GCDModel plugin
		names(MW) <- toupper(names(MW))
		names(x.atoms) <- toupper(names(x.atoms))
        z <- x[f, oxides]/MW[oxides] * x.atoms[oxides] * 1000
        return(z)
    })
    milli <- t(ee)
    milli[is.na(milli)] <- 0
    rownames(milli) <- rownames(x)
    .atoms.from.formula <- function(ox) {
        z <- gsub("[0-9]", "", ox)
        z <- sapply((strsplit(z, "O")), paste, collapse = "")
        return(z)
    }
    if (save) {
        oxides <- oxides[c(1:11, 14)]
    }
    results <- milli[, oxides, drop = FALSE]
    ee <- paste(.atoms.from.formula(oxides), ".m", sep = "")
    i <- grep("^FE.m$", ee)
    if (length(i) == 2) 
        ee[i] <- c("FE3.m", "FE2.m")
    ee <- gsub("^FET.m$", "FE.m", ee)
    colnames(results) <- ee
    # if (!getOption("gcd.shut.up") & print) 
        # print(t(round(results, precision)))
    if (save) {
        assign("milli", milli, .GlobalEnv)
        assign("results", results, .GlobalEnv)
        invisible()
    }
    invisible(milli)
}
#Adapted from GCDkit, Janousek 2006. http://www.gcdkit.org/
filterOut <- function (x, what, n = length(what)) 
{
	x <- as.matrix(x)
	x <- as.matrix(x)
    if (ncol(x) == 1) {
        x <- t(x)
    }
	# browser()
    ii <- as.vector(pmatch(what, colnames(x)))
    names(ii) <- what
    x <- subset(x, TRUE, ii, drop = FALSE)
    if (length(ii) == 1) {
        ee <- apply(x, 1, is.na)
        yy <- x
        yy <- subset(yy, !ee)
    }
    else {
        ee <- apply(x, 1, is.na)
        ee <- apply(as.matrix(ee), 2, sum)
        yy <- subset(x, ee < n, drop = FALSE)
    }
    colnames(yy) <- what
    return(yy)
}
#Adapted from saturation.r of GCDkit, Janousek 2006. http://www.gcdkit.org/
mzSaturation<-function(cats=milli,REE=filterOut(x=c0,c("La","Ce","Pr","Nd","Sm","Gd"),1),H2O=3,Xmz=1,Temp=750){ 
# browser()
    on.exit(options("show.error.messages"=TRUE))
    REE[REE==0]<-NA
    REE<-filterOut(REE,c("La","Ce","Pr","Nd","Sm","Gd"),1)
    # browser()
    if(H2O==0){
        H2O<-winDialogString("Water contents in the melt (%)","3")
        if(is.null(H2O)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
        H2O<-as.numeric(H2O)
    }
    
    if(Xmz==0){
        Xmz<-winDialogString("Molar proportion of the LREE monazite (Xmz = 0-1)","0.83")
        if(is.null(Xmz)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}    
    }
    Xmz<-as.numeric(Xmz)
    # cats<-cats[rownames(REE),]
    
    oxides<-c("SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","P2O5")
	oxides <- toupper(oxides)
    oxides<-oxides[oxides%in%colnames(cats)]
    cats<-cats[,oxides,drop=FALSE]
    
    sums<-apply(cats,1,sum,na.rm=TRUE)
    ee<-sapply(1:nrow(cats),function(i){
        z<-cats[i,]/sums[i]
        return(z)
    })
    cats<-t(ee)
    
    # Montel 1993
    MW.REE<-c(138.9055,140.12,140.9077,144.24,151.4,154.25)
    names(MW.REE)<-c("La","Ce","Pr","Nd","Sm","Gd")
    
    ree<-t(t(REE)/MW.REE)
    ree<-apply(ree,1,sum,na.rm=TRUE)
    reex<-ree/Xmz
    
	Temp <- as.numeric(Temp) + 273.15
	
    D<-(cats[,"NA2O"]+cats[,"K2O"]+2*cats[,"CAO"])/cats[,"AL2O3"]*1/(cats[,"AL2O3"]+cats[,"SIO2"]) 
    
    T.calc<-13318/(9.5+2.34*D+0.3879*sqrt(H2O)-log(reex))-273.15
    sum.ree.sat <- exp(9.5 + 2.34 * D + 0.3879 * sqrt(H2O) - 
        13318/Temp)
		
    # Kelsey et al. 2008 GIVES FUNNY RESULTS FOR ANYTHING BUT LEUCOGRANITES
    #FM<-(cats[,"Na2O"]+cats[,"K2O"]+2*(cats[,"CaO"])+cats[,"FeOt"]+cats[,"MgO"])/(cats[,"Al2O3"]*cats[,"SiO2"])
    #ree<-apply(REE,1,sum,na.rm=TRUE)
    
    #D.LREE.K<-566794/ree
    #TMnz.satK.C<--310/(log(D.LREE.K)+1.324*FM-7.5852)-273.15
    
    #y<-cbind(D,round(T.calc,3),FM,round(TMnz.satK.C,3))
    #colnames(y)<-c("Dmz","Tmz.sat.C","FM","TMnz.sat.C (Kelsey)")
    
    y<-cbind(D,round(T.calc,3),sum.ree.sat)
    colnames(y)<-c("Dmz","Tmz.sat.C","sum.ree.sat")
    
    # if(!getOption("gcd.shut.up"))print(y)
    
    y<-formatResults(y)
    assign("results",y,.GlobalEnv)
    invisible(y)
}
#Adapted from saturation.r of GCDkit, Janousek 2006. http://www.gcdkit.org/
# A general routine that solves nonlinear equation for temperature (deg C) numerically 
# (bisection method)
#(Currently not in use.)
solve.T<-function(fun,Si,ACNK,T.HW,P2O5,tmin=0,tmax=NULL){
	T.calc<-NULL
		if(ACNK>1){
			ttold<-0
			tt<-1
			if(is.null(tmax)) tt.max<-T.HW+273 else tt.max<-tmax 
				# H+W temperature is the only feasible maximum estimate
			tt.min<-tmin
			while(abs(ttold-tt)>0){
				ttold<-tt
				tt<-(tt.max-tt.min)/2+tt.min
				expr<-gsub("Si",Si,fun)
				expr<-gsub("ACNK",ACNK,expr)
				expr<-gsub("T",tt,expr)
				pp<-eval(parse(text=as.expression(expr)))
				if(pp>P2O5)tt.max<-tt else tt.min<-tt
			} 
			T.calc<-tt
		}else{
		T.calc<-NA
		}
	return(T.calc-273.15)
}

neaten_crust<-function(crust_neat=NULL,phase_aliases=NULL){
	#rename using aliases if given
	if(!phase_aliases==""){
		split_aliases<-strsplit(strsplit(phase_aliases,c(","))[[1]],"=")
		phase_aliases<-NULL
		phase_aliases_names<-NULL
		for(i in 1:length(split_aliases)){
		phase_aliases<-c(phase_aliases,split_aliases[[i]][1])
		phase_aliases_names<-c(phase_aliases_names,split_aliases[[i]][2])
		}
		names(phase_aliases)<-phase_aliases_names
		#seperate merge commands
		merge_aliases<-grep("&",phase_aliases,value=TRUE)
		if(length(merge_aliases)!=0){
		phase_aliases<-phase_aliases[-grep("&",phase_aliases)]
		}
		merge_aliases<-strsplit(merge_aliases,"&")
		if(length(merge_aliases)==0){merge_aliases<-""}
		if(length(phase_aliases)==0){phase_aliases<-""}
		if(all(phase_aliases!="")|all(merge_aliases!="")){
			for(x_i in 1:length(crust[[1]])){
				for(y_i in 1:length(crust)){
					delete<-NULL
					grab<-rownames(crust_neat[[y_i]][[x_i]])
					if(length(grab)>0){
					grab_split<-strsplit(grab,"_")
					grab_sys<-grab_name<-grab
					for(ph in 1:length(grab)){
					grab_name[ph]<-grab_split[[ph]][1]
					#if first part of system is numeric extract up to end of _
					if(!is.na(suppressWarnings(as.numeric(grab_split[[ph]][-1][1])))){
					grab_sys[ph]<-paste(grab_split[[ph]][c(-1,-2)],collapse="_")
					}else{
					grab_sys[ph]<-paste(grab_split[[ph]][-1],collapse="_")
					}
					}
					#rename using aliases
					if(all(phase_aliases!="")){
						for(ph in 1:length(grab)){
							if(any(grab_name[ph]==phase_aliases)){
							#Rename feldspars
								if(names(which(grab_name[ph]==phase_aliases))=="Pl|Kf"){
									if(length(intersect(toupper(major_elements),c("CAO","K2O")))==2){
										# fix-tag: currently rename feldspars in two places, in order ot allow phase extraction, simplify this
										# fix-tag: issue with capitalisation
										CaO_pos<-which(toupper(names(crust_neat[[y_i]][[x_i]][ph,]))=="CAO")
										K2O_pos<-which(toupper(names(crust_neat[[y_i]][[x_i]][ph,]))=="K2O")
										if(crust_neat[[y_i]][[x_i]][ph,CaO_pos]>crust_neat[[y_i]][[x_i]][ph,K2O_pos]){
										grab_name[ph]<-"Pl"
										}else{grab_name[ph]<-"Kf"}
									}
								}else{
									#Rename phases using aliases
									grab_name[ph]<-names(which(grab_name[ph]==phase_aliases)[1])
								}
								#label unwanted phases (phases of zero mass or phases labelled as "hide")
								if(crust_neat[[y_i]][[x_i]][ph,"mass"]==0|grab_name[ph]=="hide"){
									delete<-c(delete,ph)
								}
							}
						}
						#number duplicates within systems
						renamed<-paste(grab_name,grab_sys,sep="_")
						num<-1
						while(any(duplicated(renamed))){
							renamed[which(duplicated(renamed))]<-paste(grab_name[which(duplicated(renamed))],num,grab_sys[which(duplicated(renamed))],sep="_")
							num<-num+1
						}
						rownames(crust_neat[[y_i]][[x_i]])<-renamed
					}
					#merge if required
					if(all(merge_aliases!="")){
						systems<-c("rs","es","cumul")
						for(sys in systems){
							for(merge_try in 1:length(merge_aliases)){
								merge_nos<-NULL
								for(i in 1:length(merge_aliases[[merge_try]])){
									if(any(paste(grab_name,grab_sys,sep="_")==paste(merge_aliases[[merge_try]][i],sys,sep="_"))){
										merge_nos<-union(merge_nos,which(grab_name==merge_aliases[[merge_try]][i]))
									}
								}
								merger<-NULL
								if(!is.null(merge_nos)){
									for(i in merge_nos){
										merger<-rbind(merger,crust_neat[[y_i]][[x_i]][i,,drop=FALSE])
									}
									merged<-.wtd.add(merger,avname=names(merge_aliases)[merge_try])
									crust_neat[[y_i]][[x_i]][merge_nos[1],]<-merged
									rownames(crust_neat[[y_i]][[x_i]])[merge_nos[1]]<-paste(names(merge_aliases)[merge_try],sys,sep="_")
									if(length(merge_nos)>1){
										delete<-c(delete,merge_nos[-1])
									}
								}
							}
						}
					}
						#remove unwanted phases
						if(!is.null(delete)){
							crust_neat[[y_i]][[x_i]]<-crust_neat[[y_i]][[x_i]][-delete,]
						}
					}
				}
			}
		}
	}
	return(crust_neat)
}

Unlabelled fields: 2 = Ap, Bt, Grt, H2O, Ilm, Kfs, Melt, Mnz, Qz, Sil; 3 = Ap, Bt, Grt, H2O, Ilm, Kfs, Melt, Mca, Mnz, Qz; 4 = Ap, Bt, Crd, H2O, Ilm, Kfs, Melt, Mnz, Qz, Sil; 7 = Ap, Bt, H2O, Ilm, Kfs, Melt, Mca, Mnz, Qz; 10 = Ap, Bt, Crd, Ilm, Kfs, Melt, Mnz, Qz, Sil; 17 = Ap, Bt, Gt, Ilm, Melt, Mnz, Qz Sil; 21 = Ap, Bt, Crd, Ilm, Kfs, Melt, Mnz, Qz; 22 = Ap, Bt, Crd, H2O, Ilm, Melt, Mnz, Qz; 27 = Bt, H2O, Ilm, Kfs, Mca, Qz, Ru; 28 = Bt, Grt, H2O, Ilm, Kfs, Qz Sil; 29 = Bt, Grt, H2O, Ilm, Kfs, Mca, Qz; 30 = Bt, Crd, H2O, Ilm, Kfs, Qz, Sil; 31 = Bt, Crd, H2O, Ilm, Kfs, Opx, Qz; 32 = Ap, Crd, H2O, Ilm, Kfs, Melt, Opx; 34 = Ap, Bt, Grt, Ilm, Melt, Opx, Qz; 39 = Ap, Bt, Crd, H2O, Ilm, Melt, Mnz; 45 = Bt, H2O, Ilm, Kfs, Qz, Sil; 46 = Bt, H2O, Ilm, Kfs, Mca, Qz; 48 = Ap, Grt, Ilm, Melt, Opx, Qz; 52 = Ap, Bt, Grt, Ilm, Melt, Opx; 54 = Ap, Bt, Crd, Ilm, Melt, Mnz; 55 = Ap, Bt, Crd, H2O, Ilm, Melt; 59 = Ap, Bt, Crd, Ilm, Melt; 








