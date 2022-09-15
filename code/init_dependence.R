#############################################
#
#  Determine the dependence relations between points
#
############################################
# Bulk Composition *For now dependence is only on composition
#possible handled variables	 dependents = rs,fs,es,as	non-dependents = c0,x_i,y_i,x_n,y_n
#handles functions of the form"rs{x_i-1,y_i}+5-as{x_n-2,y_i*2}"
#fixtag error handling
dependencies<-rep(list(rep(list(),x_n)),y_n) 
#for each point
for(y_i in 1:y_n){
for(x_i in 1:x_n){
dependence<-NULL
# for each major element
for(maj_i in 1:length(all_elements)){
#pull out tuples that major element for that point is dependent on
  #fix-tag: look at main.r for bulk comp to allow multiple rs calls in one tuple
point<-gsub("rs","",gsub("es","",gsub("as","",gsub("fs","",input_bulk[[y_i]][[x_i]][maj_i]))))
left_brac<-gregexpr("\\{", point)[[1]]
right_brac<-gregexpr("\\}", point)[[1]]
tuples<-NULL
for(tup_i in 1:length(left_brac)){
tuples<-c(tuples,substr(point, left_brac[tup_i], right_brac[tup_i]))
}
if(tuples==""){dependence="base"}else{
split_tuples<-strsplit(gsub("\\{","",gsub("\\}","",tuples)),split=";")
for(split_i in 1:length(split_tuples)){
dependence<-c(dependence,paste(eval(parse(text=split_tuples[[split_i]][1])),eval(parse(text=split_tuples[[split_i]][2])),sep=";"))
}
}
}
#Keep common dependence
dependencies[[y_i]][[x_i]]<-union(dependence,dependence)
}
}
########################
# Determine tiers
########################
tiers<-rep(list(rep(list(NULL),x_n)),y_n)
eval_dep<-dependencies
for(y_i in 1:y_n){
  for(x_i in 1:x_n){
tier_counter<-0
while(is.null(tiers[[y_i]][[x_i]])){
#fixtag -see if this alteration works
if(eval_dep[[y_i]][[x_i]][length(eval_dep[[y_i]][[x_i]])]=="base"){
#if(eval_dep[[y_i]][[x_i]][1]=="base"){
#if(eval_dep[[y_i]][[x_i]]=="base"){
  tiers[[y_i]][[x_i]]<-1+tier_counter
}else{
tier_counter<-tier_counter+1
#Evaluate dependencies
new_dep<-NULL
for(dep_i in 1:length(eval_dep[[y_i]][[x_i]])){
  dep_pos<-strsplit(eval_dep[[y_i]][[x_i]][dep_i],split=";")[[1]]
 suppressWarnings(new_dep<-c(new_dep,dependencies[[as.numeric(dep_pos[2])]][[as.numeric(dep_pos[1])]]))
}
#keep common dependence
eval_dep[[y_i]][[x_i]]<-union(new_dep,new_dep)
}
}
}
}
########################
# Determine calculation order
########################
tiers_union<-union(unlist(tiers),unlist(tiers))
calc_order<-rep(list(rep(list(rep(list(),1)),1)),length(tiers_union))
group_num<-0
for(i in tiers_union){
  group_num<-group_num+1
  first_in_group<-TRUE
for(y_i in 1:y_n){
  for(x_i in 1:x_n){
if(tiers[[y_i]][[x_i]]==i){
if(first_in_group){
  calc_num<-0
}
calc_num<-calc_num+1
first_in_group<-FALSE
calc_order[[group_num]][[calc_num]]<-list(x_i=x_i,y_i=y_i)
}
}
}
}
#############################################
#
#  Error validation
#
############################################
for(y_i in 1:y_n){
  for(x_i in 1:x_n){
    if(is.null(dependencies[[y_i]][[x_i]])){cat("Error: No dependence defined for x_i =",x_i," y_i =",y_i," \n");stop()}
  }
}
if(!length(unlist(calc_order))/2==x_n*y_n){cat("Error: Calculation order could not be determined\n");stop()}
cat("Done with dependence determination\n")
cat("..............................................\n")