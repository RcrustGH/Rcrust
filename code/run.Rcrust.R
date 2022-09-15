#############################################
#
# Main Rcrust Loop
#
#############################################
#function-def:run.Rcrust<-function(comps,c0,press,temp,ph_extr_pnt,cumul_extract_pnt=NULL,ph_add_pnt,cumul_add_pnt=NULL)
run.Rcrust<-function(comps,c0,press,temp,ph_extr_pnt,cumul_extract_pnt=NULL,ph_add_pnt,cumul_add_pnt=NULL,...){
  #############################################
  #
  # Ancillary function to merge several lines (by mass)
  #
  ############################################
  .wtd.add<-function(thelines,prop="mass",avname="Averaged"){
    if(nrow(thelines)==1){
      foo<-thelines
    }else{
      # Two sets of cols
      if(exists_and_true(calc_mol)){
        extensive.cn<-c("wt%","vol%","mol%","mol")    # Extensive properties (mass dependant) -- add the others if required
      }else{
        extensive.cn<-c("wt%")
      }

      intensive.cn<-setdiff(colnames(thelines),c(prop,extensive.cn))
      foo<-matrix(rep(0,length(colnames(thelines))),nrow=1)
      colnames(foo)<-colnames(thelines)
      # Intensive
      foo[,intensive.cn]<-thelines[,prop]%*%thelines[,intensive.cn,drop=F]/sum(thelines[,prop])
      # Extensive
      foo[,prop]<-sum(thelines[,prop])
      if(!length(intersect(extensive.cn,thelines))==0){
        foo[,extensive.cn]<-colSums(thelines[,extensive.cn])}
    }
    rownames(foo)<-avname
    return(foo)
  }
  #exists_and_true
exists_and_true<<-function(x){
chk<-try(x,silent=TRUE)
if(class(chk)=="try-error"){return(FALSE)}else{
return(x)
}
}
  #############################################
  #
  # Expression evaluator
  #
  ############################################
  #function-def: eval_expr<-function(expr,calc_phases=calc_phases,crust=crust)
  #Evaluate expression a given calc_phases and crust
		  #() are for solution models, {} are for function terms and bodmas
		  #evaluate any { and word before it up until }
		  eval_expr<-function(expr,calc_phases=calc_phases,crust=crust){
		  #wrap outside for evaluation
		  a<-paste0("{",expr,"}")
		  while(!unlist(gregexpr("[{]",a))[1]==-1){
		  letters<-gregexpr("[a-z]",a)
		  left_bracs<-unlist(gregexpr("[{]",a))
		  names(left_bracs)<-rep("left",length(left_bracs))
		  right_bracs<-unlist(gregexpr("[}]",a))
		  names(right_bracs)<-rep("right",length(right_bracs))
		  bracs<-sort(c(left_bracs,right_bracs))
		  #Find inner brackets
		  i<-1
		  while(!(names(bracs)[i]=="left"&names(bracs)[i+1]=="right")){
		  i<-i+1
		  }
		  #evaluate inner brackets
		  #If require function call find name
		  j<-1
		  while((bracs[i]-j)%in%letters[[1]]){
		  j<-j+1
		  }
		  #if j==1 dont require function call, else j+1 is first letter of function name
		  if(!j==1){
		  funct_name<-substr(a,bracs[i]-j+1,bracs[i]-1)
		  #apply function
		  if(funct_name=="ph"){
		  #arguments of the form ph{phase;unit;x_i;y_i} where unit can be any column name in calc_phases and x_i and y_i are the current point by default
		  ph_args<-strsplit(substr(a,bracs[i]+1,bracs[i+1]-1),"[;,]")[[1]]
		  if(length(ph_args)==2){
		  #current_variable
		  chk_var<-try(eval(parse(text="calc_phases[ph_args[1],ph_args[2]]")),silent=TRUE)
            if(class(chk_var)=="try-error"){
              out<-0
            }else{
              out<-chk_var
            }
			}else{
			#past_variable
			#default look in rs
			if(length(grep("_rs", ph_args[1]))==0&length(grep("_es", ph_args[1]))==0){
			ph_args[1]<-paste0(ph_args[1],"_rs")
			}
		  chk_var<-try(eval(parse(text="crust[[eval(parse(text=ph_args[[4]]))]][[eval(parse(text=ph_args[[3]]))]][ph_args[1],ph_args[2]]")),silent=TRUE)
            if(class(chk_var)=="try-error"){
              out<-0
            }else{
              out<-chk_var
            }
			}
		  }
		  if(funct_name=="delta"){
		  skip_delta<-FALSE
		  #arguments of the form delta{ph;x_#;y_#;unit} where unit can be wt% or mass
		  delta_args<-strsplit(substr(a,bracs[i]+1,bracs[i+1]-1),"[;,]")[[1]]
		  delta_phs<-delta_args[1]
		  if(substring(delta_args[2],1,8)=="prev_ext"){
		  find_ph<-substring(delta_args[2],10)
		  pnt<-"find"
		  l<-1
		  while(pnt=="find"){
		  #find earliest of phase extraction or phase absent point, if neither show warning and take pnt=1
		  #find previous increase in _es_cumul
		  chk_pnt<-try(crust[[y_i]][[x_i-l]][paste0(find_ph,"_es"),"mass"],silent=TRUE)
		  if(class(chk_pnt)=="try-error"){
		  chk_pnt<-try(crust[[y_i]][[x_i-l]][paste0(find_ph,"_rs"),"mass"],silent=TRUE)
		  if(class(chk_pnt)=="try-error"){
		  pnt<-x_i-l
		  cat(paste0("Delta calculated to previous phase absent point at x_i = ",pnt,"\n"))
		  }
		  }else{
		  pnt<-x_i-l
		  cat(paste0("Delta calculated to previous extract at x_i = ",pnt,"\n"))
		  }
		  if((x_i-l)==1&pnt=="find"){
		  pnt<-x_i-l
		  cat(paste0("Delta calculated to first point at x_i = ",pnt,"\n"))
		  }
		  l<-l+1
		  }
		  delta_index_x<-pnt
		  }else{
		  delta_index_x<-eval(parse(text=delta_args[2]))}
		  if(!(delta_index_x<=x_n&delta_index_x>=1)){cat("Warning delta_index_x not <x_n and >=1, skipping this delta calculation\n");skip_delta<-TRUE}
		  delta_index_y<-eval(parse(text=delta_args[3]))
		  if(!(delta_index_y<=y_n&delta_index_y>=1)){cat("Warning delta_index_y not <y_n and >=1, skipping this delta calculation\n");skip_delta<-TRUE}
		  delta_unit<-delta_args[4]
		  if(!(delta_unit=="mass"|delta_unit=="wt%")){
		  cat("Error delta only currently accepted as mass or wt%\n")
		  stop()}
		  if(!skip_delta){
		    current_mode<-0
			previous_mode<-0
			#Check for plus sign (e.g. aluminosilicate given as sill+ky+and)
			delta_phases<-strsplit(delta_args[1],"+",fixed=TRUE)[[1]]
			for(each_ph in delta_phases){
			#current_mode
			chk_var<-try(eval(parse(text="calc_phases[each_ph,delta_unit]")),silent=TRUE)
			if(!class(chk_var)=="try-error"){
			current_mode<-current_mode+chk_var
			}
			#previous_mode
			delta_phase_rs<-paste(each_ph,"_rs",sep="")
			chk_var<-try(eval(parse(text="crust[[delta_index_y]][[delta_index_x]][delta_phase_rs,delta_unit]")),silent=TRUE)
            if(!class(chk_var)=="try-error"){
            previous_mode<-previous_mode+chk_var
            }
			}
		  #delta
		  if(current_mode>previous_mode){
		  out<-current_mode-previous_mode}else{
		  out<-0
		  }
		  }else{out<-0}
		  }
		  if(funct_name=="retain"){
         #c# If Retention Mode - Calculate mass of retention phases to extract
            #c# Extract till retention amount of retention phases is left
			#arguments of the form retain{amount;unit;ph} where ph can be omitted to take on the current ph
		  retain_args<-strsplit(substr(a,bracs[i]+1,bracs[i+1]-1),"[;,]")[[1]]
		  ret<-as.numeric(retain_args[1])
          retention_unit<-retain_args[2]
		  if(!is.na(retain_args[3])){ph<-retain_args[3]}
			bulk_no<-which(rownames(calc_phases)=="Bulk_rs")
			ph_no<-which(rownames(calc_phases)==ph)
            if(retention_unit=="mass"){
			system_less_ret<-calc_phases[-c(bulk_no,ph_no),"mass"]
              if(calc_phases[ph,"mass"]<=ret){
              out<-0
              }else{
              out<-calc_phases[ph,"mass"]-ret
              }
            }
            if(retention_unit=="vol%"){
			system_less_ret<-calc_phases[-c(bulk_no,ph_no),"vol%"]
			if(calc_phases[ph,"vol%"]<=ret){
            out<-0
            }else{
            new_ret_vol<-(ret*(sum(system_less_ret)))/(100-ret)
			new_ret_mass<-new_ret_vol*calc_phases[ph,"mass"]/calc_phases[ph,"vol%"]
            out<-calc_phases[ph,"mass"]-new_ret_mass
              }
            }
			            if(retention_unit=="wt%"){
			system_less_ret<-calc_phases[-c(bulk_no,ph_no),"wt%"]
              if(calc_phases[ph,"wt%"]<=ret){
              out<-0
              }else{
              new_ret_wt<-(ret*(sum(system_less_ret)))/(100-ret)
			  new_ret_mass<-new_ret_wt*calc_phases[ph,"mass"]/calc_phases[ph,"wt%"]
              out<-calc_phases[ph,"mass"]-new_ret_mass
              }
            }
		  }
		  if(funct_name=="return"){
		    #arguments of the form return{phase;amount} where amount can be % or mass
		  return_args<-strsplit(substr(a,bracs[i]+1,bracs[i+1]-1),"[;,]")[[1]]
		  #check if phase is in extract cumul subsystem
		  chk_var<-try(eval(parse(text=paste0("cumul_extract_pnt[\"",return_args[1],"_es_cumul\",\"mass\"]"))),silent=TRUE)
		  if(class(chk_var)!="try-error"){
		  if(!is.null(chk_var)){
		  percentage<-FALSE
		  		  #if percentage tag to calculate
		  if(length(grep("%",return_args[2]))!=0){
		  percentage<-TRUE
		  return_args[2]<-gsub("%","",return_args[2])
		  }
		  #evaluate for number
		    chk_num<-try(eval(parse(text=return_args[2])),silent=TRUE)
            if(class(chk_num)=="try-error"){
		  cat("Error phase addition could not evaluate isolated function correctly")
		  stop()
            }
		  #calculate
		  if(percentage){
		  take<-chk_num/100*chk_var
		  leave<-(100-chk_num)/100*chk_var
		  }else{
		  take<-chk_num
		  leave<-as.numeric(chk_var)-chk_num
		  }
		  if(leave<0){
		  take<-chk_var
		  leave<-0}
		  #transfer
		  cumul_extract_pnt[paste0(return_args[1],"_es_cumul"),"mass"]<<-leave
		  #leave calc
		  a<-paste(c(cumul_extract_pnt[paste0(return_args[1],"_es_cumul"),comps],take),collapse=",")
		  break
		  }else{
		   cat(paste0("Phase ",return_args[1]," not found in extract cumul\n"))
		   a<-paste(rep(0,length(comps)+1),collapse=",")
		   break
		   out<-""
		  }
		  }else{
		   cat(paste0("Phase ",return_args[1]," not found in extract cumul\n"))
		   a<-paste(rep(0,length(comps)+1),collapse=",")
		   break
		   out<-""
		  }
		  }
		  #replace inner brackets with function output
		  a<-paste(substr(a,1,bracs[i]-j),out,substr(a,bracs[i+1]+1,nchar(a)),sep="")
		  }else{
		  #evaluate inner brackets for bodmas output
		  out<-eval(parse(text=substr(a,bracs[i]+1,bracs[i+1]-1)))
		  a<-paste(substr(a,1,bracs[i]-j),out,substr(a,bracs[i+1]+1,nchar(a)),sep="")
		  }
		  }
		  return(a)
		  }
calc<-TRUE
pass<-1
extract<-NULL
addition_subsystem<-NULL
  while(calc){
    # Phase addition
    if(ph_add){
      if(pass==1){
      if(!is.null(ph_add_pnt)){
        # mod-tag:DFM addition (Disequilibrium Fractional Melting calculation)
        if(ph_add_pnt[[1]][[1]]=="dfm"){
          #Disequilibrium Fractional Melting calculation
          current_pnt<-run.Rcrust(comps=comps,c0=c0,
                                  press=input_pt[[y_i]][[x_i]][1],
                                  temp=input_pt[[y_i]][[x_i]][2],
                                  ph_extr_pnt=NULL,
                                  cumul_extract_pnt=NULL,
                                  ph_add_pnt=NULL
          )
          next_pnt<-run.Rcrust(comps=comps,c0=c0,
                               press=input_pt[[y_i]][[x_i+1]][1],
                               temp=input_pt[[y_i]][[x_i+1]][2],
                               ph_extr_pnt=NULL,
                               cumul_extract_pnt=NULL,
                               ph_add_pnt=NULL
          )
          #Pl dfm calculation...  X(Pl consum) = (Na(melt_2) X(melt_2) + Na(Kfs_2) X(Kfs_2) - Na(Ms_1) X(Ms_1))/Na(Pl_2)
          Na_melt_2<-try(next_pnt[paste0(melt.name,"_rs"),"NA2O"],silent=TRUE)
          if(class(Na_melt_2)=="try-error"){Na_melt_2<-0}
          X_melt_2<-try(next_pnt[paste0(melt.name,"_rs"),"wt%"],silent=TRUE)
          if(class(X_melt_2)=="try-error"){X_melt_2<-0}
          Na_Kf_2<-try(next_pnt["Kf_rs","NA2O"],silent=TRUE)
          if(class(Na_Kf_2)=="try-error"){Na_Kf_2<-0}
          X_Kf_2<-try(next_pnt["Kf_rs","wt%"],silent=TRUE)
          if(class(X_Kf_2)=="try-error"){X_Kf_2<-0}
          Na_Ms_1<-try(current_pnt["Mica(CHA)_rs","NA2O"],silent=TRUE)
          if(class(Na_Ms_1)=="try-error"){Na_Ms_1<-0}
          X_Ms_1<-try(current_pnt["Mica(CHA)_rs","wt%"],silent=TRUE)
          if(class(X_Ms_1)=="try-error"){X_Ms_1<-0}
          Na_Pl_2<-try(next_pnt["Pl_rs","NA2O"],silent=TRUE)
          if(class(Na_Pl_2)=="try-error"){cat("Error, Plagioclase not found in second point")}else{
            X_Pl_consm<-(Na_melt_2*X_melt_2+Na_Kf_2*X_Kf_2-Na_Ms_1*X_Ms_1)/Na_Pl_2
            X_Pl_1<-try(current_pnt["Pl_rs","mass"],silent=TRUE)
            if(class(X_Pl_1)=="try-error"){cat("Error, Plagioclase not found in first point")}else{
              if(X_Pl_consm<X_Pl_1){
                #Place unreacted Pl in addition subsystem
                addition_subsystem<-current_pnt["Pl_rs",,drop=FALSE]
                addition_subsystem["Pl_rs","mass"]<-addition_subsystem["Pl_rs","mass"]-X_Pl_consm
                rownames(addition_subsystem)<-"Pl_as"
                #recalculate new c0 with only the reactive Pl amount
                current_pnt["Pl_rs","mass"]<-X_Pl_consm
                current_pnt["Bulk_rs","mass"]<-0
                c0<-.wtd.add(current_pnt)[,names(c0)]
              }
            }
          }
        }
		else{
        #Standard addition
	grab<-c(all_elements,"mass")
    grab_phases<-NULL
  for(add_def_no in 1:length(ph_add_pnt)){
   #Evaluate if addition is required
	statement<-ph_add_pnt[[add_def_no]]["condition"]
  #if statement is true
  if(is.na(as.logical(eval_expr(statement,calc_phases,crust)))){cat(paste0("Error condition for phase addition definition number ",add_def_no," as ",statement," does not evaluate to logical output\n"));stop()}
  if(as.logical(eval_expr(statement,calc_phases,crust))){
    	for(add_ph in 2:length(ph_add_pnt[[add_def_no]])){
		#Evaluate expression unless numeric comma seperated vector
		if(any(is.na(suppressWarnings(as.numeric(strsplit(ph_add_pnt[[add_def_no]][add_ph],",")[[1]]))))){
		ph_add_pnt[[add_def_no]][add_ph]<-eval_expr(ph_add_pnt[[add_def_no]][add_ph],calc_phases,crust)
		}
		}
		  add_mat<-matrix(c(as.numeric(unlist(strsplit(ph_add_pnt[[add_def_no]][-1],","))),c0),length(ph_add_pnt[[add_def_no]]),byrow=TRUE)
          colnames(add_mat)<-c(all_elements,"mass")
          # Calculate new c0
          c0[c(all_elements,"mass")]<-c(.wtd.add(add_mat,prop="mass",avname="Bulk_rs"))
  }else{
    if(!silent_calc){cat("\nNo addition at this point for this condition\n")}
  }
  }
        }
          }
      }
      }
  #Sean-tag
    #Component packet. Removes and stores new component from global variables before calculations occur.
	#Mod-tag: should probably not be altering global variables
	#For new components, always isolate 100% of the mass as the remainder will be not be allocated.
    if(component_packet){
      #variable: new_components
	  #store_c0 holds the original sum of major components, mass of system, and the multiplying factor for proportional calculations
	  store_c0<-c(sum(c0[major_elements]),c0["mass"])
	  names(store_c0)[1]<-"cp_mass"
	  store_c0<-c(store_c0,as.numeric(c0["mass"]/store_c0["cp_mass"]))
	  names(store_c0)[3]<-"mass_factor"
      store_cp<-NULL
      if(!isTRUE(new_components=="")){
        for(i in 1:length(new_components)){
          #store the new component data before removing from major_elements and c0, so as to not interfere with calc_phases
          store_cp<-c(store_cp,c0[match(new_components[i],names(c0))])
          c0<-c0[-match(new_components[i],names(c0))]
		  if(!is.na(match(new_components[i],major_elements))){
			major_elements<<-major_elements[-match(new_components[i],major_elements)]
			#Fix-tag:test this with traces
			all_elements<<-all_elements[-match(new_components[i],all_elements)]
		  }
		  c0["mass"]<-c0["mass"]-store_cp[i]/store_c0["cp_mass"]*c0["mass"]
        }
      }
    }
	#Calculate phases
 comps<<-major_elements
calc_phases<-try(wrapper(comps,c(c0[1:(length(major_elements))],c0[length(c0)]),press,temp,calc_choice),silent=TRUE)
if(class(calc_phases)=="try-error"){
cat("Oops, an error occured while calculating phases\n at ",press," kbar and ",temp," C for the bulk composition:\n ",all_elements,"\n",c0,"\n")
run_errors<<-c(run_errors,y_i,x_i)
if(pause_on_error){
browser()
}else{
exit_calc<<-TRUE
}
}else{
exit_calc<<-FALSE
}
calc<-FALSE
#Sean-tag
#Component packet. Built by Sean Hoffman 2020.
if(component_packet){
	phs_missing<-c("")
	#check if all phases to be partitioned to exist.
	if(apply_trace_correction!="None"&apply_trace_correction!=""){
		if(!any(rownames(calc_phases)=="Melt")){
			phs_missing<-c("Melt")
		}
	}else{
		for(i in 1:length(cp_components)){
			cp_phases<-eval(parse(text=paste0('cp_phases_',names(cp_components[i]))))
			for(j in 1:length(cp_phases)){
				cp_phase<-names(cp_phases[j])
				if(is.na(match(cp_phase,rownames(calc_phases)))){
					phs_missing<-paste0(phs_missing," ",cp_phase)
				}
			}
		}
	}#End of phase checking
	#Mod-tag: likely make new phase decision here
	if(phs_missing==""){
	 cp_calc<-0
	 cp_component<-0
	for(i in 1:length(cp_components)){
		cp_isolate<-cp_components[[i]]
		cp_component[i]<-names(cp_components[i])
		#working with pre-existing component
		if(is.na(match(cp_component[i],new_components))){
			if(length(grep("%",cp_isolate))!=0){#calculating as percentage of component
				cp_isolate<-gsub("%","",cp_isolate)
				cp_calc[i]<-c0[cp_component[i]]*as.numeric(eval_expr(cp_isolate,calc_phases))/100
			}else{#calculate as numeric value
				cp_calc[i]<-as.numeric(eval_expr(cp_isolate,calc_phases))
			}
			#Error validation and removal from component
			if(cp_calc[i]<=0){
				stop()
			}else if(cp_calc[i]<=c0[cp_component[i]]){
				c0[cp_component[i]]<-c0[cp_component[i]]-cp_calc[i]
			}else{
				cp_calc[i]<-c0[cp_component[i]]
				c0[cp_component[i]]<-0
			}
			#adjust mass proportionately to ensure correct partioning
			#Mod-tag: should we be removing from mass at all?
			cp_calc[i]<-as.numeric(cp_calc[i]*store_c0["mass_factor"])
			c0["mass"]<-c0["mass"]-cp_calc[i]
		}else{#new components don't use c0 in calculations, component is in store_cp
			if(length(grep("%",cp_isolate))!=0){#calculating as percentage of component
			cp_isolate<-gsub("%","",cp_isolate)
			cp_calc[i]<-store_cp*as.numeric(eval_expr(cp_isolate,calc_phases))/100
			}else{#calculate as numeric value
				cp_calc[i]<-as.numeric(eval_expr(cp_isolate,calc_phases))
			}
			#Error validation and removal from component
			if(cp_calc[i]<0){
				stop()
			}else if(cp_calc[i]<=store_cp){
				store_cp<-store_cp-cp_calc[i]
			}else{
				cp_calc[i]<-store_cp
				store_cp<-0
			}
			cp_calc[i]<-as.numeric(cp_calc[i]*store_c0["mass_factor"])
		}
	}
	#recalculate calc_phases through meemum before next section, which continues after traces
	calc_phases<-try(wrapper(comps,c(c0[1:(length(major_elements))],c0[length(c0)]),press,temp,calc_choice),silent=TRUE)
	}else{
		cat("The phases,",phs_missing,", are not present at this point.\n")
	}
}#Component packet continues after traces.
#renaming felspars
	if(!exit_calc){
if(!is.na(match("Bulk_rs",rownames(calc_phases)))){
  #rename kf
  if(length(intersect(toupper(all_elements),c("CAO","K2O")))==2){
  split_names<-strsplit(rownames(calc_phases),"_","")
  first_names<-NULL
  for(i in 1:length(split_names)){
  first_names<-c(first_names,split_names[[i]][1])
  }
  for(ph in which(first_names=="Fsp")){
  CaO_pos<-which(toupper(names(calc_phases[ph,]))=="CAO")
  K2O_pos<-which(toupper(names(calc_phases[ph,]))=="K2O")
    if(calc_phases[ph,K2O_pos]<=0){calc_phases[ph,K2O_pos]<-0.0001}
    if(calc_phases[ph,CaO_pos]/calc_phases[ph,K2O_pos]>1){
      rownames(calc_phases)[ph]<-"Pl"
    }else{
    rownames(calc_phases)[ph]<-"Kf"
    }
  }
  }
  #number duplicates
    for(ph in rownames(calc_phases)[which(duplicated(rownames(calc_phases)))]){
        rownames(calc_phases)[which(rownames(calc_phases)==ph)[-1]]<-paste0(ph,"_",1:length(which(rownames(calc_phases)==ph)[-1]))
    }
}
}
#Sean-tag
#Calculate Trace element partitioning and apply Zircon, Monazite and Apatite saturation corrections
# browser()
if(calculate_traces){
	#Mod-tag: Making subsolidus routine within traces for apatite, although traces don't work subsolidus,
		#because we assume that all P2O5 will reside in apatite in subsolidus conditions .
	if(any(rownames(calc_phases)=="Melt")){
		#Apatite saturation modifies calc_phases and includes P2O5, therefore min.props, kd.ppx, bpm are calculated later.
		if(apply_trace_correction!="Apatite saturation"){
			min.props<-calc_phases[c(-which(rownames(calc_phases)=="Melt"),-which(rownames(calc_phases)=="Bulk_rs")),"wt%"]
			kd.ppx<-ppxCleanKdTable(kd.in,ppxPhases=names(min.props),interactive=FALSE)
			## BatchPM function from GCDmodel written by Jean Francois Moyen 2019
			# kd = file containing partition coefficients for trace elements in minerals, Mineral names in kd file must match abbreviation names from Perple_X
			# c0 = numeric vector containing major and trace element composition of reactive subsystem
			# pm = proportion of melt in wt.% of reactive subsystem
			# min.props = mineral proportions in wt.% of reactive subsystem
			# cmins = mineral compositions (*to be passed through function if necessary)
			# melt.args = additional parameters for example temperature dependence of kds
			# dont = phases to not be considered
			library(stats)
			bpm<-BatchPM(kd=kd.ppx,
				c0=c0[1:length(c0)-1],
				pm=calc_phases["Melt","wt%"],
				min.props=min.props,
				cmins=matrix(),melt.arg=list(),dont=character(0))
		}
		if(any(major_elements=="H2O")){major_elements_without_H2O<-major_elements[-which(major_elements=="H2O")]}
		else{major_elements_without_H2O<-major_elements}
			if(apply_trace_correction!="None"&apply_trace_correction!=""){
				# TT is in Kelvin
				if(apply_trace_correction=="Zircon Saturation (Watson & Harrison 1983)"){
					mz<-correctZrnSat(kd=kd.ppx,
						c0=c0[1:length(c0)-1],
						pm=calc_phases["Melt","wt%"],
						min.props=min.props,
						melt.arg=list(TT=temp+273.15,
									mjrs=calc_phases["Melt",major_elements_without_H2O],
									trc=bpm$cL,
									H2O=calc_phases["Melt","H2O"]
									),
						cmins=bpm$cmins,dont=character(0))
				}
				if(apply_trace_correction=="Monazite saturation (Montel 1993)"){
				# browser()
					mz<-correctMnzSat(kd=kd.ppx,
						c0=c0[1:length(c0)-1],
						pm=calc_phases["Melt","wt%"],
						min.props=min.props,
						melt.arg=list(TT=temp+273.15,
							mjrs=calc_phases["Melt",major_elements_without_H2O],
							trc=bpm$cL,
							H2O=calc_phases["Melt","H2O"]),
						cmins=bpm$cmins,dont=character(0))
				}
				if(apply_trace_correction=="Zircon and Monazite saturation"){
					mz<-correctZrnMnzSat(kd=kd.ppx,
					c0=c0[1:length(c0)-1],
					pm=calc_phases["Melt","wt%"],
					min.props=min.props,
					melt.arg=list(TT=temp+273.15,
						mjrs=calc_phases["Melt",major_elements_without_H2O],
						trc=bpm$cL,
						H2O=calc_phases["Melt","H2O"]
						),
					cmins=bpm$cmins,dont=character(0))
				}
				#Apatite saturation routine written by Sean Hoffman 2021
				#Routine can be found in apSaturation.r
				#Apatite saturation currently only deals with P2O5, CaO, H2O, and traces
				if(apply_trace_correction=="Apatite saturation"){
					#normalise c0 major elements, including P2O5.
					store_c0 <- c0
					if(!is.na(match("P2O5",names(c0)))){
						c0[c(major_elements,"P2O5")] <- c0[c(major_elements,"P2O5")]/sum(c0[c(major_elements,"P2O5")])*c0["mass"]
						calc_phases <- correctApSat(c0=c0[c(major_elements,"P2O5","mass")],
							temp=temp,
							press=press,
							calc_phases=calc_phases,
							apatite_saturation=apatite_saturation,
							major_elements=major_elements)
						#Check that Bulk_rs values match the (modified) c0.
						if(all(round(calc_phases["Bulk_rs", names(c0[c(major_elements,"P2O5")])],10) == round(c0[c(major_elements,"P2O5")],10))){
							calc_phases["Bulk_rs", names(c0[c(major_elements,"P2O5")])]<-c0[c(major_elements,"P2O5")]
						}
						#Add trace elements to calc_phases, specifically modified for ApSat.
						#Mod-tag: Duplicating code is not ideal but it works for now.
						if(!all(trace_elements=="P2O5")){
							min.props<-calc_phases[c(-which(rownames(calc_phases)=="Melt"),
								-which(rownames(calc_phases)=="Bulk_rs")),"wt%"]
							kd.ppx<-ppxCleanKdTable(kd.in,ppxPhases=names(min.props),interactive=FALSE)
							library(stats)
							bpm<-BatchPM(kd=kd.ppx,
								c0=c0[1:length(c0)-1],
								pm=calc_phases["Melt","wt%"],
								min.props=min.props,
								cmins=matrix(),melt.arg=list(),dont=character(0))
							#add trace element data to calc_phases
							trace_mat<-matrix(NA,nrow(calc_phases),length(trace_elements)-1)
							colnames(trace_mat)<-trace_elements[-which(trace_elements=="P2O5")]
							rownames(trace_mat)<-rownames(calc_phases)
							trace_mat["Melt",]<-bpm$cL
							trace_mat["Bulk_rs",]<-c0[trace_elements[-which(trace_elements=="P2O5")]]
							for(ph in rownames(bpm$cmins)){
								trace_mat[ph,]<-bpm$cmins[ph,]
							}
							#add in trace elements after majors
							last_major<-which(colnames(calc_phases)=="P2O5")
							calc_phases<-cbind(cbind(calc_phases[,1:last_major],trace_mat),
											calc_phases[,(last_major+1):ncol(calc_phases)])
						}
					} else {
						cat("\nP2O5 is not selected under trace elements\n")
						stop()
					}
				}else if(apply_trace_correction=="Apatite & Monazite Saturation"){
					store_c0 <- c0
					if(!is.na(match("P2O5",names(c0)))){
						c0[c(major_elements,"P2O5")] <- c0[c(major_elements,"P2O5")]/sum(c0[c(major_elements,"P2O5")])*c0["mass"]
						# mz<-correctMnzSat(kd=kd.ppx,
							# c0=c0[1:length(c0)-1],
							# pm=calc_phases["Melt","wt%"],
							# min.props=min.props,
							# melt.arg=list(TT=temp+273.15,
								# mjrs=calc_phases["Melt",major_elements_without_H2O],
								# trc=bpm$cL,
								# H2O=calc_phases["Melt","H2O"]),
							# cmins=bpm$cmins,dont=character(0))
						# calc_phases <- correctApMnzSatWithCa(c0=c0,
							# kd=kd.ppx,
							# temp=temp,
							# press=press,
							# calc_phases=calc_phases,
							# apatite_saturation=apatite_saturation,
							# major_elements=major_elements,
							# mz_cL=mz$cL,
							# min.props=min.props,
							# cmins=bpm$cmins,
							# melt.arg=list(TT=temp+273.15,
								# mjrs=calc_phases["Melt",major_elements_without_H2O],
								# trc=bpm$cL,
								# H2O=calc_phases["Melt","H2O"]),
							# dont=character(0))
						calc_phases <- correctApMnzSatWithCa(c0=c0,
							kd=kd.ppx,
							temp=temp,
							press=press,
							calc_phases=calc_phases,
							apatite_saturation=apatite_saturation,
							major_elements=major_elements,
							Xmz=as.numeric(Xmz))
						c0 <- store_c0
					} else {
						cat("\nP2O5 is not selected under trace elements\n")
						stop()
					}
				}else{
					#Add in monazite and zircon
					trace_mat<-matrix(NA,length(mz$min.props)+2,length(trace_elements))
					colnames(trace_mat)<-trace_elements
					rownames(trace_mat)<-c("Melt",names(mz$min.props),"Bulk_rs")
					trace_mat["Melt",]<-mz$cL
					trace_mat["Bulk_rs",]<-c0[trace_elements]
					for(ph in rownames(mz$cmins)){
						trace_mat[ph,]<-mz$cmins[ph,]
					}
					#Expand trace mat
					last_major<-which(colnames(calc_phases)==major_elements[length(major_elements)])
					prenames<-colnames(calc_phases)[1:last_major]
					postnames<-colnames(calc_phases)[-(1:last_major)]
					premat<-matrix(NA,nrow(trace_mat),length(prenames))
					postmat<-matrix(NA,nrow(trace_mat),length(postnames))
					full_mat<-cbind(cbind(premat,trace_mat),postmat)
					colnames(full_mat)<-c(prenames,colnames(trace_mat),postnames)
					for(row_x in rownames(calc_phases)){
						for(col_x in colnames(calc_phases)){
							full_mat[row_x,col_x]<-calc_phases[row_x,col_x]
						}
					}
					#replace melt composition with corrected composition (not sure if needed here)
					full_mat["Melt",trace_elements]<-mz$cL
					#renormalise solids on wt.%
					normalised_solids<-mz$min.props*(100-full_mat["Melt","wt%"])
					full_mat[names(mz$min.props),"wt%"]<-round(normalised_solids,2)
					calc_phases<-full_mat
					calc_phases[,"mass"]<-calc_phases[,"wt%"]/100*calc_phases["Bulk_rs","mass"]
					#we assume 0 for major elements for zrc and mnz
					if(any(rownames(calc_phases)=="Zrn")){
						calc_phases["Zrn",major_elements]<-0
					}
					if(any(rownames(calc_phases)=="Mnz")){
						calc_phases["Mnz",major_elements]<-0
					}
				}#End: saturation routines other than apatite saturation
			}else{#Mod-tag: to get traces into apatite, apatite saturation needs to come in here or duplicate the code.
				#This section triggers is if trace correction is not selected. But apatite is part of trace corrections..
				#Trace element data without saturation corrections
				#add trace element data to calc_phases
				trace_mat<-matrix(NA,nrow(calc_phases),length(trace_elements))
				colnames(trace_mat)<-trace_elements
				rownames(trace_mat)<-rownames(calc_phases)
				trace_mat["Melt",]<-bpm$cL
				trace_mat["Bulk_rs",]<-c0[trace_elements]
				for(ph in rownames(bpm$cmins)){
					trace_mat[ph,]<-bpm$cmins[ph,]
				}
				#add in trace elements after majors
				last_major<-which(colnames(calc_phases)==major_elements[length(major_elements)])
				calc_phases<-cbind(cbind(calc_phases[,1:last_major],trace_mat),calc_phases[,(last_major+1):ncol(calc_phases)])
			}

	}else{ #subsolidus routines
		#Mod-tag: subsolidus is mutually exclusive with trace element partitioning.
		if(apply_trace_correction=="Apatite saturation"){
			#subsolidus apatite routine. P2O5 assumed to occur in apatite.
			if(!is.na(match("P2O5",names(c0)))){
				store_c0 <- c0
				c0[c(major_elements,"P2O5")] <- c0[c(major_elements,"P2O5")]/sum(c0[c(major_elements,"P2O5")])*c0["mass"]
				calc_phases <- correctApSat(c0=c0[c(major_elements,"P2O5","mass")],
					temp=temp,
					press=press,
					calc_phases=calc_phases,
					apatite_saturation=apatite_saturation,
					major_elements=major_elements)
				#Check that Bulk_rs values match the (modified) c0.
				if(all(round(calc_phases["Bulk_rs", names(c0[c(major_elements,"P2O5")])],10) == round(c0[c(major_elements,"P2O5")],10))){
					calc_phases["Bulk_rs", names(c0[c(major_elements,"P2O5")])]<-c0[c(major_elements,"P2O5")]
				}
				#There is an edge case where isolating CaO from the bulk and then recalculating stable phases
					#will increase hydrate bulk sufficiently to allow Melt to form. Thus traces must be partitioned.
					#P2O5 is still assumed to only occur in apatite.
				if(any(rownames(calc_phases)=="Melt")){
					#Add trace elements to calc_phases, specifically modified for ApSat.
					#Mod-tag: Duplicating code is not ideal but it works for now.
					if(!all(trace_elements=="P2O5")){
						min.props<-calc_phases[c(-which(rownames(calc_phases)=="Melt"),
												-which(rownames(calc_phases)=="Bulk_rs")),"wt%"]
						kd.ppx<-ppxCleanKdTable(kd.in,ppxPhases=names(min.props),interactive=FALSE)
						library(stats)
						bpm<-BatchPM(kd=kd.ppx,
								c0=c0[1:length(c0)-1],
								pm=calc_phases["Melt","wt%"],
								min.props=min.props,
								cmins=matrix(),melt.arg=list(),dont=character(0))
						#add trace element data to calc_phases
						trace_mat<-matrix(NA,nrow(calc_phases),length(trace_elements)-1)
						colnames(trace_mat)<-trace_elements[-which(trace_elements=="P2O5")]
						rownames(trace_mat)<-rownames(calc_phases)
						trace_mat["Melt",]<-bpm$cL
						trace_mat["Bulk_rs",]<-c0[trace_elements[-which(trace_elements=="P2O5")]]
						for(ph in rownames(bpm$cmins)){
							trace_mat[ph,]<-bpm$cmins[ph,]
						}
						#add in trace elements after majors
						last_major<-which(colnames(calc_phases)=="P2O5")
						calc_phases<-cbind(cbind(calc_phases[,1:last_major],trace_mat),
										calc_phases[,(last_major+1):ncol(calc_phases)])
					}
				}
			} else {
				cat("\nP2O5 is not selected under trace elements\n")
				stop()
			}
		}#end of subsolidus apatite saturation
	}#end of subsolidus trace/saturation routines
}#end of traces & saturation.

#Renaming of feldspars to be usable by component packet and phase extraction.
if(!exit_calc){
if(!is.na(match("Bulk_rs",rownames(calc_phases)))){
  #rename kf
  if(length(intersect(toupper(all_elements),c("CAO","K2O")))==2){
  split_names<-strsplit(rownames(calc_phases),"_","")
  first_names<-NULL
  for(i in 1:length(split_names)){
  first_names<-c(first_names,split_names[[i]][1])
  }
  for(ph in which(first_names=="Fsp")){
  CaO_pos<-which(toupper(names(calc_phases[ph,]))=="CAO")
  K2O_pos<-which(toupper(names(calc_phases[ph,]))=="K2O")
    if(calc_phases[ph,K2O_pos]<=0){calc_phases[ph,K2O_pos]<-0.0001}
    if(calc_phases[ph,CaO_pos]/calc_phases[ph,K2O_pos]>1){
      rownames(calc_phases)[ph]<-"Pl"
    }else{
    rownames(calc_phases)[ph]<-"Kf"
    }
  }
  }
  #number duplicates
    for(ph in rownames(calc_phases)[which(duplicated(rownames(calc_phases)))]){
        rownames(calc_phases)[which(rownames(calc_phases)==ph)[-1]]<-paste0(ph,"_",1:length(which(rownames(calc_phases)==ph)[-1]))
    }
}
}

#Component packet partition section
if(component_packet){
	#component packet part 2
	if(phs_missing==""){	
	#adds new row cp_packet to calc_phases
	calc_phases<-rbind(calc_phases,0)
	rownames(calc_phases)[length(rownames(calc_phases))]<-"cp_packet"
	for(i in 1:length(cp_components)){
		#add new component to calc_phases matrix
		if(!is.na(match(cp_component[i],new_components))){
			aa<-matrix(0,nrow(calc_phases),1)
			colnames(aa)<-names(cp_components[i])
			last_major<-which(colnames(calc_phases)==major_elements[length(major_elements)])
			calc_phases<-cbind(cbind(calc_phases[,1:last_major],aa),calc_phases[,(last_major+1):ncol(calc_phases)])
			major_elements<-c(major_elements,names(cp_components[i]))
		}
		cp_phases<-eval(parse(text=paste0('cp_phases_',names(cp_components[i]))))
		calc_phases["cp_packet","mass"]<-cp_calc[i]
		calc_phases["cp_packet",names(cp_components[i])]<-100
		for(j in 1:length(cp_phases)){
			#new row of values to be calculated from cp_packet
			cp_add_to_phase<-calc_phases["cp_packet",]
			if(cp_phases[[j]]=="excess"){
				cp_add_to_phase["mass"]<-calc_phases["cp_packet","mass"]
			}else if(length(grep("%",cp_phases[[j]]))!=0){
				#calculating as percentage of amount in cp_packet, first remove the % symbol
				cp_phases[[j]]<-gsub("%","",cp_phases[[j]])
				#calculate mass to be added to phase from cp_packet
				cp_add_to_phase["mass"]<-cp_calc[i]*(as.numeric(eval_expr(cp_phases[[j]],calc_phases))/100)
			}else{#mass to be added to phase is calculated as numeric value of amount in component packet
				if(grepl("[A-Za-z]",cp_phases[[j]])){
					cp_add_to_phase["mass"]<-as.numeric(eval_expr(cp_phases[[j]],calc_phases))
				} else{
					cp_add_to_phase["mass"]<-as.numeric(eval_expr(cp_phases[[j]],calc_phases))*mass_factor
				}
			}
			#partition remainder of packet if attempting to partition more than available
			if(round(cp_add_to_phase["mass"],10)>round(calc_phases["cp_packet","mass"],10)){
					cat("Warning: attempted to partition portion of ",names(cp_components[i])," to ",names(cp_phases[j])," that is greater than remainder in packet.\n")
					cp_add_to_phase["mass"]<-calc_phases["cp_packet","mass"]
			}
			#adds mass to phase as weighted calculation
			calc_phases[match(names(cp_phases[j]),rownames(calc_phases)),]<-.wtd.add(rbind(cp_add_to_phase,
				calc_phases[match(names(cp_phases[j]),rownames(calc_phases)),]))
			#transfer mass from cp_packet to Bulk_rs.
			#Mod-tag: should we be modifying the bulk mass?
			calc_phases["cp_packet","mass"]<-round(calc_phases["cp_packet","mass"]-cp_add_to_phase["mass"],8)
			calc_phases["Bulk_rs","mass"]<-round(calc_phases["Bulk_rs","mass"]+cp_add_to_phase["mass"],5)
			calc_phases["Bulk_rs",cp_component[i]]<-calc_phases["Bulk_rs",cp_component[i]]+cp_add_to_phase["mass"]/store_c0["mass"]*100
		}
		#normalise wt% of phases and proportions of major components
		mass_phases<-calc_phases[1:which(rownames(calc_phases)=="Bulk_rs")-1,"mass"]
		calc_phases[1:which(rownames(calc_phases)=="Bulk_rs")-1,"wt%"]<-mass_phases/calc_phases["Bulk_rs","mass"]*100
		mass_majors<-calc_phases["Bulk_rs",major_elements]
		calc_phases["Bulk_rs",major_elements] <- mass_majors/sum(mass_majors)*100
		#set the portion of component to 0 if nothing left in packet.
		#Mod-tag: should the remainder in packet always go into some phase?
		if(calc_phases["cp_packet","mass"]<=0){
			calc_phases["cp_packet",names(cp_components[i])]<-0
		}
	}
	#remove cp_packet from calc_phases if no mass left
	if(calc_phases["cp_packet","mass"]<=0){
		calc_phases<-calc_phases[-which(rownames(calc_phases)=="cp_packet"),]
	}
	}
}#end of component packet

if(!exit_calc){
  #Phase Extraction
if(ph_extr){
    #pass_1 extraction
if(!is.null(ph_extr_pnt)){
if(pass==1){
	grab<-c(all_elements,"mass")
    grab_phases<-NULL
  for(ex_def_no in 1:length(ph_extr_pnt)){
   #Evaluate if extraction is required
	statement<-ph_extr_pnt[[ex_def_no]]["condition"]
  #if statement is true
  if(is.na(as.logical(eval_expr(statement,calc_phases,crust)))){cat(paste0("Error condition for phase extraction definition number ",ex_def_no," as ",statement," does not evaluate to logical output\n"));stop()}
  if(as.logical(eval_expr(statement,calc_phases,crust))){
	#must evaluate conditions in order then phases in order, retain phases must be last, cannot have overlapping extract conditions that target the same phase
    extr_phases<-names(ph_extr_pnt[[ex_def_no]])[2:length(ph_extr_pnt[[ex_def_no]])]
	if(any(extr_phases=="any_phase")){
	#old_extr_phases<-extr_phases
	while(length(which(extr_phases=="any_phase"))!=0){
	rep_no<-which(extr_phases=="any_phase")[1]
	#Remove bulk
	cur_phases<-rownames(calc_phases)[-nrow(calc_phases)]
	#Remove phases with existing definitions
	existing_defs<-NULL
	for(i in 1:length(ph_extr_definitions)){
	existing_defs<-union(existing_defs,names(ph_extr_definitions[[i]]))
	}
	cur_phases<-setdiff(cur_phases,existing_defs)
	x<-as.list(extr_phases)
	x[[rep_no]]<-cur_phases
	extr_phases<-unlist(x)
	}
	extr_vals<-unlist(lapply(1:length(extr_phases),function(i){
	gsub("any_phase",extr_phases[i],ph_extr_pnt[[ex_def_no]]["any_phase"])}
	))
	names(extr_vals)<-extr_phases
	ex_def<-names(ph_extr_pnt[[ex_def_no]])[-1]
	ex_def<-ex_def[-which(ex_def=="any_phase")]
	extr_vals[ex_def]<-ph_extr_pnt[[ex_def_no]][ex_def]
	}else{
	extr_vals<-ph_extr_pnt[[ex_def_no]][extr_phases]
	}
    for(ph in extr_phases){
	#If extract phase is present
      if(length(which(rownames(calc_phases)==ph))>0){
	  #grab phase details
        chk<-try(calc_phases[ph,grab,drop=FALSE],silent=TRUE)
        if(!class(chk)=="try-error"){
          grab_phases<-rbind(grab_phases,calc_phases[ph,grab,drop=FALSE])
        }
		#Evaluate for extraction value unless ending in % sign
		a<-extr_vals[ph]
		if(!substring(extr_vals[ph],nchar(extr_vals[ph]),nchar(extr_vals[ph]))=="%"){
		a<-eval_expr(a,calc_phases,crust)
		}
		#Evaluate unless extraction value is non numeric
		if(!is.na(suppressWarnings(as.numeric(gsub("%","",a))))){
		#If number then take smallest of number or mass present
		if(regexpr("%",a)[1]==-1){
		class(a)<-"numeric"
		 if(a<calc_phases[ph,"mass"]){
                  grab_phases[ph,"mass"]<-as.numeric(a)
                  calc_phases[ph,"mass"]<-calc_phases[ph,"mass"]-grab_phases[ph,"mass"]
				}else{
                  calc_phases[ph,"mass"]<-0
                }
		}else{
		 #If percentage then take as proportion
		grab_phases[ph,"mass"]<-round(calc_phases[ph,"mass"]*as.numeric(gsub("%","",a))/100,6)
        calc_phases[ph,"mass"]<-round(calc_phases[ph,"mass"],6)-grab_phases[ph,"mass"]
		}
		}else{
		cat("Error in expression used for extraction\n")
        stop()
		}
     }
    }
  }else{
    if(!silent_calc){cat("\nNo extraction at this point for this condition\n")}
  }
  }
   #Check for negatives
	if(any(calc_phases[,"mass"]<0)|any(grab_phases[,"mass"]<0)){
    cat("Error negative mass produced upon extraction\n")
    }
  #h# Create Extract
  if(!sum(grab_phases[,"mass"])==0){
    #c# Place extracts into extract and calculate Bulk_es
    extr_bulk<-.wtd.add(grab_phases,prop="mass",avname="Bulk_es")
    rownames(grab_phases)<-paste0(rownames(grab_phases),"_es")
    extract<-rbind(grab_phases,extr_bulk)

    #c# Calculate new c0
	c0[c(all_elements,"mass")]<-c(.wtd.add(calc_phases[-which(rownames(calc_phases)=="Bulk_rs"),c(all_elements,"mass")]))
	#Recalculate mass dependent properties in calc_phases if not re-equilibrating
	if(!reequilibrate_steps){
	calc_phases[,"mol"][-which(rownames(calc_phases)=="Bulk_rs")]<-calc_phases[,"mass"][-which(rownames(calc_phases)=="Bulk_rs")]/(calc_phases[,"wt%"][-which(rownames(calc_phases)=="Bulk_rs")]/calc_phases[,"mol"][-which(rownames(calc_phases)=="Bulk_rs")])
	calc_phases[,"mol%"][-which(rownames(calc_phases)=="Bulk_rs")]<-calc_phases[,"mol"][-which(rownames(calc_phases)=="Bulk_rs")]/sum(calc_phases[,"mol"][-which(rownames(calc_phases)=="Bulk_rs")])*100
	calc_phases[,"wt%"][-which(rownames(calc_phases)=="Bulk_rs")]<-calc_phases[,"mass"][-which(rownames(calc_phases)=="Bulk_rs")]/sum(calc_phases[,"mass"][-which(rownames(calc_phases)=="Bulk_rs")])*100
	vols<-calc_phases[,"mass"][-which(rownames(calc_phases)=="Bulk_rs")]/calc_phases[,"Density(kg/m3)"][-which(rownames(calc_phases)=="Bulk_rs")]
	calc_phases[,"vol%"][-which(rownames(calc_phases)=="Bulk_rs")]<-vols/sum(vols)*100
	calc_phases["Bulk_rs",names(c0)]<-c0
	#mod-tag: check that changed all mass dependent properties
	}
	}
  #if(changed c0 comp) then trigger for recalculation
  if(reequilibrate_steps){
  if(!all(calc_phases["Bulk_rs",names(c0)]==c0)){
    calc<-TRUE
    }
	}
}
}
#end extraction
}
  pass<-pass+1
}
 }
#end calculation
#if successful calc phases
if(!exit_calc){
if(!is.na(match("Bulk_rs",rownames(calc_phases)))){
  #rename kf
  if(length(intersect(toupper(all_elements),c("CAO","K2O")))==2){
  #for(ph in which(rownames(calc_phases)=="Fsp")){
	#Mod-tag: changed to search for rownames containing Fsp instead of matching.
	for(ph in grep("Fsp",rownames(calc_phases))){
  CaO_pos<-which(toupper(names(calc_phases[ph,]))=="CAO")
  K2O_pos<-which(toupper(names(calc_phases[ph,]))=="K2O")
    if(calc_phases[ph,K2O_pos]<=0){calc_phases[ph,K2O_pos]<-0.0001}
    if(calc_phases[ph,CaO_pos]/calc_phases[ph,K2O_pos]>1){
      rownames(calc_phases)[ph]<-"Pl"
    }else{
    rownames(calc_phases)[ph]<-"Kf"
    }
  }
  }
  #number duplicates
    for(ph in rownames(calc_phases)[which(duplicated(rownames(calc_phases)))]){
        rownames(calc_phases)[which(rownames(calc_phases)==ph)[-1]]<-paste0(ph,"_",1:length(which(rownames(calc_phases)==ph)[-1]))
    }
  #add rs labels
  rownames(calc_phases)[-match("Bulk_rs",rownames(calc_phases))]<-paste0(rownames(calc_phases)[-match("Bulk_rs",rownames(calc_phases))],"_rs")
}
#Update combined extracts
combined_extracts<-NULL
#Expand extract
if(!is.null(extract)){
expanded<-matrix(0,nrow(extract),ncol(calc_phases))
rownames(expanded)<-rownames(extract)
colnames(expanded)<-colnames(calc_phases)
expanded[rownames(extract),colnames(extract)]<-extract
extract<-expanded
}
#Calculate new cumul_extract_pnt
if(is.null(cumul_extract_pnt)){
  new_cumul_extract<-extract
}else{
rownames(cumul_extract_pnt)<-gsub("_cumul","",rownames(cumul_extract_pnt))
all_extr<-union(rownames(extract),rownames(cumul_extract_pnt))
new_cumul_extract<-matrix(0,length(all_extr),ncol(calc_phases))
rownames(new_cumul_extract)<-all_extr
for(nam_i in all_extr){
  #In extract and cumul - .wtd.add
  if(length(which(rownames(extract)==nam_i))>0&length(which(rownames(cumul_extract_pnt)==nam_i))>0){
    new_cumul_extract[nam_i,]<-.wtd.add(rbind(cumul_extract_pnt[nam_i,],extract[nam_i,]))
  }
  #In extract- extract
  if(length(which(rownames(extract)==nam_i))>0&!length(which(rownames(cumul_extract_pnt)==nam_i))>0){
    new_cumul_extract[nam_i,]<-extract[nam_i,]
  }
  #In cumul- cumul
  if(!length(which(rownames(extract)==nam_i))>0&length(which(rownames(cumul_extract_pnt)==nam_i))>0){
    new_cumul_extract[nam_i,]<-cumul_extract_pnt[nam_i,]
  }
}
#Calculate new Bulk_es_cumul
colnames(new_cumul_extract)<-colnames(cumul_extract_pnt)
new_cumul_extract<-rbind(new_cumul_extract[-(which(rownames(new_cumul_extract)=="Bulk_es")),,drop=FALSE],.wtd.add(new_cumul_extract[-(which(rownames(new_cumul_extract)=="Bulk_es")),,drop=FALSE],,"Bulk_es"))
}
if(!is.null(new_cumul_extract)){
rownames(new_cumul_extract)<-paste0(rownames(new_cumul_extract),"_cumul")
}
combined_extracts<-rbind(extract,new_cumul_extract)
#Update IS
#Calculate new cumul_add_pnt
if(is.null(cumul_add_pnt)){
  new_cumul_add<-addition_subsystem
}else{
  rownames(cumul_add_pnt)<-gsub("_cumul","",rownames(cumul_add_pnt))
  all_add<-union(rownames(addition_subsystem),rownames(cumul_add_pnt))
  new_cumul_add<-matrix(0,length(all_add),ncol(calc_phases))
  rownames(new_cumul_add)<-all_add
  for(nam_i in all_add){
    #In add and cumul - .wtd.add
    if(length(which(rownames(addition_subsystem)==nam_i))>0&length(which(rownames(cumul_add_pnt)==nam_i))>0){
      new_cumul_add[nam_i,]<-.wtd.add(rbind(cumul_add_pnt[nam_i,],addition_subsystem[nam_i,]))
    }
    #In addition_subsystem           - addition_subsystem
    if(length(which(rownames(addition_subsystem)==nam_i))>0&!length(which(rownames(cumul_add_pnt)==nam_i))>0){
      new_cumul_add[nam_i,]<-addition_subsystem[nam_i,]
    }
    #In cumul             - cumul
    if(!length(which(rownames(addition_subsystem)==nam_i))>0&length(which(rownames(cumul_add_pnt)==nam_i))>0){
      new_cumul_add[nam_i,]<-cumul_add_pnt[nam_i,]
    }
  }
}
if(!is.null(new_cumul_add)){
  rownames(new_cumul_add)<-paste0(rownames(new_cumul_add),"_cumul")
}
combined_add<-rbind(addition_subsystem,new_cumul_add)
  #Compile Full System
full_system<-rbind(calc_phases,combined_extracts,combined_add)
}else{
# Return Blank Comp
full_system<-matrix(0,2,length(c0)+1,dimnames=list(c("Error","Bulk_error"),c("wt%",names(c0))))
}
  #Output results to console
if(silent_calc){
#flush.console()
}else{
cat("-------------------------------\n","Input:\n",paste0("x_i = ",x_i,"y_i = ",y_i,"P = ",press," kbar"," ;     T = ",temp," C\n"))
print(c0)
cat("-------------------------------\n","Output:\n")
print(full_system[,c(all_elements,"mass")])
flush.console()
}
return(full_system)
}