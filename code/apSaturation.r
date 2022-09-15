#apSaturation.r written by Sean Hoffman, 2021-2022
#Mod-tag: .wtd.add function should be taken out of this file and placed into a useful functions file
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
#Mod-tag: calcACNK function should be taken out of this file and placed into a useful functions file
#Coded by Sean Hoffman, 2021
#General function to return the molar A/CNK ratio of comp.
#Function accepts a vector object that must be named as component oxides (including Al2O3, CaO, Na2O, K2O), wt.% units
calcACNK<-function(comp){
	mwAl2O3 <- 101.96128; mwCaO <- 56.0774; mwNa2O <- 61.97894; mwK2O <- 94.19600
	# Resolving case-sensitive components.
	Al <- comp[[which(toupper(names(comp))=="AL2O3")]]
	Ca <- comp[[which(toupper(names(comp))=="CAO")]]
	Na <- comp[[which(toupper(names(comp))=="NA2O")]]
	K <- comp[[which(toupper(names(comp))=="K2O")]]
	ACNK <- (Al/mwAl2O3)/((Ca/mwCaO)+(Na/mwNa2O)+(K/mwK2O))
	return(ACNK)
}
#Apatite saturation functions. Each returns the saturation concentration of P2O5 (wt.%) in melt.
#Functions correctApSatHW, correctApSatBea, correctApSatPich, correctApSatWolfLondon
	#adapted from saturation.r of GCDkit, Janousek 2006 http://www.gcdkit.org/
#Apatite Saturation model for metaluminous systems by Harrison & Watson, 1984.
correctApSatHW<-function(Cmelt,T){
	#Equation uses Si wt% of melt, Temperature (K)
	if(!is.na(match("SiO2",names(Cmelt)))){Si<-Cmelt["SiO2"]}else{Si<-Cmelt["SIO2"]}
	Si<-Si/100
	# Harrison and Watson (1984)
	A<-8400+(Si[[1]]-0.5)*26400
	B<-3.1+12.4*(Si[[1]]-0.5)
	D.HW<-exp(A/T-B)
	P2O5.HW<-42/D.HW
	return(P2O5.HW)
}
#Correction for peraluminous systems by Bea et al., 1992
correctApSatBea<-function(Cmelt,T){
	#equation uses mol A/CNK of melt, Si wt% of melt, Temperature (C)
	ACNK<-calcACNK(Cmelt)
	# Harrison and Watson (1984)
	P2O5.HW<-correctApSatHW(Cmelt,T)
	#Bea et al. 1992
	P2O5.Bea<-NULL
	E<-(ACNK-1)*6429
	P2O5.Bea<-P2O5.HW*exp(E/(T-273.15))
	return(P2O5.Bea)
}
#Correction for peraluminous systems by Pichavant et al., 1992
correctApSatPich<-function(Cmelt,T){
	#equation uses mol A/CNK of melt, Si wt% of melt, Temperature
	ACNK<-calcACNK(Cmelt)
	# Harrison and Watson (1984)
	P2O5.HW<-correctApSatHW(Cmelt,T)
	if(!is.na(match("SiO2",names(Cmelt)))){Si<-Cmelt["SiO2"]}else{Si<-Cmelt["SIO2"]}
	Si<-Si/100
	# Pichavant et al. (1992)  
	P2O5.PV<-NULL
	C<--5900
	D<--3.22*Si[[1]]+9.31  
	P2O5.PV<-P2O5.HW+(ACNK-1)*exp(C/T+D)   
	
	#T.PV<-solve.T("42/exp((8400+(Si-0.5)*26400)/T-3.1-12.4*(Si-0.5))+(ACNK-1)*exp(-5900/T-3.22*Si+9.31)")
	return(P2O5.PV)
}
#Correction for peraluminous systems by Wolf & London, 1994
correctApSatWolfLondon <- function (Cmelt){
	#WL94 equation uses ASI of melt
	ACNK<-calcACNK(Cmelt)
	P2O5.WL <- -3.4+3.1*ACNK
	return(P2O5.WL)
}
#Monazite saturation of Stepanov, 2012.
#applicable to melts with low CaO, MgO, FeO and ASI >0.85
#temperature is in K. pressure is in kbar.
#Xmz is the mole ratio of LREEmnz / (LREEmnz + Ymnz + Thmnz +Umnz). default is 0.83
correctMnzSatStepanov <- function(press, temp, Cmelt, Xmz=0.83){
	H2O <- as.numeric(Cmelt[["H2O"]])
	temp <- as.numeric(temp)
	press <- as.numeric(press)
	#Stepanov, 2012
	return(exp(16.16 + (0.23 * sqrt(H2O)) - (11494 / temp) - (19.4 * (press / temp)) + log(Xmz)))
}
#correctApSat written by Sean Hoffman, 2021
correctApSat <- function(c0,temp,press,calc_phases,apatite_saturation,major_elements,override=FALSE,Old_P_Melt=0){
#Mod-tag: Should change variable names to those used in description in Thesis/article.
#Naming protocol:
	#[Phase]_[Component] is the wt.% of [Component] in the [Phase]
	#[Component]_[Phase] is the mass of [Component] in the [Phase]. Or wt.% of the component as a proportion of the bulk composition.
	#First part solves for P_Melt after isolating portion of CaO for apatite.
		#P_Melt is the P2O5 wt% of bulk (or P2O5 mass) to saturate melt with respect to apatite.
	continue <- FALSE
	#Apatite component proportions. Average values of silicate magmas from Webster & Piccoli (2015)
	propCaO <- 54; propP2O5 <- 41
	counter <- 1
	if(any(rownames(calc_phases)=="Melt")){
		#Melt_P is melt P2O5 wt% required to saturate with respect to apatite
		if(apatite_saturation=="Harrison & Watson 1984"){
			Melt_P <- correctApSatHW(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="H&W with Bea et al. 1992"){
			Melt_P <- correctApSatBea(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="H&W with Pichavant et al. 1992"){
			Melt_P <- correctApSatPich(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="Wolf & London 1994"){
			Melt_P <- correctApSatWolfLondon(Cmelt=calc_phases["Melt",major_elements])
		}
		#1.2
		#Multiply Melt_P/100 by Melt wt% to get P2O5 wt% of bulk to saturate Melt with respect to apatite
		P_Melt <-calc_phases["Melt","wt%"]*Melt_P/100
	} else { #subsolidus. All P goes to apatite.
		P_Melt <- 0
	}
	#iterate to solve for P_Melt.
	#mid_P_Melt is initialized before the loop but only used after the 2nd loop if the itn_threshold is not yet achieved.
	mid_P_Melt <- 0
	while(!continue){
		#the midpoint between last 2 values (mid_P_Melt) is used to solve for the correct P_Melt (of melt) below a certain threshold.
			#a lower threshold (itn_threshold) takes longer to calculate. <10% may take 2-3 iterations. <1% may take 2-5 iterations, etc.
			#Limit to iterations has been set to 5, as there can be scenarios where an infinite loop occurs.
			#This is especially the case with W&L94 equation.
		store_c0 <- c0
		#P_Bulk is concentration of P2O5 in bulk composition. This value never changes.
		#P_Ap is concentration of P2O5 allocated to apatite. Value can only be <= P_Bulk and >= 0
		P_Bulk <- c0[["P2O5"]]
		#Sometimes a negative saturation value is found, treated as 0.
		if(P_Melt<0 || mid_P_Melt<0){ P_Melt <- 0; mid_P_Melt <- 0; continue <- TRUE}
		#2
		if(counter>2){
			if(mid_P_Melt >= P_Bulk){
				P_Ap <- 0
				continue <- TRUE
			} else {
				P_Ap <- P_Bulk - mid_P_Melt
			}
		} else {
			if(P_Melt >= P_Bulk){
				P_Ap <- 0
				continue <- TRUE
			} else {
				P_Ap <- P_Bulk - P_Melt
			}
		}
		#Mod-tag: should continue case be here? To skip calc_phases if apatite does not saturate at all?
		#calculate required CaO for apatite
		Ca_Ap <- P_Ap / propP2O5 * propCaO
		if(!is.na(match("CaO",names(c0)))){Ca <- "CaO"}else{Ca <- "CAO"}
		if(Ca_Ap <= c0[Ca]){
			#removing portions of components from c0, which are reserved for apatite and P2O5 in melt
			#then recalculating calc_phases
			c0[Ca] <- c0[Ca]-Ca_Ap
			c0["mass"] <- c0["mass"]-sum(Ca_Ap,c0["P2O5"])
			calc_phases <- try(wrapper(comps,c(c0[1:(length(major_elements))],c0[length(c0)]),press,temp,calc_choice),silent=TRUE)
		} else {
			cat("Error: Not enough calcium to form apatite")
			stop()
		}
		#store P_Melt as Old_P_Melt to be compared after isolating components and recalculating stable phases.
		Old_P_Melt <- P_Melt
		#continue if negative saturation value or if apatite does not saturate.
		if(!continue){
			if(any(rownames(calc_phases)=="Melt")){
				#calc P_Melt
				if(apatite_saturation=="Harrison & Watson 1984"){
					Melt_P <- correctApSatHW(Cmelt=calc_phases["Melt",major_elements],
						temp+273.15)
				} else if(apatite_saturation=="H&W with Bea et al. 1992"){
					Melt_P <- correctApSatBea(Cmelt=calc_phases["Melt",major_elements],
						temp+273.15)
				} else if(apatite_saturation=="H&W with Pichavant et al. 1992"){
					Melt_P <- correctApSatPich(Cmelt=calc_phases["Melt",major_elements],
						temp+273.15)
				} else if(apatite_saturation=="Wolf & London 1994"){
					Melt_P <- correctApSatWolfLondon(Cmelt=calc_phases["Melt",major_elements])
				}
				#Multiply Melt_P/100 by Melt wt% to get P2O5 wt% of bulk to saturate Melt
				P_Melt <-calc_phases["Melt","wt%"]*Melt_P/100
			} else { #subsolidus calculation
				#All P2O5 goes to apatite
				P_Melt <- 0
			}
			#exit iteration if saturation is 0 (or melt not stable)
			if(P_Melt == 0 && Old_P_Melt == 0){
				continue <- TRUE
			} else {
				#Change itn_threshold value to set iteration threshold (%)
				itn_threshold <- 1
				#counter > 2 uses midpoint value as calculation.
				if(counter>2){
					if(isTRUE(all.equal(P_Melt,mid_P_Melt))){
						continue <- TRUE
					} else if(abs((P_Melt-mid_P_Melt)/P_Melt*100) < itn_threshold){
						continue <- TRUE
					} else if(counter >4){
						print(paste0("Did not resolve to itn_theshold (",itn_threshold,
							"%) in 5 iterations. ApSat iteration %change ",abs((P_Melt-mid_P_Melt)/P_Melt*100)))
						continue <- TRUE
					}else if(counter > 1){
						#if counter is on 2 or more and has not been solved to satisfactory itn_threshold, iteration needs to be assisted.
						#Take last 2 P_Melt values and find halfway point.
						mid_P_Melt <- (P_Melt + Old_P_Melt)/2
					}
				}else{
					if(isTRUE(all.equal(P_Melt,Old_P_Melt))){
						continue <- TRUE
					} else if(abs((P_Melt-Old_P_Melt)/P_Melt*100) < itn_threshold){
						continue <- TRUE
					} else if(counter > 1){
						#if counter is on 2 or more and has not been solved to satisfactory itn_threshold, iteration needs to be assisted.
						#Take last 2 P_Melt values and find halfway point.
						mid_P_Melt <- (P_Melt + Old_P_Melt)/2
					}
				}
			}
		}
		c0 <- store_c0
		counter <- counter + 1
		# browser()
	}
	if(counter > 3){ final_P_Sat <- mid_P_Melt } else { final_P_Sat <- Old_P_Melt }
	#This section partitions components to melt and/or adds apatite to calc_phases and formats it
	#Ca_Ap, P_Ap = mass of each component allocated to apatite. P_Melt = mass of P2O5 allocated melt. P_Bulk is bulk P2O5 mass
	store_c0 <- c0
	if(any(rownames(calc_phases)=="Melt")){
		last_major <- which(colnames(calc_phases)==major_elements[length(major_elements)])
		if(P_Ap==0){#Apatite doesn't form, allocate all P2O5 to melt
			#Add P2O5 column to calc_phases
			aa <- matrix(0,nrow(calc_phases),1)
			colnames(aa) <- "P2O5"
			calc_phases <- cbind(cbind(calc_phases[,1:last_major],aa),calc_phases[,(last_major+1):ncol(calc_phases)])
			#Add P2O5 to melt
			aa <- matrix(0,nrow = 1, ncol(calc_phases))
			colnames(aa) <- colnames(calc_phases)
			rownames(aa) <- "Melt"
			aa["Melt","mass"] <- P_Bulk
			aa["Melt","P2O5"] <- 100
			calc_phases["Melt",] <- .wtd.add(rbind(aa,calc_phases["Melt",]))
			calc_phases["Bulk_rs",] <- .wtd.add(rbind(aa,calc_phases["Bulk_rs",]))
		}else{#2. apatite forms. Melt holds P2O5 up to P_Melt. P_Ap for apatite.
			#Add P2O5 column to calc_phases
			aa <- matrix(0,nrow(calc_phases),1)
			colnames(aa) <- "P2O5"
			calc_phases <- cbind(cbind(calc_phases[,1:last_major],aa),calc_phases[,(last_major+1):ncol(calc_phases)])
			#Adds P2O5 to Melt up to calculated saturation
			aa <- matrix(0,nrow = 1, ncol(calc_phases))
			colnames(aa) <- colnames(calc_phases)
			rownames(aa) <- "Melt"
			aa["Melt","mass"] <- P_Melt
			aa["Melt","P2O5"] <- 100
			calc_phases["Melt",] <- .wtd.add(rbind(aa,calc_phases["Melt",]))
			calc_phases["Bulk_rs",] <- .wtd.add(rbind(aa,calc_phases["Bulk_rs",]))
			#Adds phase, apatite ("Ap"), to calc_phases.
			aa <- matrix(0,nrow = 1, ncol(calc_phases))
			colnames(aa) <- colnames(calc_phases)
			rownames(aa) <- "Ap"
			aa[,Ca] <- propCaO
			aa[,"P2O5"] <- propP2O5
			aa[,"mass"] <- Ca_Ap / 54 * 100
			br <- which(rownames(calc_phases)=="Bulk_rs")
			calc_phases <- rbind(rbind(calc_phases[1:br-1,],aa),calc_phases[br:nrow(calc_phases),,drop=FALSE])
			calc_phases["Bulk_rs",] <- .wtd.add(rbind(aa,calc_phases["Bulk_rs",]))
			calc_phases["Ap",(last_major+2):length(aa)]<-NaN
			calc_phases["Ap","Density(kg/m3)"] <- 3190
			calc_phases["Ap","mol"] <- aa[,"mass"]/(aa[,Ca]/100*56.0774 + aa[,"P2O5"]/100*141.9445)
		}
	} else {#subsolidus routine
		last_major <- which(colnames(calc_phases)==major_elements[length(major_elements)])
		#Add P2O5 column to calc_phases
		aa <- matrix(0,nrow(calc_phases),1)
		colnames(aa) <- "P2O5"
		calc_phases <- cbind(cbind(calc_phases[,1:last_major],aa),calc_phases[,(last_major+1):ncol(calc_phases)])
		#Adds phase, apatite ("Ap"), to calc_phases.
		aa <- matrix(0,nrow = 1, ncol(calc_phases))
		colnames(aa) <- colnames(calc_phases)
		rownames(aa) <- "Ap"
		aa[,Ca] <- propCaO
		aa[,"P2O5"] <- propP2O5
		aa[,"mass"] <- Ca_Ap / 54 * 100
		br <- which(rownames(calc_phases)=="Bulk_rs")
		calc_phases <- rbind(rbind(calc_phases[1:br-1,],aa),calc_phases[br:nrow(calc_phases),,drop=FALSE])
		calc_phases["Bulk_rs",] <- .wtd.add(rbind(aa,calc_phases["Bulk_rs",]))
		calc_phases["Ap",(last_major+2):length(aa)]<-NaN
		calc_phases["Ap","Density(kg/m3)"] <- 3190
		calc_phases["Ap","mol"] <- aa[,"mass"]/(aa[,Ca]/100*56.0774 + aa[,"P2O5"]/100*141.9445)
	}
	# browser()
	#floating point errors may accumulate, neatening output:
	br <- which(rownames(calc_phases)=="Bulk_rs")-1
	if(sum(calc_phases[1:br,"mass"])!=as.numeric(c0["mass"])){
		calc_phases[1:br,"mass"] <- calc_phases[1:br,"mass"]/sum(calc_phases[1:br,"mass"])*c0["mass"]
	}
	calc_phases[1:br,"mass"] <- calc_phases[1:br,"mass"]/sum(calc_phases[1:br,"mass"])*c0["mass"]
	calc_phases["Bulk_rs","mass"] <- 100
	#wt%
	calc_phases[1:br,"wt%"] <- calc_phases[1:br,"mass"]/sum(calc_phases[1:br,"mass"])*100
	calc_phases["Bulk_rs","wt%"] <- 100
	#vol%
	volume <- calc_phases[1:br,"mass"]/calc_phases[1:br,"Density(kg/m3)"]
	calc_phases[1:br,"vol%"] <- volume/sum(volume)*100
	calc_phases["Bulk_rs","vol%"] <- 100
	#mol%
	calc_phases[1:br,"mol%"] <- calc_phases[1:br,"mol"]/sum(calc_phases[1:br,"mol"])*100
	calc_phases["Bulk_rs","mol%"] <- 100
	#Mod-tag: renaming feldspars. Might be useful to make it into a function since feldspars need to renamed
		#at multiple points due to recalling meemum.
		#Need all_elements, calc_phases. Return calc_phases.
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
	#Adding P2O5 saturation value and accuracy value of P2O5 saturation in melt
	#calc P_Sat
	if(any(rownames(calc_phases)=="Melt")){
		#calc P_Melt
		if(apatite_saturation=="Harrison & Watson 1984"){
			Melt_P <- correctApSatHW(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="H&W with Bea et al. 1992"){
			Melt_P <- correctApSatBea(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="H&W with Pichavant et al. 1992"){
			Melt_P <- correctApSatPich(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="Wolf & London 1994"){
			Melt_P <- correctApSatWolfLondon(Cmelt=calc_phases["Melt",major_elements])
		}
		#Multiply Melt_P/100 by Melt wt% to get P2O5 wt% of bulk to saturate Melt
		P_Melt <-calc_phases["Melt","wt%"]*Melt_P/100
	} else {
		final_P_Sat <- 0
		P_Melt <- 0
	}
	#should be slightly off due to coupling apatite and monazite saturation
	#normalisation may also affect the accuracy slightly
	aa <- matrix(NaN,nrow(calc_phases),2)
	colnames(aa) <- c("final_Sat","final_Sat accuracy")
	calc_phases <- cbind(calc_phases,aa)
	if(!is.na(any(match(rownames(calc_phases),"Ap")))){
		calc_phases["Ap","final_Sat"] <- signif(final_P_Sat,4)
		if(final_P_Sat == 0 && P_Melt == 0){
			calc_phases["Ap","final_Sat accuracy"] <- 100
		}else if(P_Melt != 0){
			calc_phases["Ap","final_Sat accuracy"] <- round(100-abs((P_Melt-final_P_Sat)/P_Melt*100),3)
		}
	}
	return(calc_phases)
}
correctApMnzSat <- function(c0,kd,temp,press,calc_phases,apatite_saturation,major_elements,mz_cL,min.props,cmins,melt.arg,dont) {
#Code for calculating distribution of P and LREE between melt/apatite/monazite
#Mod-tag: There is no subsolidus routine. For apatite alone, we assume that all P is in apatite. We probably can't assume this for Ap-Mnz routine.
#Mod-tag: isolating calcium makes the function un-detachable from Rcrust. Can make a more basic function that does not isolate calcium.
#First need to find Ap saturation value. For now, not going to isolate components to form Ap & Mnz to keep same equilibrium from start.
	#Not isolating or accounting for Ca-in-apatite has some effect on phase stabilities, mostly at points just above solidus. Otherwise it 
		#has a slight affect on plagioclase feldspar mode and composition.
	if(apatite_saturation=="Harrison & Watson 1984"){
		P_Sat <- correctApSatHW(Cmelt=calc_phases["Melt",major_elements],
			temp+273.15)
	} else if(apatite_saturation=="H&W with Bea et al. 1992"){
		P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt",major_elements],
			temp+273.15)
	} else if(apatite_saturation=="H&W with Pichavant et al. 1992"){
		P_Sat <- correctApSatPich(Cmelt=calc_phases["Melt",major_elements],
			temp+273.15)
	} else if(apatite_saturation=="Wolf & London 1994"){
		P_Sat <- correctApSatWolfLondon(Cmelt=calc_phases["Melt",major_elements])
	}
	#Multiply Melt_P/100 by Melt wt% to get P2O5 wt% of bulk to saturate Melt
	#P_Melt <-calc_phases["Melt","wt%"]*Melt_P/100
# browser()
#LREE are La, Ce, Pr, Nd, Pm, Sm, Gd
#P is worked with as P2O5 wt%
#LREE_Sat is the sum concentration of LREE in LREE saturated Melt (ppm)
	#LREE_Sat is first calculated from a monazite saturation equation, such as Montel, (1993), or Stepanov et al., (2012).
	#in run.Rcrust.R, monazite saturation function returns the trace concentration of the liquid, mz$cL.
	#The sum of LREE can be calculated from this variable, assuming we used Montel, (1993) and not Stepanov et al., (2012).
#P_Sat is the saturation concentration of phosphorus in melt.
	#P_Sat is calculated by apatite saturation equations such as Harrison & Watson, (1984) and others
	#Yakymchuck, (2017) uses Wolf and London (1994). Wolf and London (1994) produces similar results to Pichavant et al. (1992).
#Ap_LREE is a fixed calculation of the concentration of LREE in apatite (ppm)
	#To calculate Ap_LREE, a partition coefficient of LREE for apatite/melt is needed. D_ApMelt_LREE
	#Yakymchuck, (2017) uses 10, the average partition coefficient of Prowatke & Klemme, (2006), for La partitioning between apatite/melt with concentrations of SiO2.
	#It seems ideal to be able to enter a partition coefficient value then, and perhaps an option to automatically calculate it from the K-d file.
	#For granites the value would be an average of 11.9(1.9) & 8.0(1.7). D_ApMelt_LREE equal to 9.95
#LREE_Mnz is the sum of LREE in Mnz (ppm).
	#LREE_Mnz is intitially set to LREE_Bulk. After 4th step in iteration loop, the new value is passed onto next loop.
#LREE_Melt is the sum of LREE in Melt (ppm)
	#Set to LREE_Sat which is calculated from monazite saturation equations.
#P_Bulk is the total P2O5 wt% of bulk.
	#This is passed to function from c0 or trace elements.
#Mnz_P and Mnz_LREE are fixed stoichiometry of 29 wt.% P2O5 : 566794 ppm LREE in Mnz
#Ap_P is a fixed stoichiometry of P2O5 in apatite, 41 wt%.
	LREE_names <- c("Th", "La", "Ce", "Pr", "Nd", "Sm", "Gd")
	LREE_Sat <- sum(mz_cL[LREE_names])	#ppm
	LREE_Melt <- LREE_Sat	#ppm
	D_ApMelt_LREE <- 9.95
	Ap_LREE <- D_ApMelt_LREE*LREE_Melt	#ppm
	#LREE_Mnz is first set to the sum of LREE of the Bulk.
	LREE_Bulk <- sum(c0[LREE_names])	#ppm
	LREE_Mnz <- LREE_Bulk	#ppm
	P_Bulk<-c0[["P2O5"]]	#wt% (P2O5)
	Mnz_P <- 29	#wt%
	Mnz_LREE <- 566794	#ppm
	Ap_P <- 41	#wt%
	Ap_Ca <- 54 #wt%
# browser()
#If LREE_Mnz takes all LREE, exit iteration.
	cum_LREE_Mnz <- c()
	i <- 1
####
# browser()
##Mod-tag: this section is affected by LREE_Sat, calculated by Montel '93 equation.
	repeat{
		P_Mnz <- Mnz_P/Mnz_LREE * LREE_Mnz
		if(P_Mnz >= P_Bulk) {
			#Set P_Mnz to max P available
			P_Mnz <- P_Bulk
			#Reduce LREE_Mnz to match P_Mnz
			LREE_Mnz <- P_Mnz*Mnz_LREE/Mnz_P
		}
		P_Ap <- P_Bulk - P_Mnz
		LREE_Ap <- Ap_LREE / Ap_P * P_Ap
		LREE_Mnz <- LREE_Bulk - LREE_Ap
		cum_LREE_Mnz <- c(cum_LREE_Mnz, LREE_Mnz)
		if(i > 1 && abs((cum_LREE_Mnz[i]-cum_LREE_Mnz[i-1])/cum_LREE_Mnz[i]*100) < 0.01) {
			break
		}
		if(i > 5) {print("was >5"); print(abs((cum_LREE_Mnz[i]-cum_LREE_Mnz[i-1])/cum_LREE_Mnz[i]*100)); break}
		i <- i + 1
	}
# browser()
#P saturation is controlled by Ap. LREE saturation is controlled by Mnz.
	Ap_Diss <- P_Sat / P_Ap * calc_phases["Melt","wt%"]/100
	Mnz_Diss <- LREE_Sat / LREE_Mnz * calc_phases["Melt","wt%"]/100

	Ap_DissPct <- (P_Sat / (P_Ap + (Mnz_Diss * P_Mnz))) * calc_phases["Melt","wt%"]
	Mnz_DissPct <- (LREE_Sat / (LREE_Mnz + (Ap_Diss * LREE_Ap))) * calc_phases["Melt","wt%"]
# browser()
	if(any(rownames(calc_phases)=="Melt")){
		# last_major <- which(colnames(calc_phases)==major_elements[length(major_elements)])
		# Add P2O5 column to calc_phases
		# aa <- matrix(0,nrow(calc_phases),1)
		# colnames(aa) <- "P2O5"
		# calc_phases <- cbind(cbind(calc_phases[,1:last_major],aa),calc_phases[,(last_major+1):ncol(calc_phases)])
		Ap <- matrix(0, nrow = 1, ncol = 3)
		colnames(Ap) <- c("P2O5","LREE (ppm)","mass")
		rownames(Ap) <- "Ap"
		Ap[,"P2O5"] <- 41
		Ap[,"LREE (ppm)"] <- Ap_LREE
		Ap[,"mass"] <- P_Ap/41*100
		#Mod-tag: CaO = 54%, not included within phase but considered for total mass.
		#Mod-tag: Would have to isolate from c0 and redo calculations to include CaO as a component in Ap.
		Mnz <- matrix(0, nrow = 1, ncol = 3)
		colnames(Mnz) <- c("P2O5","LREE (ppm)","mass")
		rownames(Mnz) <- "Mnz"
		Mnz[,"P2O5"] <- Mnz_P
		Mnz[,"LREE (ppm)"] <- Mnz_LREE
		Mnz[,"mass"] <- P_Mnz/Mnz[,"P2O5"]*100
		# browser()
		#combined is used to calculate the proportions of Ap and Mnz for melt.
		combined <- rbind(Ap,Mnz)
		if (Ap_Diss >= 1){
			Ap_Diss <- 1
		}
		combined["Ap","mass"] <- combined["Ap","mass"]*Ap_Diss
		Ap[,"mass"] <- Ap[,"mass"] - combined["Ap","mass"]
		if (Mnz_Diss >= 1){
			Mnz_Diss <- 1
		}
		combined["Mnz","mass"] <- combined["Mnz","mass"]*Mnz_Diss
		Mnz[,"mass"] <- Mnz[,"mass"] - combined["Mnz","mass"]
		melt <- .wtd.add(combined)
		# colnames(melt) <- c("P2O5","LREE (ppm)","mass")
		rownames(melt) <- "melt"
		#Mod-tag: everything beyond this point could be moved into run.Rcrust.R. Then the function could be called and be independent of Rcrust.
			#Would need to return Ap, Mnz, melt as a list.
		#Formatting of calc_phases before trace elements are partitioned.
		# browser()
		br <- which(rownames(calc_phases)=="Bulk_rs")
		last_major <- which(colnames(calc_phases)==major_elements[length(major_elements)])
		#Add P2O5 column to calc_phases
		aa <- matrix(0,nrow(calc_phases),1)
		colnames(aa) <- "P2O5"
		calc_phases <- cbind(cbind(calc_phases[,1:last_major],aa),calc_phases[,(last_major+1):ncol(calc_phases)])
		calc_phases[1:br-1,"mass"] <- calc_phases[1:br-1,"mass"]/sum(calc_phases[1:br-1,"mass"]) * (100 - sum(Ap[,"mass"], Mnz[,"mass"], melt[,"mass"]))
		calc_phases["Bulk_rs",major_elements] <- calc_phases["Bulk_rs",major_elements]/sum(calc_phases["Bulk_rs",major_elements]) * (100 - c0["P2O5"])
		if(Ap_Diss != 1 || Ap[,"mass"] != 0){
			#Adds phase, apatite "Ap", to calc_phases.
			aa <- matrix(0,nrow = 1, ncol(calc_phases))
			colnames(aa) <- colnames(calc_phases)
			rownames(aa) <- "Ap"
			#Mod-tag: CaO 54 wt% is not isolated from bulk composition, remains in other phases. There is an empty mass for it in apatite
			# if(!is.na(any(match(major_elements,"CaO")))) { Ca <- "CaO" } else { Ca <- "CAO"}
			# aa[,Ca] <- 54
			aa[,"P2O5"] <- Ap[,"P2O5"]
			aa[,"mass"] <- Ap[,"mass"]
			br <- which(rownames(calc_phases)=="Bulk_rs")
			#there is an unaccounted-for mass of calcium that remains in phases other than apatite. The total mass of apatite includes the calcium portion and the mass of possible OH/Cl/F which can not be handled currently.
			calc_phases <- rbind(rbind(calc_phases[1:br-1,],aa),calc_phases[br:nrow(calc_phases),,drop=FALSE])
			last_major <- which(colnames(calc_phases)=="P2O5")
			calc_phases["Ap",(last_major+1):length(aa)]<-NaN
			#Average value of Ap density used to estimate vol%
			calc_phases["Ap","Density(kg/m3)"] <- 3190
			#mol of Ap is not a precise calculation as the density is estimated and the third site of F/Cl/OH is not included.
			calc_phases["Ap","mol"] <- calc_phases["Ap","mass"]/(54/100*56.0774 + 41/100*141.9445)
		}
		#Only add Mnz if it saturates
		if(Mnz_Diss != 1 || Mnz[,"mass"] != 0){
			#Adds phase, apatite "Mnz", to calc_phases.
			aa <- matrix(0,nrow = 1, ncol(calc_phases))
			colnames(aa) <- colnames(calc_phases)
			rownames(aa) <- "Mnz"
			aa[,"P2O5"] <- Mnz[,"P2O5"]
			aa[,"mass"] <- Mnz[,"mass"]
			#Mod-tag: LREE are partitioned using kd values. This does not necessarily equate to the assumed stoichiometry of 566794 ppm
			br <- which(rownames(calc_phases)=="Bulk_rs")
			calc_phases <- rbind(rbind(calc_phases[1:br-1,],aa),calc_phases[br:nrow(calc_phases),,drop=FALSE])
			calc_phases["Mnz",(last_major+1):length(aa)]<-NaN
			#not feasible to calculate mol of Mnz as it is largely made up of LREE. Will have insignificant effect on mol% of rock-forming phases.
			#Average value of Mnz density used to estimate vol%
			calc_phases["Mnz","Density(kg/m3)"] <- 5100
			calc_phases["Mnz","mol"] <- 0
		}
		aa <- matrix(0,nrow = 1, ncol(calc_phases))
		colnames(aa) <- colnames(calc_phases)
		rownames(aa) <- "melt"
		aa[,"P2O5"] <- melt[,"P2O5"]
		aa[,"mass"] <- melt[,"mass"]
		calc_phases["Melt",] <- .wtd.add(rbind(aa,calc_phases["Melt",]))
		calc_phases["Bulk_rs","P2O5"] <- c0["P2O5"]
		#Total masses of phases and components balance out if the empty mass of Ap, Mnz, melt is accounted for.
		#total mass + Ap[,"mass"]*0.59 + Mnz[,"mass"]*0.71 + melt[,"mass"]*(1-melt[,"P2O5"]/100)
		#formatting of output
		br <- which(rownames(calc_phases)=="Bulk_rs")
		#Fixing wt%
		calc_phases[1:br-1,"wt%"] <- calc_phases[1:br-1,"mass"]/sum(calc_phases[1:br-1,"mass"])*100
		calc_phases["Bulk_rs","wt%"] <- 100
		#Fixing vol%
		volume <- calc_phases[1:br-1,"mass"]/calc_phases[1:br-1,"Density(kg/m3)"]
		calc_phases[1:br-1,"vol%"] <- volume/sum(volume)*100
		calc_phases["Bulk_rs","vol%"] <- 100
		#Fixing mol%
		calc_phases[1:br-1,"mol%"] <- calc_phases[1:br-1,"mol"]/sum(calc_phases[1:br-1,"mol"])*100
		calc_phases["Bulk_rs","mol%"] <- 100
		calc_phases[1:br,1:which(colnames(calc_phases) == "mol")] <- round(calc_phases[1:br,1:which(colnames(calc_phases) == "mol")],4)
		min.props<-calc_phases[c(-which(rownames(calc_phases)=="Melt"),
									-which(rownames(calc_phases)=="Bulk_rs")),"wt%"]
		bpm <- BatchPM(kd = kd, c0=c0[1:length(c0)-1], pm = calc_phases["Melt","wt%"], min.props = min.props)
		# browser()
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
		#Mod-tag: How does LREE get trasferred to melt? Need to have specific LREE trace values to be able to transfer to melt
		#... Ap[LREE]
	}
	# print(paste0("Press ",press,"\n"))
	# print(paste0("Ap_Diss ",Ap_Diss,"\n"))
	# print(paste0("Mnz_Diss ",Mnz_Diss,"\n"))
	# print(paste0("Ap_DissPct ",Ap_DissPct,"\n"))
	# print(paste0("Mnz_DissPct ",Mnz_DissPct,"\n"))
	# return(paste0("Temp ",temp,"\n","Press ",press,"\n","Ap_Diss ",Ap_Diss,"\n","Mnz_Diss ",Mnz_Diss,"\n","Ap_DissPct ",Ap_DissPct,"\n","Mnz_DissPct ",Mnz_DissPct,"\n"))
	# printthis <- paste0("Temp ",temp,"\n","Press ",press,"\n","Ap_Diss ",Ap_Diss,"\n","Mnz_Diss ",Mnz_Diss,"\n","Ap_DissPct ",Ap_DissPct,"\n","Mnz_DissPct ",Mnz_DissPct,"\n")
	# return(printthis)
	#Mod-tag: It might be more useful to return an object other than calc_phases so that the function could be used universally and manipulated without being attached to Rcrust. Calc_phases could then be formatted in run.Rcrust.
	return(calc_phases)
}
#correctApMnzSatWithCa written by Sean Hoffman, 2022
correctApMnzSatWithCa <- function(c0,kd,temp,press,calc_phases,apatite_saturation,major_elements,Xmz) {
#Code for calculating distribution of P and LREE between melt/apatite/monazite, Ca to apatite.
#Mod-tag: There is no subsolidus routine. For apatite alone, we assume that all P is in apatite. We probably can't assume this for Ap-Mnz routine.
#Mod-tag: isolating calcium makes the function un-detachable from Rcrust. Can make a more basic function that does not isolate calcium.
	LREE_Sat <- correctMnzSatStepanov(press=press, temp=(temp+273.15), Cmelt=calc_phases["Melt",], Xmz=Xmz)
	if(apatite_saturation=="Harrison & Watson 1984"){
		P_Sat <- correctApSatHW(Cmelt=calc_phases["Melt",major_elements],
			temp+273.15)
	} else if(apatite_saturation=="H&W with Bea et al. 1992"){
		P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt",major_elements],
			temp+273.15)
	} else if(apatite_saturation=="H&W with Pichavant et al. 1992"){
		P_Sat <- correctApSatPich(Cmelt=calc_phases["Melt",major_elements],
			temp+273.15)
	} else if(apatite_saturation=="Wolf & London 1994"){
		P_Sat <- correctApSatWolfLondon(Cmelt=calc_phases["Melt",major_elements])
	}
#LREE are La, Ce, Pr, Nd, Pm, Sm, Gd
#P is worked with as P2O5 wt%
#LREE_Sat is the sum concentration of LREE in LREE saturated Melt (ppm)
	#LREE_Sat is first calculated from a monazite saturation equation, such as Montel, (1993), or Stepanov et al., (2012).
	#in run.Rcrust.R, monazite saturation function returns the trace concentration of the liquid, mz$cL.
	#The sum of LREE can be calculated from this variable, assuming we used Montel, (1993) and not Stepanov et al., (2012).
#P_Sat is the saturation concentration of phosphorus in melt. Could be named Melt_P, but kept to name used in Yakymchuck, (2017).
	#P_Sat is calculated by apatite saturation equations such as Harrison & Watson, (1984) and others
	#Yakymchuck, (2017) uses Wolf and London (1994). Wolf and London (1994) produces similar results to Pichavant et al. (1992).
#Ap_LREE is a fixed calculation of the concentration of LREE in apatite (ppm)
	#To calculate Ap_LREE, a partition coefficient of LREE for apatite/melt is needed. D_ApMelt_LREE
	#Yakymchuck, (2017) uses 10, the average partition coefficient of Prowatke & Klemme, (2006), for La partitioning between apatite/melt with concentrations of SiO2.
#LREE_Mnz is the sum of LREE in Mnz (ppm).
	#LREE_Mnz is intitially set to LREE_Bulk. After 4th step in iteration loop, the new value is passed onto next loop.
#LREE_Melt is the sum of LREE in Melt (ppm)
	#Set to LREE_Sat which is calculated from monazite saturation equations.
#P_Bulk is the total P2O5 wt% of bulk.
	#This is passed to function from c0 or trace elements.
#Mnz_P and Mnz_LREE are fixed stoichiometry of 29 wt.% P2O5 : 566794 ppm LREE in Mnz
#Ap_P is a fixed stoichiometry of P2O5 in apatite, 41 wt%.
	# LREE_names <- c("Th", "La", "Ce", "Pr", "Nd", "Sm", "Gd")
	LREE_names <- c("La", "Ce", "Pr", "Nd", "Sm")
	# LREE_Sat <- sum(mz_cL[LREE_names])	#ppm
	D_ApMelt_LREE <- 10
	#LREE_Mnz is first set to the sum of LREE of the Bulk.
	LREE_Bulk <- sum(c0[LREE_names])	#ppm
	LREE_Mnz <- LREE_Bulk	#ppm
	P_Bulk<-c0[["P2O5"]]	#wt% (P2O5)
	Mnz_P <- 29	#wt%
	Mnz_LREE <- 566794	#ppm
	Ap_P <- 41	#wt%
	Ap_Ca <- 54 #wt%
	#values used for solving for P_Sat when isolating CaO
	counter <- 1
	mid_P_Sat <- 0
	#loop that solves for P_Sat when isolating CaO
	repeat{
		LREE_Melt <- LREE_Sat	#ppm
		Ap_LREE <- D_ApMelt_LREE*LREE_Melt	#ppm
		cum_LREE_Mnz <- c()
		i <- 1
		#LREE and P mass components of Ap, Mnz are iterated for.
		repeat{
			P_Mnz <- Mnz_P/Mnz_LREE * LREE_Mnz
			if(P_Mnz >= P_Bulk) {
				#Set P_Mnz to max P available
				P_Mnz <- P_Bulk
				#Reduce LREE_Mnz to match P_Mnz
				LREE_Mnz <- P_Mnz*Mnz_LREE/Mnz_P
			}
			P_Ap <- P_Bulk - P_Mnz
			LREE_Ap <- Ap_LREE / Ap_P * P_Ap
			LREE_Mnz <- LREE_Bulk - LREE_Ap
			cum_LREE_Mnz <- c(cum_LREE_Mnz, LREE_Mnz)
			if(i > 1 && abs((cum_LREE_Mnz[i]-cum_LREE_Mnz[i-1])/cum_LREE_Mnz[i]*100) < 0.01) {
				break
			}
			if(i > 5) {print("iteration count >5"); print(abs((cum_LREE_Mnz[i]-cum_LREE_Mnz[i-1])/cum_LREE_Mnz[i]*100)); break}
			i <- i + 1
		}
		if(P_Sat<0 || mid_P_Sat<0){ P_Sat <- 0; mid_P_Sat <- 0}
		if(counter <= 2){
		#P saturation is controlled by Ap. LREE saturation is controlled by Mnz.
			Ap_Diss <- P_Sat / P_Ap * calc_phases["Melt","wt%"]/100
			Mnz_Diss <- LREE_Sat / LREE_Mnz * calc_phases["Melt","wt%"]/100
		#Step is repeated considering the amount of LREE contributed by the breakdown of apatite and the amount of phosphorus contributed by the breakdown of monazite.
			Ap_DissPct <- (P_Sat / (P_Ap + (Mnz_Diss * P_Mnz))) * calc_phases["Melt","wt%"]
			Mnz_DissPct <- (LREE_Sat / (LREE_Mnz + (Ap_Diss * LREE_Ap))) * calc_phases["Melt","wt%"]
		} else {
		#P saturation is controlled by Ap. LREE saturation is controlled by Mnz.
			Ap_Diss <- mid_P_Sat / P_Ap * calc_phases["Melt","wt%"]/100
			Mnz_Diss <- LREE_Sat / LREE_Mnz * calc_phases["Melt","wt%"]/100
		#Step is repeated considering the amount of LREE contributed by the breakdown of apatite and the amount of phosphorus contributed by the breakdown of monazite.
			Ap_DissPct <- (mid_P_Sat / (P_Ap + (Mnz_Diss * P_Mnz))) * calc_phases["Melt","wt%"]
			Mnz_DissPct <- (LREE_Sat / (LREE_Mnz + (Ap_Diss * LREE_Ap))) * calc_phases["Melt","wt%"]
		}
		Ca_Ap <- P_Ap / Ap_P * Ap_Ca
	#The amount of Ca that would be dissolved can be excluded from isolation, can go back into use for phase equilibria calculations
		if(Ap_DissPct >= 100){ Ap_DissPct <- 100}
		isolate_Ca_Ap <- Ca_Ap * (100-Ap_DissPct)/100
		Ca <- which(toupper(names(c0))=="CAO")
		if(isolate_Ca_Ap <= c0[Ca]){
			#removing portions of components from c0, which are reserved for apatite, Mnz and Melt
			#then recalculating calc_phases
			store_c0 <- c0
			c0[Ca] <- c0[Ca]-isolate_Ca_Ap
			c0["mass"] <- c0["mass"]-sum(isolate_Ca_Ap,c0["P2O5"])
			old_calc_phases <- calc_phases
			calc_phases <- try(wrapper(comps,c(c0[1:(length(major_elements))],c0[length(c0)]),press,temp,calc_choice),silent=TRUE)
		} else {
			cat("Error: Not enough calcium to form apatite\n")
			stop()
		}
		if(counter <= 2){
			Old_P_Sat <- P_Sat
		} else {
			Old_P_Sat <- mid_P_Sat
		}
		Old_LREE_Sat <- LREE_Sat
		if(any(rownames(calc_phases)=="Melt")){
			#calc P_Melt
			if(apatite_saturation=="Harrison & Watson 1984"){
				P_Sat <- correctApSatHW(Cmelt=calc_phases["Melt",major_elements],
					temp+273.15)
			} else if(apatite_saturation=="H&W with Bea et al. 1992"){
				P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt",major_elements],
					temp+273.15)
			} else if(apatite_saturation=="H&W with Pichavant et al. 1992"){
				P_Sat <- correctApSatPich(Cmelt=calc_phases["Melt",major_elements],
					temp+273.15)
			} else if(apatite_saturation=="Wolf & London 1994"){
				P_Sat <- correctApSatWolfLondon(Cmelt=calc_phases["Melt",major_elements])
			}
			LREE_Sat <- correctMnzSatStepanov(press=press, temp=(temp+273.15), Cmelt=calc_phases["Melt",], Xmz=Xmz)
		} else {
			#melt is not a stable phase after recalculation, return old_calc_phases
			cat("melt is not a stable phase after isolating ",isolate_Ca_Ap," CaO, skipped Ap_Mnz saturation calculation.\n")
			return(old_calc_phases)
		}
		#P saturation is controlled by Ap. LREE saturation is controlled by Mnz.
		Ap_Diss <- P_Sat / P_Ap * calc_phases["Melt","wt%"]/100
		Mnz_Diss <- LREE_Sat / LREE_Mnz * calc_phases["Melt","wt%"]/100
		Ap_DissPct <- (P_Sat / (P_Ap + (Mnz_Diss * P_Mnz))) * calc_phases["Melt","wt%"]
		Mnz_DissPct <- (LREE_Sat / (LREE_Mnz + (Ap_Diss * LREE_Ap))) * calc_phases["Melt","wt%"]
		#Change itn_threshold value to set iteration threshold (%)
		itn_threshold <- 1
		#counter > 2 uses midpoint value as calculation.
		if(counter > 2){
			if(isTRUE(all.equal(P_Sat,mid_P_Sat)) || (abs((P_Sat - mid_P_Sat) / P_Sat * 100) < itn_threshold)){
				c0 <- store_c0 ; break
			}
			if(counter >4){
				print(paste0("Did not resolve to itn_theshold (",itn_threshold,
					"%) in 5 iterations. ApSat iteration %change ",abs((P_Sat - mid_P_Sat) / P_Sat * 100)))
				c0 <- store_c0 ; break
			}
			if(counter > 1){
				#if counter is on 2 or more and has not been solved to satisfactory itn_threshold, iteration needs to be assisted.
				#Take last 2 P_Sat values and find halfway point.
				mid_P_Sat <- (P_Sat + Old_P_Sat) / 2
			}
		}else{
			if(isTRUE(all.equal(P_Sat,Old_P_Sat)) || (abs((P_Sat - Old_P_Sat) / P_Sat * 100) < itn_threshold)){
				c0 <- store_c0
				break
			} else if(counter > 1){
				#if counter is on 2 or more and has not been solved to satisfactory itn_threshold, iteration needs to be assisted.
				#Take last 2 P_Sat values and find halfway point.
				mid_P_Sat <- (P_Sat + Old_P_Sat) / 2
			}
		}
		c0 <- store_c0
		counter <- counter + 1
	}
	#final_LREE_Sat set to the value used in recalculation of stable phases
	final_LREE_Sat<-Old_LREE_Sat
	#final_P_Sat set to the value used to isolate CaO and recalculate stable phases
	if(counter > 3){ final_P_Sat <- mid_P_Sat } else { final_P_Sat <- Old_P_Sat }
	#This section partitions components to melt and/or adds apatite & monazite to calc_phases and formats it
	store_c0 <- c0
	if(any(rownames(calc_phases)=="Melt")){
		Ap <- matrix(0, nrow = 1, ncol = 3)
		colnames(Ap) <- c("P2O5","LREE (ppm)","mass")
		rownames(Ap) <- "Ap"
		Ap[,"P2O5"] <- 41
		Ap[,"LREE (ppm)"] <- Ap_LREE
		Ap[,"mass"] <- P_Ap/41*100
		Mnz <- matrix(0, nrow = 1, ncol = 3)
		colnames(Mnz) <- c("P2O5","LREE (ppm)","mass")
		rownames(Mnz) <- "Mnz"
		Mnz[,"P2O5"] <- Mnz_P
		Mnz[,"LREE (ppm)"] <- Mnz_LREE
		Mnz[,"mass"] <- P_Mnz/Mnz[,"P2O5"]*100
		#combined is used to calculate the proportions of Ap and Mnz for melt.
		combined <- rbind(Ap,Mnz)
		if (Ap_DissPct >= 100){ Ap_DissPct <- 100 }
		combined["Ap","mass"] <- combined["Ap","mass"] * Ap_DissPct/100
		if(Ap_DissPct == 100){ Ap[,"mass"] <- 0 } else {
		Ap[,"mass"] <- Ap[,"mass"] - combined["Ap","mass"] }
		if (Mnz_DissPct >= 100){ Mnz_DissPct <- 100 }
		combined["Mnz","mass"] <- combined["Mnz","mass"]*Mnz_DissPct/100
		if(Mnz_DissPct == 100) { Mnz[,"mass"] <- 0 } else {
		Mnz[,"mass"] <- Mnz[,"mass"] - combined["Mnz","mass"] }
		melt <- .wtd.add(combined)
		rownames(melt) <- "melt"
		#Mod-tag: everything beyond this point could be moved into run.Rcrust.R or a trace element partitioning function
		#Then the function could be called and be independent of Rcrust.
			#Would need to return Ap, Mnz, melt as a list.
		#Formatting of calc_phases before trace elements are partitioned.
		br <- which(rownames(calc_phases)=="Bulk_rs")
		last_major <- which(colnames(calc_phases)==major_elements[length(major_elements)])
		#Add P2O5 column to calc_phases
		aa <- matrix(0,nrow(calc_phases),1)
		colnames(aa) <- "P2O5"
		calc_phases <- cbind(cbind(calc_phases[,1:last_major],aa),calc_phases[,(last_major+1):ncol(calc_phases)])
		calc_phases[1:br-1,"mass"] <- calc_phases[1:br-1,"mass"]/sum(calc_phases[1:br-1,"mass"]) * (calc_phases["Bulk_rs","mass"] - sum(Ap[,"mass"]*0.05, Mnz[,"mass"]*0.71, melt[,"mass"]*(100-melt[,"P2O5"])/100))
		calc_phases["Bulk_rs",major_elements] <- calc_phases["Bulk_rs",major_elements]/sum(calc_phases["Bulk_rs",major_elements]) * (100 - c0["P2O5"] - isolate_Ca_Ap)
		if(Ap_DissPct != 100 || Ap[,"mass"] != 0){
			#Adds phase, apatite "Ap", to calc_phases.
			aa <- matrix(0,nrow = 1, ncol(calc_phases))
			colnames(aa) <- colnames(calc_phases)
			rownames(aa) <- "Ap"
			if(!is.na(any(match(major_elements,"CaO")))) { Ca <- "CaO" } else { Ca <- "CAO"}
			aa[,Ca] <- 54
			aa[,"P2O5"] <- Ap[,"P2O5"]
			aa[,"mass"] <- Ap[,"mass"]
			br <- which(rownames(calc_phases)=="Bulk_rs")
			calc_phases <- rbind(rbind(calc_phases[1:br-1,],aa),calc_phases[br:nrow(calc_phases),,drop=FALSE])
			last_major <- which(colnames(calc_phases)=="P2O5")
			calc_phases["Ap",(last_major+1):length(aa)]<-NaN
			#Average value of Ap density used to estimate vol%
			calc_phases["Ap","Density(kg/m3)"] <- 3190
			#mol of Ap is not a precise calculation as the density is estimated and the third site of F/Cl/OH is not included.
			calc_phases["Ap","mol"] <- calc_phases["Ap","mass"]/(54/100*56.0774 + 41/100*141.9445)
			#must re-add the isolated calcium to bulk_rs rownames
			calc_phases["Bulk_rs",Ca] <- calc_phases["Bulk_rs",Ca] + isolate_Ca_Ap
			calc_phases["Bulk_rs","P2O5"] <- calc_phases["Bulk_rs","P2O5"] + (calc_phases["Ap","P2O5"]/100)*calc_phases["Ap","mass"]
		}
		#Only add Mnz if it saturates
		if(Mnz_DissPct != 100 || Mnz[,"mass"] != 0){
			#Adds phase, apatite "Mnz", to calc_phases.
			aa <- matrix(0,nrow = 1, ncol(calc_phases))
			colnames(aa) <- colnames(calc_phases)
			rownames(aa) <- "Mnz"
			aa[,"P2O5"] <- Mnz[,"P2O5"]
			aa[,"mass"] <- Mnz[,"mass"]
			#Mod-tag: LREE are partitioned using kd values. This does not necessarily equate to the assumed stoichiometry of 566794 ppm
			br <- which(rownames(calc_phases)=="Bulk_rs")
			calc_phases <- rbind(rbind(calc_phases[1:br-1,],aa),calc_phases[br:nrow(calc_phases),,drop=FALSE])
			calc_phases["Mnz",(last_major+1):length(aa)]<-NaN
			#not feasible to calculate mol of Mnz as it is largely made up of LREE. Will have insignificant effect on mol% of rock-forming phases.
			#Average value of Mnz density used to estimate vol%
			calc_phases["Mnz","Density(kg/m3)"] <- 5100
			calc_phases["Mnz","mol"] <- 0
			calc_phases["Bulk_rs","P2O5"] <- calc_phases["Bulk_rs","P2O5"] + (calc_phases["Mnz","P2O5"]/100)*calc_phases["Mnz","mass"]
		}
		calc_phases["Melt","mass"] <- calc_phases["Melt","mass"] + melt[,"mass"]
		calc_phases["Melt","P2O5"] <- ((melt[,"mass"] * (melt[,"P2O5"]/100)) / calc_phases["Melt","mass"]) * 100
		calc_phases["Bulk_rs","P2O5"] <- calc_phases["Bulk_rs","P2O5"] + ((calc_phases["Melt","P2O5"]/100) * calc_phases["Melt","mass"])
		#Total masses of phases and components balance out if the empty mass of Ap, Mnz is accounted for.
		#total mass + calc_phases["Ap","mass"]*0.05 + calc_phases["Mnz","mass"]*0.71
		#formatting of output
		br <- which(rownames(calc_phases)=="Bulk_rs")
		calc_phases[1:br-1,"mass"] <- calc_phases[1:br-1,"mass"]/sum(calc_phases[1:br-1,"mass"])*100
		calc_phases["Bulk_rs","mass"] <- 100
		#Fixing wt%
		calc_phases[1:br-1,"wt%"] <- calc_phases[1:br-1,"mass"]/sum(calc_phases[1:br-1,"mass"])*100
		calc_phases["Bulk_rs","wt%"] <- 100
		#Fixing vol%
		volume <- calc_phases[1:br-1,"mass"]/calc_phases[1:br-1,"Density(kg/m3)"]
		calc_phases[1:br-1,"vol%"] <- volume/sum(volume)*100
		calc_phases["Bulk_rs","vol%"] <- 100
		#Fixing mol%
		calc_phases[1:br-1,"mol%"] <- calc_phases[1:br-1,"mol"]/sum(calc_phases[1:br-1,"mol"])*100
		calc_phases["Bulk_rs","mol%"] <- 100
		calc_phases[1:br,1:which(colnames(calc_phases) == "mol")] <- round(calc_phases[1:br,1:which(colnames(calc_phases) == "mol")],4)
		#Mod-tag: renaming feldspars. Might be useful to make it into a function since feldspars need to renamed
		#at multiple points due to recalling meemum.
		#Need all_elements, calc_phases. Return calc_phases.
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
		if(!is.na(any(match(rownames(calc_phases),"H2O")))) {
			min.props <- calc_phases[c(-which(rownames(calc_phases)=="Melt"),
										-which(rownames(calc_phases)=="Bulk_rs"),
										-which(rownames(calc_phases)=="H2O")),"wt%"]
		} else {
			min.props <- calc_phases[c(-which(rownames(calc_phases)=="Melt"),
									-which(rownames(calc_phases)=="Bulk_rs")),"wt%"]
		}
		bpm <- BatchPM(kd = kd, c0=c0[1:length(c0)-1], pm = calc_phases["Melt","wt%"], min.props = min.props)
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
		#Adding P2O5 saturation value and accuracy value of P2O5 saturation in melt
		#calc P_Sat
		if(apatite_saturation=="Harrison & Watson 1984"){
			P_Sat <- correctApSatHW(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="H&W with Bea et al. 1992"){
			P_Sat <- correctApSatBea(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="H&W with Pichavant et al. 1992"){
			P_Sat <- correctApSatPich(Cmelt=calc_phases["Melt",major_elements],
				temp+273.15)
		} else if(apatite_saturation=="Wolf & London 1994"){
			P_Sat <- correctApSatWolfLondon(Cmelt=calc_phases["Melt",major_elements])
		}
		LREE_Sat <- correctMnzSatStepanov(press=press, temp=(temp+273.15), Cmelt=calc_phases["Melt",], Xmz=Xmz)
		#P_Sat and LREE_Sat accuracy
		#should be slightly off due to coupling apatite and monazite saturation
		#normalisation may also affect the accuracy slightly
		aa <- matrix(NaN,nrow(calc_phases),2)
		colnames(aa) <- c("final_Sat","final_Sat accuracy")
		calc_phases <- cbind(calc_phases,aa)
		if(!is.na(any(match(rownames(calc_phases),"Ap")))){
			calc_phases["Ap","final_Sat"] <- signif(final_P_Sat,4)
			if(final_P_Sat == 0 && P_Sat == 0){
				calc_phases["Ap","final_Sat accuracy"] <- 100
			}else if(P_Sat != 0){
				calc_phases["Ap","final_Sat accuracy"] <- round(100-abs((P_Sat-final_P_Sat)/P_Sat*100),3)
			}
		}
		if(!is.na(any(match(rownames(calc_phases),"Mnz")))){
			calc_phases["Mnz","final_Sat"] <- signif(final_LREE_Sat,4)
			if(final_LREE_Sat == 0 && LREE_Sat == 0){
				calc_phases["Mnz","final_Sat accuracy"] <- 100
			}else if(LREE_Sat != 0){
				calc_phases["Mnz","final_Sat accuracy"] <- round(100-abs((LREE_Sat-final_LREE_Sat)/LREE_Sat*100),3)
			}
		}
	}
	#Mod-tag: It might be more useful to return an object other than calc_phases so that the function could be used universally and manipulated without being attached to Rcrust. Calc_phases could then be formatted in run.Rcrust.
	return(calc_phases)
}

#Mod-tag: write a general function to add phases to calc_phases
#return a formatted calc_phases.