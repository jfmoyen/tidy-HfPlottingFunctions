if(length(grep("^176Hf/177Hf",colnames(WR)))>0){
# Constants ##########################
#lambda<-1.93e-11) # 176Lu Blichert-Toft and Albarede 1997

lambda<-1.867e-11 # 176Lu  Söderlund et al. 2004 # used by Yingde

#lambda<-1.865e-11  # 176Lu Scherer et al. (2001) - used by CERCAMS 
names(lambda)<-c("Hf")
e<-2.718281828 #exp(1)

# DM - present-day Hf isotopic composition

#DM<-c(0.0384,0.283251) # based on 1.93e-11, 0.28323 Belousova et al. 2010 and references therein
DM<-c(0.0384,0.28325) # CERCAMS: assumes 0.279718 176/177 initial 4.56 Ga ago, Yingde uses the same (Griffin et al. 2000) 
names(DM)<-c("176Lu/177Hf","176Hf/177Hf")
DM<-t(as.matrix(DM))
rownames(DM)<-"DM"

# CHUR - present-day Hf isotopic composition
#CHUR<-c(0.0336,0.282785) # Bouvier et al. (2008)
CHUR<-c(0.0332,0.282772) # CERCAMS: based on 1.93e-11!!!, values of Blichert-Toft et al. (1997), Yingde uses the same 
names(CHUR)<-c("176Lu/177Hf","176Hf/177Hf")
CHUR<-t(as.matrix(CHUR))
rownames(CHUR)<-"CHUR"

# 
# 176Lu/177Hf ratio of the source is unknown and thus needs to be assumed

# 176Lu/177Hf of average continetal crust
#RCC<-0.015 # 176Lu/177Hf in Bulk CC; Rudnick Gao 2003 nebo Griffin et al., 2004 in Belousova 2010 
RCC<-0.022 # basaltic 176Lu/177Hf ratio: Lancaster et al. 2011
#RCC<-0.0093 # upper continetal crust, Amelin et al. 1999


# Important functions #
########################################################################

# Calculates Hf initials ****** 
Hfinitial<-function(x,age){
    R<-x[,"176Lu/177Hf",drop=FALSE]
    I<-x[,"176Hf/177Hf",drop=FALSE]
    y<-I-(R*(e^(lambda["Hf"]*age*10^6)-1))
    age[is.na(age)]<--1
    y[age==0]<-I
    age[age==-1]<-NA
    return(y)
}

# Calculates initial Eps Hf values ****** 
Hfepsilon<-function(x,age){
    if(any(colnames(x)=="176Hf/177Hfi")){
        ee<-x[,"176Hf/177Hfi",drop=FALSE]
    }else{   
        ee<-Hfinitial(x,age)
    }
    X<-round((ee/Hfinitial(CHUR,age)-1)*10^4,2)
    return(X)
}

# Calculates single-stage Hf model ages ****** 
HfDMage<-function(x){
    #R<-x[,"176Lu/177Hf"]
    #I<-x[,"176Hf/177Hf"]
    X<-round(1/lambda["Hf"]*log(((I-DM[,"176Hf/177Hf"])/(R-DM[,"176Lu/177Hf"]))+1)/10^9,3)
    X
}

# Calculates THf two/stage DM model ages TO BE CHECKED
HfDM2stgAge<-function(I,age){
    #R<-x[,"176Lu/177Hf"]
    #I<-x[,"176Hf/177Hfi"]
    #age<-x[,"Age"]
    #citatel<-I-(exp(lambda["Hf"]*age*1e6)-1)*(R-RCC)-DM[,"176Hf/177Hf"]
    #jmenovatel<-RCC-DM[,"176Lu/177Hf"]
    
    X<-round(1/lambda["Hf"]*log(((I-Hfinitial(DM,age))/(RCC-DM[,"176Lu/177Hf"]))+1)/10^9,3)
    #cbind(WR[,"TDM"],HfDM2stgAge()-Hfinit[,"Age (Ma)"])
    return(X+age/1e3)
}


Hfiso<-function(age=NULL){
    if(length(grep("^176Hf/177Hf",colnames(WR)))==0)return()
    on.exit(options("show.error.messages"=TRUE))
    
    prd<-which(toupper(colnames(labels))=="SAMPLE")
    if(length(prd)>0){
        sample<-labels[,prd[1]]
    }else{
        sample<-rep(NA,nrow(labels))
    }
    
    if(!any(toupper(colnames(WR))=="AGE")){
        if(is.null(age)){
            age<-winDialogString("Age (Ma)",as.character(340))
            if(is.null(age)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
            age<-as.numeric(age)
            age<-rep(age,times=nrow(WR))
        }
    }else{
        age<-WR[,toupper(colnames(WR))=="AGE"]
    }        
    
    if(any(colnames(WR)=="176Hf/177Hfi")){ # Use initial Hf isotopic ratios if present in te original file
        Hfinit<-WR[,"176Hf/177Hfi"]
    }else{   
        Hfinit<-Hfinitial(WR,age)
    }
    
    eps<-Hfepsilon(WR,age)
    #TDM<-HfDMage(WR)
    
    TDM2stg<-HfDM2stgAge(Hfinit,age) # Use the two-stage CC age
    Hfinit<-data.frame(age,as.character(sample),round(Hfinit,6),eps,TDM2stg,check.names=FALSE)
    
    rownames(Hfinit)<-rownames(WR)
    #colnames(init)<-c("Age (Ma)","176Hf/177Hfi","EpsHfi","HfTDM","HfTDM.2stg")
    colnames(Hfinit)<-c("Age (Ma)","Sample","176Hf/177Hfi","EpsHfi","HfTDM.2stg")
     
    ee<-filterOut(Hfinit,colnames(Hfinit),n=6)
    if(!getOption("gcd.shut.up")){cat("\n");print(ee)}
    
    # Set up labels for plotting with correct formatting
    ii<-unique(age)
    if(length(ii)==1){
        #b<-paste("epsilon*Hf[",ii,"]",sep="")
        #b<-as.expression(b)
        #Hfepslabi<<-parse(text=b)
        b<-paste("epsilon[",ii,"]^Hf",sep="")
        b<-as.expression(b)
        epsHflabi<<-parse(text=b)
        
        b<-paste("\" \"^176*Hf/\" \"^177*Hf[",ii,"]",sep="")
        b<-as.expression(b)
        hflab<<-parse(text=b)

    }else{
        epsHflabi<<-expression(epsilon[i]^Hf)
        hflab<<-expression(" "^176*Hf/" "^177*Hf[i])
    }
    epsHflab<<-expression(epsilon[Hf])
    age<<-age
    return(Hfinit)
}
Hfinit<-Hfiso()

winMenuAdd("Plugins/Hf isotopes")
winMenuAddItem("Plugins/Hf isotopes","Recalculate to new age","Hfinit<<-Hfiso()")
winMenuAddItem("Plugins/Hf isotopes","Save Hf isotopic data","ee<-HfsaveResultsIso(digits=6)")
winMenuAddItem("Plugins/Hf isotopes","Append Hf isotopic data","ee<-HfaddResultsIso()")
winMenuAddItem("Plugins/Hf isotopes","------------------------------------    ","none")
winMenuAddItem("Plugins/Hf isotopes","Hf growth lines","ee<-HfageEps()")
winMenuAddItem("Plugins/Hf isotopes","------------------------------------        ","none")
winMenuAddItem("Plugins/Hf isotopes","Boxplot isotopic ratios/model ages","boxplotIso()")
winMenuAddItem("Plugins/Hf isotopes","Stripplot isotopic ratios/model ages","stripplotIso()")

# Age vs EpsNdi GENERAL ENTRY POINT ##############################################################
HfageEps<-function(){    
    on.exit(options("show.error.messages"=TRUE))
    ee<-select.list(c("Hf two stage"),preselect="Hf two stage",multiple=FALSE) 
    if(nchar(ee)==0){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
    
    if(ee=="Hf two stage"){
        ee<-HfageEps2();return()
    }
    # Only two-stage model available so far
    ##### Single stage model still TODO !!!!!!!!!!!!!!!!!!!! #####
}


# Age vs EpsNdi 2 stg ######################################################################
HfageEps2<-function(which=rownames(Hfinit),xmin=0,ymax=20,ymin=-15,xmax=3,pch=labels[which,"Symbol"],col=labels[which,"Colour"],cex=labels[which,"Size"],new=TRUE,main="Hf isotopic growth diagram (two-stage)",...){
    if(length(pch)==1&length(which)>1){
        pch<-rep(pch,length(which))
    }
    names(pch)<-which
    
    if(length(col)==1&length(which)>1){
        col<-rep(col,length(which))
    }
    names(col)<-which
    
    if(length(cex)==1&length(which)>1){
        cex<-rep(cex,length(which))
    }
    names(cex)<-which
    
    x<-age[which]/1000
    names(x)<-rownames(Hfinit[which,])
    y<-Hfinit[which,"EpsHfi"]
    names(y)<-rownames(Hfinit[which,])
    plotWithLimits(x,y,xmin=xmin,ymax=ymax,ymin=ymin,xmax=xmax,xlab="Age (Ga)",ylab=epsHflab,col=col,pch=pch,cex=cex,new=new,main=main,...)
    sheet$demo$template$abline1<-list("abline",h=0,col="darkgrey",lwd=3)
    
    Q<--CHUR[,"176Lu/177Hf"]*lambda["Hf"]/CHUR[,"176Hf/177Hf"]*1e4*1e9
    
    # DM
    sheet$demo$template$abline2<-list("abline",a=(DM[,"176Hf/177Hf"]/CHUR[,"176Hf/177Hf"]-1)*1e4,b=Q*(DM[,"176Lu/177Hf"]-CHUR[,"176Lu/177Hf"])/CHUR[,"176Lu/177Hf"],col="darkgrey",lwd=3)
    
    xmax<-sheet$demo$call$xlim[2]
    sheet$demo$template$text1<-list("text",x=xmax-xmax/20,y=0.5,text="CHUR",size=1.5)
    sheet$demo$template$text2<-list("text",x=xmax-xmax/20,y=15-1/.456*xmax-1,text="DM",size=1.5)
    
    # Hf isotopic development lines for our samples
    temp.lines<-list()

    for (i in 1:length(which)){
            if(!is.na(Hfinit[which[i],"176Hf/177Hfi"])){
                x<-c(Hfinit[which[i],"Age (Ma)"]/1e3,Hfinit[which[i],"HfTDM.2stg"])
                y<-c(Hfinit[which[i],"EpsHfi"],Hfepsilon(DM,Hfinit[which[i],"HfTDM.2stg"]*10^3))
                ee<-eval(parse(text=paste("temp.lines$lines",i,"<-list(\"lines\",x=c(", paste(x,collapse=","),"),y=c(", paste(y,collapse=","),"),col= \"",col[which[i]],"\",lty=\"solid\")",sep="")))
            }
    }
    
    sheet$demo$call$col<-col#labels[names(x.data),"Colour"]
    sheet$demo$call$pch<-pch#labels[names(x.data),"Symbol"]
    sheet$demo$call$cex<-cex#labels[names(x.data),"Size"]
    sheet$demo$template<-c(sheet$demo$template,temp.lines)
    sheet$demo$template$rug1<-list("rug",x=Hfinit[which,"HfTDM.2stg"],ticksize=-0.01,side=1,lwd=1,col="darkred")
    sheet$demo$template$rug2<-list("rug",x=Hfinit[which,"EpsHfi"],ticksize=-0.01,side=2,lwd=1,col="darkred")
    assign("sheet",sheet,.GlobalEnv)
    
    pp<<-figaro(demo,prefix="sheet")
    assign("pp",pp,.GlobalEnv)
    
    #sheet<<-sheet
    #pp<<-figaro(demo,prefix="sheet")
    figRedraw()
    if(new){
        figaroOn()
    }
    invisible()
}

.HfselectVariable<-function(){
    on.exit(options("show.error.messages"=TRUE))
    vars<-c("176Hf/177Hfi","EpsHfi","2 stg DM model ages")
      
    what<-select.list(vars,preselect="EpsHfi",multiple=FALSE)
    if(nchar(what)==0){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
    
    ee<-c("176Hf/177Hfi","EpsHfi","HfTDM.2stg")[which(vars==what)]
    return(c(ee,what))
}

boxplotIso<-function(){
    what<-.HfselectVariable()
    names(groups)<-rownames(WR)
    if(substr(what,1,5)!="delta"){
        x<-factor(groups[rownames(Hfinit)],ordered=TRUE)
    }else{
        x<-factor(groups[rownames(WR)],ordered=TRUE)
    }
    windows(width = 5, height = 5, pointsize = 10,title=paste("Boxplot of ",what[2]))
    par(las=1)
    par(oma=c(1,5,1,1))
    
    ee<-grep("model",what[2])
    if(length(ee)>0){
        xlab<-paste(what[2]," (Ga)",sep="")
    }else{
        xlab<-annotate(what[2])
    }
    
    if(substr(what,1,5)!="delta"){
        boxplot(Hfinit[,what[1]]~x,ylab="",xlab=xlab,varwidth=TRUE,cex.axis=0.7,col="lightblue",horizontal=TRUE,xaxs="i")
        ee<-Hfinit[,what[1]]
    }else{
        boxplot(WR[,what[1]]~x,ylab="",xlab=xlab,varwidth=TRUE,cex.axis=0.7,col="lightblue",horizontal=TRUE,xaxs="i")
        ee<-WR[,what[1]]
    }
    ee<-ee[!is.na(ee)]
    i<-pretty(ee)
    
    abline(v=i,col="gray",lty="dotted")
    abline(v=0,lty="dashed")
    abline(h=(0.5:35.5),lty="dotted",col="gray")
    box()
    figaroOff()
}

stripplotIso<-function(){
    strip.main<-function(elem,xlab="",title="Stripplot"){    
        names(groups)<-rownames(WR)
        if(substr(elem,1,5)!="delta"){
            ee<-Hfinit[,elem,drop=FALSE]
        }else{
            ee<-WR[,elem,drop=FALSE]
        }
        ii<-!is.na(ee)
        which<-rownames(ee)[ii]
        
        ee<-ee[which,,drop=FALSE]
        x<-factor(groups[which],ordered=TRUE)
        trellis.device(bg = "white",new = TRUE,retain = FALSE,title=title)
        trellis.par.set("background",list(col="white"))
        stripplot(x ~ ee,aspect = 1,jitter = TRUE, pch=labels[which,"Symbol"],col=labels[which,"Colour"],xlab=xlab,cex=labels[which,"Size"])
    }
    
    what<-.HfselectVariable()   
    ee<-grep("model",what[2])
    if(length(ee)>0){
        xlab<-paste(what[2]," (Ga)",sep="")
    }else{
        xlab<-annotate(what[2])
    }
     figaroOff()
    strip.main(what[1],xlab=xlab,title=paste("Stripplot of ",what[2],sep=""))
}


HfaddResultsIso<-function(){
    if(is.na(all( match(colnames(Hfinit),colnames(WR))))){
        cat("Appending...\n")
        addResults("Hfinit")
    }else{
        cat("Replacing...\n")
        ee<-WR
        ee[rownames(Hfinit),colnames(Hfinit)]<-Hfinit
        assign("WR",ee,.GlobalEnv)
    }    
}

HfsaveResultsIso<-function(digits = 6){
    saveResults(Hfinit,digits=digits)
}

}else{
    if(!getOption("gcd.shut.up"))cat("WARNING: Skipping: No Hf isotopic data present!")
}
