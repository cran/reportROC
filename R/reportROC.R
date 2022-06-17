reportROC=function(gold,
                   predictor=NULL,
                   predictor.binary=NULL,
                   important="se",
                   positive="l",
                   exact=NULL,
                   plot=TRUE,xlab="1-Specificity",ylab="Sensitivity"){

  ### if the predictor is continuous

  if(is.null(predictor.binary)){
    data=data.frame(gold,predictor)
    nrow1=nrow(data)
    data=na.omit(data)
    nrow2=nrow(data)
    if(nrow1!=nrow2){
      message(paste(nrow1-nrow2,"sample(s) with missing data was removed!"))
    }
    gold=data$gold
    predictor=data$predictor
  }

  if(is.null(predictor)){
    data=data.frame(gold,predictor.binary)
    nrow1=nrow(data)
    data=na.omit(data)
    nrow2=nrow(data)
    if(nrow1!=nrow2){
      message(paste(nrow1-nrow2,"sample(s) with missing data was removed!"))
    }
    gold=data$gold
    predictor.binary=data$predictor.binary
  }

  if(length(unique(gold))==2){
    error=FALSE
    table.gold=2-as.numeric(as.factor(gold))#'0' is positive

    if(!is.null(predictor) & is.null(predictor.binary)){

      if(length(unique(gold))!=2){
        message("Error! The 'gold' variable must be binomial!");error=TRUE
      }
      if(!("numeric" %in% is(predictor))){
        message("Error! The 'predictor' variable must be numeric!");error=TRUE
      }

      if(!error){
        if(positive=="l"){roc.rst=roc(gold~predictor,auc=TRUE,ci=TRUE,direction = "<")}
        if(positive=="s"){roc.rst=roc(gold~predictor,auc=TRUE,ci=TRUE,direction = ">")}
        sensitivities=roc.rst$sensitivities
        specificities=roc.rst$specificities
        cutpoints=roc.rst$thresholds
        yuedeng=sensitivities+specificities-1
        cut.s=cutpoints[yuedeng==max(yuedeng)]
        se.s=sensitivities[yuedeng==max(yuedeng)]
        sp.s=specificities[yuedeng==max(yuedeng)]
        if(important=="se"){
          cut=cut.s[se.s==max(se.s)]
          se=se.s[se.s==max(se.s)]
          sp=sp.s[se.s==max(se.s)]
        }
        if(important=="sp"){
          cut=cut.s[sp.s==max(sp.s)]
          se=se.s[sp.s==max(sp.s)]
          sp=sp.s[sp.s==max(sp.s)]
        }

        AUC=roc.rst$auc
        AUC.SE=(roc.rst$ci[3]-roc.rst$ci[2])/1.96
        AUC.low=roc.rst$ci[1]
        AUC.up=roc.rst$ci[3]
        wilc.t=wilcox.test(predictor[table.gold == 1], predictor[table.gold == 0], alternative = "great", exact=exact)
        P=wilc.t$p.value
        if(AUC.low>0.5 & P>=0.05){
          wilc.t=wilcox.test(predictor[table.gold == 1], predictor[table.gold == 0], alternative = "less", exact=exact)
          P=wilc.t$p.value
        }
        if(AUC.low<=0.5 & P<0.05){
          wilc.t=wilcox.test(predictor[table.gold == 1], predictor[table.gold == 0], alternative = "less", exact=exact)
          P=wilc.t$p.value
        }

        predictor.binary=rep(0,length(roc.rst$predictor))
        predictor.binary[roc.rst$predictor>=cut]=1
        predictor.binary=as.factor(predictor.binary)
        levels(predictor.binary)=c("0","1")
        if(positive=='l'){predictor.binary=factor(predictor.binary,levels=c(1,0))}
        if(positive=='s'){predictor.binary=factor(predictor.binary,levels=c(0,1))}
        pre.tbl=table(predictor.binary,table.gold)
        pre=as.vector(pre.tbl)
        a=tp=pre[1]#a
        b=fp=pre[3]#b
        c=fn=pre[2]#c
        d=tn=pre[4]#d

        acc=(tp+tn)/sum(pre)
        acc.low=acc-1.96*(acc*(1-acc)/sum(pre))
        acc.up=acc+1.96*(acc*(1-acc)/sum(pre))

        se.low=se-1.96*sqrt(se*(1-se)/(tp+fn))
        se.up=se+1.96*sqrt(se*(1-se)/(tp+fn))

        sp.low=sp-1.96*sqrt(sp*(1-sp)/(fp+tn))
        sp.up=sp+1.96*sqrt(sp*(1-sp)/(fp+tn))

        PLR=se/(1-sp)
        PLR.low=PLR*exp(-1.96*sqrt((1-se)/tp+sp/fp))
        PLR.up=PLR*exp(1.96*sqrt((1-se)/tp+sp/fp))

        NLR=(1-se)/sp
        NLR.low=NLR*exp(-1.96*sqrt(se/fn+(1-sp)/tn))
        NLR.up=NLR*exp(1.96*sqrt(se/fn+(1-sp)/tn))

        PPV=tp/(tp+fp)
        PPV.low=PPV-1.96*sqrt(PPV*(1-PPV)/(tp+fp))
        PPV.up=PPV+1.96*sqrt(PPV*(1-PPV)/(tp+fp))

        NPV=tn/(fn+tn)
        NPV.low=NPV-1.96*sqrt(NPV*(1-NPV)/(tn+fn))
        NPV.up=NPV+1.96*sqrt(NPV*(1-NPV)/(tn+fn))

        PPA=a/(a+c)
        PPA.low=PPA-1.96*sqrt(PPA*(1-PPA)/(a+c))
        PPA.up=PPA+1.96*sqrt(PPA*(1-PPA)/(a+c))

        NPA=d/(b+d)
        NPA.low=NPA-1.96*sqrt(NPA*(1-NPA)/(b+d))
        NPA.up=NPA+1.96*sqrt(NPA*(1-NPA)/(b+d))

        TPA=(a+d)/sum(pre)
        TPA.low=TPA-1.96*sqrt(TPA*(1-TPA)/sum(pre))
        TPA.up=TPA+1.96*sqrt(TPA*(1-TPA)/sum(pre))

        KAPPA.rst=Kappa(pre.tbl)
        KAPPA=KAPPA.rst$Unweighted["value"]
        KAPPA.low=confint(KAPPA.rst)[1,1]
        KAPPA.up=confint(KAPPA.rst)[1,2]

        if(AUC.up>1){AUC.up=1}
        if(acc.up>1){acc.up=1}
        if(se.up>1){se.up=1}
        if(sp.up>1){sp.up=1}

        rst=round(data.frame(
          Cutoff=cut,
          AUC,AUC.SE,AUC.low,AUC.up,P,
          ACC=acc,ACC.low=acc.low,ACC.up=acc.up,
          SEN=se,SEN.low=se.low,SEN.up=se.up,
          SPE=sp,SPE.low=sp.low,SPE.up=sp.up,
          PLR,PLR.low,PLR.up,
          NLR,NLR.low,NLR.up,
          PPV,PPV.low,PPV.up,
          NPV,NPV.low,NPV.up,
          PPA,PPA.low,PPA.up,
          NPA,NPA.low,NPA.up,
          TPA,TPA.low,TPA.up,
          KAPPA,KAPPA.low,KAPPA.up),3)
        rst[1,]=sprintf("%.3f",rst[1,])
        row.names(rst)=c("")

        if(plot){
          par(mai=c(1,1,0.3,0.3))
          plot(1-specificities,sensitivities,col="black",xlab=xlab,
               ylab=ylab,type='l',lwd=1,lty=1,cex.lab=1,cex.axis=1)
          points(1-sp,se,col="grey",cex=1,pch=16)
          text(1-sp,se,pos=4,paste("(",round(1-sp,2),", ",round(se,2),")",sep=""))
          legend("bottomright",paste0("AUC = ",sprintf("%.2f",AUC),"(",sprintf("%.2f",AUC.low),"-",sprintf("%.2f",AUC.up),")"),col=1,bty="n",cex=1)
          lines(c(0,1),c(0,1),lty=2)
        }
        return(rst)
      }
    }

    ### if the predictor is binomial

    if(is.null(predictor) & !is.null(predictor.binary)){

      predictor.binary=as.factor(predictor.binary)
      levels(predictor.binary)=c("0","1")
      if(length(unique(predictor.binary))!=2){
        message("Error! The 'predictor.binay' variable must be binomial!");error=TRUE
      }

      if(!error){

        if(positive=='l'){
          predictor.binary=factor(predictor.binary,levels=c(1,0))
        }
        if(positive=='s'){
          predictor.binary=factor(predictor.binary,levels=c(0,1))
        }
        pre.tbl=table(predictor.binary,table.gold)
        pre=as.vector(pre.tbl)
        a=tp=pre[1]
        b=fp=pre[3]
        c=fn=pre[2]
        d=tn=pre[4]

        acc=(tp+tn)/sum(pre)
        acc.low=acc-1.96*(acc*(1-acc)/sum(pre))
        acc.up=acc+1.96*(acc*(1-acc)/sum(pre))

        se=tp/(tp+fn)
        sp=tn/(fp+tn)

        se.low=se-1.96*sqrt(se*(1-se)/(tp+fn))
        se.up=se+1.96*sqrt(se*(1-se)/(tp+fn))

        sp.low=sp-1.96*sqrt(sp*(1-sp)/(fp+tn))
        sp.up=sp+1.96*sqrt(sp*(1-sp)/(fp+tn))

        PLR=se/(1-sp)
        PLR.low=PLR*exp(-1.96*sqrt((1-se)/tp+sp/fp))
        PLR.up=PLR*exp(1.96*sqrt((1-se)/tp+sp/fp))

        NLR=(1-se)/sp
        NLR.low=NLR*exp(-1.96*sqrt(se/fn+(1-sp)/tn))
        NLR.up=NLR*exp(1.96*sqrt(se/fn+(1-sp)/tn))

        PPV=tp/(tp+fp)
        PPV.low=PPV-1.96*sqrt(PPV*(1-PPV)/(tp+fp))
        PPV.up=PPV+1.96*sqrt(PPV*(1-PPV)/(tp+fp))

        NPV=tn/(fn+tn)
        NPV.low=NPV-1.96*sqrt(NPV*(1-NPV)/(tn+fn))
        NPV.up=NPV+1.96*sqrt(NPV*(1-NPV)/(tn+fn))

        PPA=a/(a+c)
        PPA.low=PPA-1.96*sqrt(PPA*(1-PPA)/(a+c))
        PPA.up=PPA+1.96*sqrt(PPA*(1-PPA)/(a+c))

        NPA=d/(b+d)
        NPA.low=NPA-1.96*sqrt(NPA*(1-NPA)/(b+d))
        NPA.up=NPA+1.96*sqrt(NPA*(1-NPA)/(b+d))

        TPA=(a+d)/sum(pre)
        TPA.low=TPA-1.96*sqrt(TPA*(1-TPA)/sum(pre))
        TPA.up=TPA+1.96*sqrt(TPA*(1-TPA)/sum(pre))

        KAPPA.rst=Kappa(pre.tbl)
        KAPPA=KAPPA.rst$Unweighted["value"]
        KAPPA.low=confint(KAPPA.rst)[1,1]
        KAPPA.up=confint(KAPPA.rst)[1,2]

        AUC=se*(1-sp)/2+(se+1)*sp/2
        AUC.low=se.low*(1-sp.low)/2+(se.low+1)*sp.low/2
        AUC.up=se.up*(1-sp.up)/2+(se.up+1)*sp.up/2
        AUC.SE=(AUC.up-AUC.low)/(2*1.96)
        AUC.up=ifelse(AUC.up>1,1,AUC.up)
        wilc.t=wilcox.test(as.numeric(predictor.binary[table.gold == 1]),
                           as.numeric(predictor.binary[table.gold == 0]), alternative = "great", exact=exact)
        P=wilc.t$p.value
        if(AUC.low>0.5 & P>=0.05){
          wilc.t=wilcox.test(as.numeric(predictor.binary[table.gold == 1]),
                             as.numeric(predictor.binary[table.gold == 0]), alternative = "less", exact=exact)
          P=wilc.t$p.value
        }
        if(AUC.low<=0.5 & P<0.05){
          wilc.t=wilcox.test(as.numeric(predictor.binary[table.gold == 1]),
                             as.numeric(predictor.binary[table.gold == 0]), alternative = "less", exact=exact)
          P=wilc.t$p.value
        }

        if(AUC.up>1){AUC.up=1}
        if(acc.up>1){acc.up=1}
        if(se.up>1){se.up=1}
        if(sp.up>1){sp.up=1}

        rst=round(data.frame(
          AUC,AUC.SE,AUC.low,AUC.up,P,
          ACC=acc,ACC.low=acc.low,ACC.up=acc.up,
          SEN=se,SEN.low=se.low,SEN.up=se.up,
          SPE=sp,SPE.low=sp.low,SPE.up=sp.up,
          PLR,PLR.low,PLR.up,
          NLR,NLR.low,NLR.up,
          PPV,PPV.low,PPV.up,
          NPV,NPV.low,NPV.up,
          PPA,PPA.low,PPA.up,
          NPA,NPA.low,NPA.up,
          TPA,TPA.low,TPA.up,
          KAPPA,KAPPA.low,KAPPA.up),3)
        rst[1,]=sprintf("%.3f",rst[1,])
        row.names(rst)=c("")

        if(plot){
          par(mai=c(1,1,0.3,0.3))
          plot(c(0,1-sp,1),c(0,se,1),col="black",xlab=xlab,
               ylab=ylab,type='l',lwd=1,lty=1,cex.lab=1,cex.axis=1)
          points(1-sp,se,col="grey",cex=1,pch=16)
          text(1-sp,se,pos=4,paste("(",round(1-sp,2),", ",round(se,2),")",sep=""))
          legend("bottomright",paste0("AUC = ",sprintf("%.2f",AUC),"(",sprintf("%.2f",AUC.low),"-",sprintf("%.2f",AUC.up),")"),col=1,bty="n",cex=1)
          lines(c(0,1),c(0,1),lty=2)
        }

        return(rst)
      }
    }
  }

}
