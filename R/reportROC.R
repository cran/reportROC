reportROC=function(gold,
                   predictor=NULL,
                   predictor.binary=NULL,
                   important="se",
                   plot=TRUE,positive='l'){

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
    table.gold=2-as.numeric(as.factor(gold))#'0' is case

    if(!is.null(predictor) & is.null(predictor.binary)){

      if(length(unique(gold))!=2){
        message("Error! The 'gold' variable must be binomial!");error=TRUE
      }
      if(class(predictor)!="numeric"){
        message("Error! The 'predictor' variable must be numeric!");error=TRUE
      }

      if(!error){
        if(positive=="l"){
          roc.rst=roc(gold~predictor,auc=TRUE,ci=TRUE,direction = "<")
        }
        if(positive=="s"){
          roc.rst=roc(gold~predictor,auc=TRUE,ci=TRUE,direction = ">")
        }
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
        AUC.low=roc.rst$ci[2]
        AUC.up=roc.rst$ci[3]

        predictor.binary=rep(0,length(roc.rst$predictor))
        predictor.binary[roc.rst$predictor>=cut]=1
        predictor.binary=as.factor(predictor.binary)
        levels(predictor.binary)=c("0","1")
        if(positive=='l'){
          predictor.binary=factor(predictor.binary,levels=c(1,0))
        }
        if(positive=='s'){
          predictor.binary=factor(predictor.binary,levels=c(0,1))
        }
        pre=table(predictor.binary,table.gold)
        pre=as.vector(pre)
        tp=pre[1]
        fp=pre[3]
        fn=pre[2]
        tn=pre[4]

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

        PPV=pre[1]/(pre[1]+pre[3])
        PPV.low=PPV-1.96*sqrt(PPV*(1-PPV)/(tp+fp))
        PPV.up=PPV+1.96*sqrt(PPV*(1-PPV)/(tp+fp))

        NPV=pre[4]/(pre[4]+pre[2])
        NPV.low=NPV-1.96*sqrt(NPV*(1-NPV)/(tn+fn))
        NPV.up=NPV+1.96*sqrt(NPV*(1-NPV)/(tn+fn))

        rst=round(data.frame(
          Cutoff=cut,
          AUC,AUC.SE,AUC.low,AUC.up,
          ACC=acc,ACC.low=acc.low,ACC.up=acc.up,
          SEN=se,SEN.low=se.low,SEN.up=se.up,
          SPE=sp,SPE.low=sp.low,SPE.up=sp.up,
          PLR,PLR.low,PLR.up,
          NLR,NLR.low,NLR.up,
          PPV,PPV.low,PPV.up,
          NPV,NPV.low,NPV.up),3)

        if(plot){
          par(mai=c(1,1,0.3,0.3))
          plot(1-specificities,sensitivities,col="black",xlab="1-Specificity",
               ylab="Sensetivity",type='l',lwd=1,lty=1,cex.lab=1,cex.axis=1)
          points(1-sp,se,col="grey",cex=1,pch=16)
          text(1-sp,se,pos=4,paste("(",round(1-sp,2),", ",round(se,2),")",sep=""))
          legend("bottomright",paste("AUC =",round(AUC,2)),col=1,bty="n",cex=1)
        }
        return(rst)
      }
    }

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
        pre=table(predictor.binary,table.gold)
        pre=as.vector(pre)
        tp=pre[1]
        fp=pre[3]
        fn=pre[2]
        tn=pre[4]

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

        PPV=pre[1]/(pre[1]+pre[3])
        PPV.low=PPV-1.96*sqrt(PPV*(1-PPV)/(tp+fp))
        PPV.up=PPV+1.96*sqrt(PPV*(1-PPV)/(tp+fp))

        NPV=pre[4]/(pre[4]+pre[2])
        NPV.low=NPV-1.96*sqrt(NPV*(1-NPV)/(tn+fn))
        NPV.up=NPV+1.96*sqrt(NPV*(1-NPV)/(tn+fn))

        AUC=se*(1-sp)/2+(se+1)*sp/2
        AUC.low=se.low*(1-sp.low)/2+(se.low+1)*sp.low/2
        AUC.up=se.up*(1-sp.up)/2+(se.up+1)*sp.up/2
        AUC.SE=(AUC.up-AUC.low)/(2*1.96)
        AUC.up=ifelse(AUC.up>1,1,AUC.up)

        rst=round(data.frame(
          AUC,AUC.SE,AUC.low,AUC.up,
          ACC=acc,ACC.low=acc.low,ACC.up=acc.up,
          SEN=se,SEN.low=se.low,SEN.up=se.up,
          SPE=sp,SPE.low=sp.low,SPE.up=sp.up,
          PLR,PLR.low,PLR.up,
          NLR,NLR.low,NLR.up,
          PPV,PPV.low,PPV.up,
          NPV,NPV.low,NPV.up),3)

        if(plot){
          par(mai=c(1,1,0.3,0.3))
          plot(c(0,1-sp,1),c(0,se,1),col="black",xlab="1-Specificity",
               ylab="Sensetivity",type='l',lwd=1,lty=1,cex.lab=1,cex.axis=1)
          points(1-sp,se,col="grey",cex=1,pch=16)
          text(1-sp,se,pos=4,paste("(",round(1-sp,2),", ",round(se,2),")",sep=""))
          legend("bottomright",paste("AUC =",round(AUC,2)),col=1,bty="n",cex=1)
        }

        return(rst)
      }
    }
  }

}
