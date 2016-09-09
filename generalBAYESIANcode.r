splt=function(x,split,num) {

  x<-as.character(x)
  p1=strsplit(x,split=split,fixed=T)
  p2=unlist(lapply(p1,"[",num))
  return(p2)
}

draws=function(n.iter,n.burnin,n.thin,n.chains=3) {
  return(((n.iter-n.burnin)/n.thin)*n.chains)
}


makeTransparent<-function(someColor, alpha=100)
{
  newColor<-someColor + alpha #I wish
  return(newColor)
}


rmse=function(m,x)
{
  sqrt(mean((x-m)^2))
}

mae=function(m,x)
{
  mean(abs(((x-m))))
}

cloglog=function (theta, bvalue = NULL, inverse = FALSE, deriv = 0, short = TRUE, 
          tag = FALSE) 
{
  if (is.character(theta)) {
    string <- if (short) 
      paste("cloglog(", theta, ")", sep = "")
    else paste("log(-log(1-", theta, "))", sep = "")
    if (tag) 
      string <- paste("Complementary log-log:", string)
    return(string)
  }
  if (!inverse && length(bvalue)) {
    theta[theta <= 0] <- bvalue
    theta[theta >= 1] <- 1 - bvalue
  }
  if (inverse) {
    if (deriv > 0) {
      1/Recall(theta = theta, bvalue = bvalue, inverse = FALSE, 
               deriv = deriv)
    }
    else {
      junk <- exp(theta)
      -expm1(-junk)
    }
  }
  else {
    switch(deriv + 1, {
      log(-log1p(-theta))
    }, -(1 - theta) * log1p(-theta), {
      junk <- log1p(-theta)
      -(1 - theta) * (1 + junk) * junk
    }, stop("argument 'deriv' unmatched"))
  }
}


binplot<-
function (size,surv, ncuts = 70,...,addy="NO",lcol="darkgreen")
{
    os <- order(size)
    os.surv <- (surv)[os]
    os.size <- (size)[os]
    psz <- tapply(os.size, as.numeric(cut(os.size, ncuts)), mean,
        na.rm = TRUE)
    ps <- tapply(os.surv, as.numeric(cut(os.size, ncuts)), mean,
        na.rm = TRUE)
    
    if(addy=="NO") {
        plot(as.numeric(psz), as.numeric(ps),...) }
    
    else{ points(as.numeric(psz), as.numeric(ps),...) }
        lines(as.numeric(psz), as.numeric(ps),...,col=lcol,lwd=2)
        return(cbind(as.numeric(psz),as.numeric(ps)))}


#convert between population and individual data
expand.dft <- function(x, na.strings = "NA", as.is = FALSE, dec = ".") {
    # Take each row in the source data frame table and replicate it
    # using the Freq value
    DF <- sapply(1:nrow(x), 
                 function(i) x[rep(i, each = x$Freq[i]), ],
                 simplify = FALSE)

    # Take the above list and rbind it to create a single DF
    # Also subset the result to eliminate the Freq column
    DF <- subset(do.call("rbind", DF), select = -Freq)

    # Now apply type.convert to the character coerced factor columns  
    # to facilitate data type selection for each column 
    for (i in 1:ncol(DF)) {
        DF[[i]] <- type.convert(as.character(DF[[i]]),
                                na.strings = na.strings,
                                as.is = as.is, dec = dec)
    }

    DF
}



#NOTE THAT FOR GLMS YOU CAN USE PREDICT INSTEAD OF X
  lz<-function(x) {log(x+1)}
  
  RSS1<-function(m,x) {sum((lz((m))-lz(x))^2)}
  TSS1<-function(x) {sum((lz(x)-lz(mean(x)))^2)}
  
  
  RSS<-function(m,x) {sum((((m))-(x))^2)}
  TSS<-function(x) {sum(((x)-(mean(x)))^2)}
  
  
  R2<-function(m,x,logo=F) {
    mors<-na.omit(data.frame(m,x))
    m<-mors$m
    x<-mors$x
    
    if(logo==F) {return(1-RSS(m,x)/TSS(x))} else {return(1-RSS1(m,x)/TSS1(x))}
  }


#summary[[2]]
gisum<-function(model,Start=1,Finish=20) {return(model[[2]]$summary[Start:Finish,])}

gipar<-function(model) {return(model[2]$BUGSoutput$sims.list)}

Rhat<-function(model) {
return(model[[2]]$summary[,8])
}

Q<-function(model,variable,colo="light gray") {
simlist<-gipar(model)

simvar<-which(names(simlist)==variable)

plot(density(simlist[[simvar]]),main=variable)
polygon(density(simlist[[simvar]]),col=colo)

abline(v=quantile(simlist[[simvar]],c(0.025,0.5,0.975)),lwd=2)

return(quantile(simlist[[simvar]],c(0.025,0.5,0.975)))}


pplotCOEF<-function(sims,thickness=2,lab1="",lab2="Values",listorno="T",...,cexX=1.1) {      #this is designed to plot the 95% confidence intervals
if(listorno=="T") {simlist<-sims

simlist$deviance<-NULL

} else{
simlist<-gipar(sims)
dumb<-which(names(simlist)=="deviance")

simlist<-simlist[-dumb]
}
variablenames<-names(simlist)
xo<-c(1:length(names(simlist)))

upo<-max(unlist(lapply(simlist,ups)))
downo<-min(unlist(lapply(simlist,lows)))

plot(0,xo[1],pch="",ylim=c(downo,upo),xlim=c(1,max(xo)),xlab=lab1,ylab=lab2,xaxt="n",...)
 axis(side=1,at=c(1:max(xo)),labels=variablenames,cex.axis=cexX,las=2)
grid()
 abline(h=0,lty=2,col="gray",lwd=thickness)
for(i in 1:length(xo)) {
valz<-simlist[[i]]
m<-median(valz)
u<-ups(valz)
l<-lows(valz)

u5<-ups5(valz)
l5<-lows5(valz)


segments(x0=xo[i],y0=l5,y1=u5,lty=1,col="gray",lwd=10)
points(m~xo[i],pch=19,col="black",cex=1.5)
segments(x0=xo[i],y0=l,y1=u,lty=1,col="black",lwd=3)

}
}


pplot2<-function(input,pvalz,lab1="Predictor",cexo=0.5,lab2="Response",ylimo=c(min(l),max(u)),
xlimo=c(min(input),max(input)),maino="",meormed="med") {      #this is designed to plot the 95% confidence intervals

if(meormed=="me") {
m<-apply(pvalz,2,mean)
}
else{m<-apply(pvalz,2,median)}

u<-apply(pvalz,2,ups)
l<-apply(pvalz,2,lows)
have<-data.frame(l,m,u,input)
colnames(have)<-c("lower","mean","upper","input")
have<-have[order(have$input),]
m<-have$m
l<-have$l
u<-have$u
ip<-have$input

plot(m~ip,ylim=ylimo,xlim=xlimo,pch=19,cex=cexo,xlab=lab1,ylab=lab2,col="black",main=maino,cex.lab=1.8,cex.axis=1.8)
#lines(m~ip,lty=2)
for(i in 1:length(m)) {
segments(x0=ip[i],y0=l[i],y1=u[i],lty=2,col="gray")
#points(m~ip,pch=19,cex=cexo,col="black")
}
}

ups<-function(listie) {
want<-length(listie)*0.975
take<-sort(listie)[want]
return(take)
}

ups5<-function(listie) {
  want<-length(listie)*0.75
  take<-sort(listie)[want]
  return(take)
}



lows5<-function(listie) {
  want<-length(listie)*0.25
  take<-sort(listie)[want]
  return(take)
}


lows<-function(listie) {
want<-length(listie)*0.025
take<-sort(listie)[want]
return(take)
}

ppointFX<-function(input,predicted) { #this is designed for plotting the 95% confidence intervals around a functional relationship
pvalz<-predicted
m<-apply(pvalz,2,median)
u<-apply(pvalz,2,ups)
l<-apply(pvalz,2,lows)
have<-data.frame(l,m,u,input)
colnames(have)<-c("lower","mean","upper","input")
have<-have[order(have$input),]
m<-have$m
l<-have$l
u<-have$u
ip<-have$input
points(m~ip,ylim=c(min(l),max(u)),pch=19,col="gray")
lines(m~ip,lty=2)
#points(u~ip,pch=19,col="gray") points is to confusing
lines(u~ip,lty=2,col="red") #the lines might be good if
#points(l~ip,pch=19,col="gray")
lines(l~ip,lty=2,col="red")
return(have)
}



pplotFX<-function(input,predicted,lab1,lab2) { #this is designed for plotting the 95% confidence intervals around a functional relationship
pvalz<-predicted
m<-apply(pvalz,2,median)
u<-apply(pvalz,2,ups)
l<-apply(pvalz,2,lows)
have<-data.frame(l,m,u,input)
colnames(have)<-c("lower","mean","upper","input")
have<-have[order(have$input),]
m<-have$m
l<-have$l
u<-have$u
ip<-have$input
plot(m~ip,ylim=c(min(l),max(u)),pch=19,col="gray",xlab=lab1,ylab=lab2)
lines(m~ip,lty=2)
#points(u~ip,pch=19,col="gray") points is to confusing
lines(u~ip,lty=2,col="red") #the lines might be good if
#points(l~ip,pch=19,col="gray")
lines(l~ip,lty=2,col="red")
return(have)
}



pplot<-function(input,predicted,lab1,lab2) {      #this is designed to plot the 95% confidence intervals
pvalz<-predicted
m<-apply(pvalz,2,median)
u<-apply(pvalz,2,ups)
l<-apply(pvalz,2,lows)
have<-data.frame(l,m,u,input)
colnames(have)<-c("lower","mean","upper","input")
have<-have[order(have$input),]
m<-have$m
l<-have$l
u<-have$u
ip<-have$input
plot(m~ip,ylim=c(min(l),max(u)),pch=19,col="gray",xlab=lab1,ylab=lab2)
#lines(m~ip,lty=2)
for(i in 1:length(m)) {
segments(x0=ip[i],y0=l[i],y1=u[i],lty=2,col="gray70")
}
#points(u~ip,pch=19,col="gray") points is to confusing
#lines(u~ip,lty=2,col="red") #the lines might be good if
#points(l~ip,pch=19,col="gray")
#lines(l~ip,lty=2,col="red")
return(have)
}



ppoint<-function(input,predicted,lab1,lab2) {      #this is designed to plot the 95% confidence intervals
pvalz<-predicted
m<-apply(pvalz,2,median)
u<-apply(pvalz,2,ups)
l<-apply(pvalz,2,lows)
have<-data.frame(l,m,u,input)
colnames(have)<-c("lower","mean","upper","input")
have<-have[order(have$input),]
m<-have$m
l<-have$l
u<-have$u
ip<-have$input
points(m~ip,ylim=c(min(l),max(u)),pch=19,col="gray")
#lines(m~ip,lty=2)
for(i in 1:length(m)) {
segments(x0=ip[i],y0=l[i],y1=u[i],lty=2,col="gray70")
}
#points(u~ip,pch=19,col="gray") points is to confusing
#lines(u~ip,lty=2,col="red") #the lines might be good if
#points(l~ip,pch=19,col="gray")
#lines(l~ip,lty=2,col="red")
return(have)
}


#I'm pretty sure this is wrong-o, wrong-o!
binMSPE<-function(y,sims.p) {
p<-apply(sims.p,2,mean)
mspe<-sum((p-y)^2)+sum(p*(1-p))
m<-sum((p-y)^2)
great<-data.frame(mspe,m)
names(great)<-c("MSPE","mean only")
return(great)
}



#I'm pretty sure this is wrong-o, wrong-o!
binMSPE<-function(y,n=1,sims.p) {
  p<-apply(sims.p,2,mean)
  mspe<-sum((n*p-n*y)^2)+sum((n*p)*(1-p))
  m<-sum((p-y)^2)
  great<-data.frame(mspe,m)
  names(great)<-c("MSPE","mean only")
  return(great)
}


#(2/(length(y)*(length(y)+1))) #an alternative. I'm not sure which of these to use

poisMSPE<-function(y,sims.p,its=length(sims.p[1,])) {


runny<-rep(0,times=its)
for(i in 1:its) {

runny[i]<-(1/(length(y)))*sum((log(sims.p[,i]+0.5)-log(y+0.5))^2)
}
return(runny)
}


poisMSPE1it<-function(y,sims.p) {


runny<-(1/(length(y)))*sum((log(sims.p+0.5)-log(y+0.5))^2)

return(runny)
}


normalMSPE<-function(y,sims.p) {
p<-apply(sims.p,2,mean)
v<-apply(sims.p,2,var)
mspe<-sum((p-y)^2)+sum(v)
m<-sum((p)^2)
great<-data.frame(mspe,m)
names(great)<-c("MSPE","mean only")
return(great)
}


#OR IF SIGMA IS VARIABLE:

normalMSPEvar<-function(y,sims.p,sigma) {
p<-apply(sims.p,2,mean)
sigma.m<-apply(sigma,2,mean)
mspe<-sum((p-y)^2)+sum(sigma.m)
m<-sum((p)^2)
great<-data.frame(mspe,m)
names(great)<-c("MSPE","mean only")
return(great)
}


#what's next:
#1. experiment with value for total seeds
#2. put in covariates for final model! yay!
#3.. dealing with missing values

###NOTES:
#the lmer model for site with a binomial model is 0 for all random effects...seems problematic

credint<-function(x) {
gret<-sort(x)
up<-round(length(x)*0.95)
down<-round(length(x)*0.05)
final<-c(gret[down],median(gret),mean(gret),gret[up])
names(final)<-c("lower 5%","median","mean","upper 95%")
return(final)
}

stdize<-function(x) {(x-mean(x))/(2*sd(x))}


transdis<-function(x,y) {sqrt((521995-x)^2+(1727781-y)^2)
}




facTOlevs<-function(fac) {
levs<-levels(fac)
mats<-matrix(NA,nrow=length(fac),ncol=length(levs))
for(i in 1:length(levs)) {
mats[,i]<-ifelse(fac==levs[i],1,0)
}
colnames(mats)<-levs
return(data.frame(mats))
}


recodeFACTOR<-function(FACTOR) {
number<-as.numeric(as.factor(as.numeric(FACTOR)))
should<-c(1:max(number))
if(length(which(should %in% number==F))>0) {
return("BROKE BROKE BROKE")
}
else {return(number)}
}

Rhatty<-function(listo) {
  rhat<-Rhat(listo)
  hist(rhat,breaks=100)
  print(max(rhat))
  return(listo[[2]]$summary[which(listo[[2]]$summary[,8]>1.1),])
}