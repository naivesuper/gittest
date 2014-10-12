## This script takes ol type as command line arguments and produces the output from sankoff algorithm (including ancestor states and newick tree based on state change) and can be called using $ Rscript sankoff.R oltype
args <- commandArgs(trailingOnly = TRUE)

library("scales")
## Choose cost matrix ##########################################################
cost=read.table("~/research/phylogeny/data/cost-eq.txt",sep=" ",header=F)


## Input the data ##############################################################
leaf=read.table(file=paste("sankoffin-",args,".txt",sep=""),sep=" ",header=F)

library("matrixStats")
colProds(leaf)

leaf=leaf[,which(colProds(leaf) != 1)]

coln=ncol(leaf)



## Sankoff algorithm inplementation ############################################
ances=function(site){
#initializing terminal nodes  
s1=site[1];s2=site[2];s3=site[3];s4=site[4];s5=site[5];s6=site[6];s7=site[7];s8=site[8];s9=site[9];s10=site[10];s11=site[11];s12=site[12];
#initializing internal nodes
s13=s14=s15=s16=s17=s18=s19=s20=s21=s22=s23=rep(0,6)
## 2nd node
for(i in 1:6) {s13[i]=cost[i,s2]+cost[i,s3]}
for(i in 1:6) {s14[i]=cost[i,s4]+cost[i,s5]}
for(i in 1:6) {s15[i]=cost[i,s7]+cost[i,s8]}
for(i in 1:6) {s16[i]=cost[i,s10]+cost[i,s11]}
S13=matrix(rep(t(s13),6),byrow=F,6,6)
S14=matrix(rep(t(s14),6),byrow=F,6,6)
S15=matrix(rep(t(s15),6),byrow=F,6,6)
S16=matrix(rep(t(s16),6),byrow=F,6,6)
## 3rd node
S17.13=cost+S13
S17.13.min=apply(S17.13,2,function(m){c=min(m);return(c)})
for(i in 1:6) {s17[i]=S17.13.min[i]+cost[i,s1]}
S18.16=cost+S16
S18.16.min=apply(S18.16,2,function(m){c=min(m);return(c)})
for(i in 1:6) {s18[i]=S18.16.min[i]+cost[i,s12]}
S17=matrix(rep(t(s17),6),byrow=F,6,6)
S18=matrix(rep(t(s18),6),byrow=F,6,6)
## 4th node
S19.17=cost+S17
S19.14=cost+S14
S19.17.min=apply(S19.17,2,function(m){c=min(m);return(c)})
S19.14.min=apply(S19.14,2,function(m){c=min(m);return(c)})
s19=S19.17.min+S19.14.min
S19=matrix(rep(t(s19),6),byrow=F,6,6)

S20.19=cost+S19
S20.19.min=apply(S20.19,2,function(m){c=min(m);return(c)})
for(i in 1:6) {s20[i]=S20.19.min[i]+cost[i,s6]}
S20=matrix(rep(t(s20),4),byrow=F,6,6)

S21.20=cost+S20
S21.15=cost+S15
S21.20.min=apply(S21.20,2,function(m){c=min(m);return(c)})
S21.15.min=apply(S21.15,2,function(m){c=min(m);return(c)})
s21=S21.20.min+S21.15.min
S21=matrix(rep(t(s21),6),byrow=F,6,6)

S22.21=cost+S21
S22.21.min=apply(S22.21,2,function(m){c=min(m);return(c)})
for(i in 1:6) {s22[i]=S22.21.min[i]+cost[i,s9]}
S22=matrix(rep(t(s22),6),byrow=F,6,6)

S23.22=cost+S22
S23.18=cost+S18
S23.22.min=apply(S23.22,2,function(m){c=min(m);return(c)})
S23.18.min=apply(S23.18,2,function(m){c=min(m);return(c)})
s23=S23.22.min+S23.18.min

## termination
s23.min=min(s23)
s23.state=which(s23==min(s23))
#print(s23);
#print(s23.min);
#print(s23.state);
                  
## Backtracing
if (length(s23.state)==1){
  slice=S23.22[,s23.state]
  s22.state=which(slice==min(slice))

  slice=S23.18[,s23.state]
  s18.state=which(slice==min(slice))

  if (length(s18.state)==1){ #if 18 is determined consider 16
    slice=S18.16[,s18.state]
    s16.state=which(slice==min(slice))
  }else {s16.state=7
       } #if 18 is ambiguous 16 is ambiguous (encoded as 7)

  if (length(s22.state)==1){
    slice=S22.21[,s22.state]
    s21.state=which(slice==min(slice))
    if (length(s21.state)==1){
      slice=S21.15[,s21.state]
      s15.state=which(slice==min(slice))
      
      slice=S21.20[,s21.state]
      s20.state=which(slice==min(slice))
      if (length(s20.state)==1){
        slice=S20.19[,s20.state]
        s19.state=which(slice==min(slice))
        if (length(s19.state)==1){
          slice=S19.14[,s19.state]
          s14.state=which(slice==min(slice))

          slice=S19.17[,s19.state]
          s17.state=which(slice==min(slice))
          if (length(s17.state)==1){
            slice=S17.13[,s17.state]
            s13.state=which(slice==min(slice))
          }else{s13.state=7
          }
        }else{s17.state=s13.state=s14.state=7
        }
      }else {s19.state=s17.state=s13.state=s14.state=7
      }
    }else{s20.state=s19.state=s17.state=s13.state=s14.state=s15.state=7
    }
  }else{s21.state=s20.state=s19.state=s17.state=s13.state=s14.state=s15.state=7
  } 
}else{s22.state=s21.state=s18.state=s16.state=s21.state=s20.state=s19.state=s17.state=s13.state=s14.state=s15.state=7
}

all.state=list(n23=s23.state,n18=s18.state,n16=s16.state,n22=s22.state,n21=s21.state,n15=s15.state,n20=s20.state,n19=s19.state,n14=s14.state,n17=s17.state,n13=s13.state)
return(all.state)
}

## Calculate ancestral states using the ances function #########################
ances.state=apply(leaf[,1:coln],2,ances)

## Slice the leaf nodes:
n1=leaf[1,];n1=t(n1)
n2=leaf[2,];n2=t(n2)
n3=leaf[3,];n3=t(n3)
n4=leaf[4,];n4=t(n4)
n5=leaf[5,];n5=t(n5)
n6=leaf[6,];n6=t(n6)
n7=leaf[7,];n7=t(n7)
n8=leaf[8,];n8=t(n8)
n9=leaf[9,];n9=t(n9)
n10=leaf[10,];n10=t(n10)
n11=leaf[11,];n11=t(n11)
n12=leaf[12,];n12=t(n12)

## Retrieve the states for internal nodes:
n23=sapply(ances.state,"[",i=1) ##to slice the list into seperate nodes
n18=sapply(ances.state,"[",i=2)
n16=sapply(ances.state,"[",i=3)
n22=sapply(ances.state,"[",i=4)
n21=sapply(ances.state,"[",i=5)
n15=sapply(ances.state,"[",i=6)
n20=sapply(ances.state,"[",i=7)
n19=sapply(ances.state,"[",i=8)
n14=sapply(ances.state,"[",i=9)
n17=sapply(ances.state,"[",i=10)
n13=sapply(ances.state,"[",i=11)

## Find parsimony informative sites
l.n23=l.n13=l.n14=l.n15=l.n16=l.n17=l.n18=l.n19=l.n20=l.n21=l.n22=rep(0,coln)
for (i in 1:coln){
  l.n23[i]=length(unlist(n23[i]))
  l.n13[i]=length(unlist(n13[i]))
  l.n14[i]=length(unlist(n14[i]))
  l.n15[i]=length(unlist(n15[i]))
  l.n16[i]=length(unlist(n16[i]))
  l.n17[i]=length(unlist(n17[i]))
  l.n18[i]=length(unlist(n18[i]))
  l.n19[i]=length(unlist(n19[i]))
  l.n20[i]=length(unlist(n20[i]))
  l.n21[i]=length(unlist(n21[i]))
  l.n22[i]=length(unlist(n22[i]))
}
sure=intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(which(l.n23==1), which(l.n13==1)),which(l.n14==1)),which(l.n15==1)),which(l.n16==1)),which(l.n17==1)),which(l.n18==1)),which(l.n19==1)),which(l.n20==1)),which(l.n21==1)),which(l.n22==1))
l.n23.1=unlist(n23[sure])
l.n13.1=unlist(n13[sure])
l.n14.1=unlist(n14[sure])
l.n15.1=unlist(n15[sure])
l.n16.1=unlist(n16[sure])
l.n17.1=unlist(n17[sure])
l.n18.1=unlist(n18[sure])
l.n19.1=unlist(n19[sure])
l.n20.1=unlist(n20[sure])
l.n21.1=unlist(n21[sure])
l.n22.1=unlist(n22[sure])

l.n1.1=n1[sure]
l.n2.1=n2[sure]
l.n3.1=n3[sure]
l.n4.1=n4[sure]
l.n5.1=n5[sure]
l.n6.1=n6[sure]
l.n7.1=n7[sure]
l.n8.1=n8[sure]
l.n9.1=n9[sure]
l.n10.1=n10[sure]
l.n11.1=n11[sure]
l.n12.1=n12[sure]



## Calculate branch length as the change steps
n23.18=l.n23.1-l.n18.1;#n23.18[which(n23.18!=0)]
n18.16=l.n18.1-l.n16.1;#n18.16[which(n18.16!=0)]
n18.12=l.n18.1-l.n12.1;#n18.12[which(n18.12!=0)]
n23.22=l.n23.1-l.n22.1;#n23.22[which(n23.22!=0)]
n22.21=l.n22.1-l.n21.1;#n22.21[which(n22.21!=0)]
n22.9=l.n22.1-l.n9.1;#n22.9[which(n22.9!=0)]
n21.15=l.n21.1-l.n15.1;#n21.15[which(n21.15!=0)]
n21.20=l.n21.1-l.n20.1;#n21.20[which(n21.20!=0)]
n20.19=l.n20.1-l.n19.1;#n20.19[which(n20.19!=0)]
n20.6=l.n20.1-l.n6.1;#n20.6[which(n20.6!=0)]
n19.14=l.n19.1-l.n14.1;#n19.14[which(n19.14!=0)]
n19.17=l.n19.1-l.n17.1;#n19.17[which(n19.17!=0)]
n17.13=l.n17.1-l.n13.1;#n17.13[which(n17.13!=0)]
n17.1=l.n17.1-l.n1.1;#n17.1[which(n17.1!=0)]
n13.2=l.n13.1-l.n2.1;#n13.2[which(n13.2!=0)]
n13.3=l.n13.1-l.n3.1;#n13.3[which(n13.3!=0)]
n14.4=l.n14.1-l.n4.1;#n14.4[which(n14.4!=0)]
n14.5=l.n14.1-l.n5.1;#n14.5[which(n14.5!=0)]
n15.7=l.n15.1-l.n7.1;#n15.7[which(n15.7!=0)]
n15.8=l.n15.1-l.n8.1;#n15.8[which(n15.8!=0)]
n16.10=l.n16.1-l.n10.1;#n16.10[which(n16.10!=0)]
n16.11=l.n16.1-l.n11.1;#n16.11[which(n16.11!=0)]


## Output the transition matrix ###############################################
x=c(l.n23.1,l.n18.1,l.n18.1,l.n23.1,l.n22.1,l.n22.1,l.n21.1,l.n21.1,l.n20.1,l.n20.1,l.n19.1,l.n19.1,l.n17.1,l.n17.1,l.n13.1,l.n13.1,l.n14.1,l.n14.1,l.n15.1,l.n15.1,l.n16.1,l.n16.1)
y=c(l.n18.1,l.n16.1,l.n12.1,l.n22.1,l.n21.1,l.n9.1,l.n15.1,l.n20.1,l.n19.1,l.n6.1,l.n14.1,l.n17.1,l.n13.1,l.n1.1,l.n2.1,l.n3.1,l.n4.1,l.n5.1,l.n7.1,l.n8.1,l.n10.1,l.n11.1)
z=paste(x,y,sep="")
trans.mat=paste(args[1],".transmat",sep="")
write.table(table(z),file=trans.mat,row.names=FALSE,col.names=FALSE,quote=FALSE)

#parsimony informative site:
total.unambi=length(sure);total.unambi

bl1.nol=0.2155
bl2.nol=0.0756
bl3.nol=0.0387
bl4.nol=0.248
bl5.nol=0.0264
bl6.nol=0.0651
bl7.nol=0.0932
bl8.nol=0.0167
bl9.nol=0.1354
bl10.nol=0.0704
bl11.nol=0.0616
bl12.nol=0.0466
bl13.nol=0.1258
bl14.nol=0.0844
bl15.nol=0.0598
bl16.nol=0.1522
bl17.nol=0
bl18.nol=0.058
bl19.nol=0.0862
bl20.nol=0.0273
bl21.nol=0.1029
bl22.nol=0



##not-normalize by nol overlapping pairs
#bl22=round(length(n23.18[which(n23.18!=0)])/total.unambi/bl22.nol,4) ##bl indicates branch length
bl22=0
bl20=round(length(n18.16[which(n18.16!=0)])/total.unambi,4)
bl21=round(length(n18.12[which(n18.12!=0)])/total.unambi,4)
#bl17=round(length(n23.22[which(n23.22!=0)])/total.unambi/bl17.nol,4)
bl17=0
bl15=round(length(n22.21[which(n22.21!=0)])/total.unambi,4)
bl16=round(length(n22.9[which(n22.9!=0)])/total.unambi,4)
bl14=round(length(n21.15[which(n21.15!=0)])/total.unambi,4)
bl11=round(length(n21.20[which(n21.20!=0)])/total.unambi,4)
bl9=round(length(n20.19[which(n20.19!=0)])/total.unambi,4)
bl10=round(length(n20.6[which(n20.6!=0)])/total.unambi,4)
bl8=round(length(n19.14[which(n19.14!=0)])/total.unambi,4)
bl5=round(length(n19.17[which(n19.17!=0)])/total.unambi,4)
bl3=round(length(n17.13[which(n17.13!=0)])/total.unambi,4)
bl4=round(length(n17.1[which(n17.1!=0)])/total.unambi,4)
bl1=round(length(n13.2[which(n13.2!=0)])/total.unambi,4)
bl2=round(length(n13.3[which(n13.3!=0)])/total.unambi,4)
bl6=round(length(n14.4[which(n14.4!=0)])/total.unambi,4)
bl7=round(length(n14.5[which(n14.5!=0)])/total.unambi,4)
bl12=round(length(n15.7[which(n15.7!=0)])/total.unambi,4)
bl13=round(length(n15.8[which(n15.8!=0)])/total.unambi,4)
bl18=round(length(n16.10[which(n16.10!=0)])/total.unambi,4)
bl19=round(length(n16.11[which(n16.11!=0)])/total.unambi,4)



##normalize by nol overlapping pairs
#bl22=round(length(n23.18[which(n23.18!=0)])/total.unambi/bl22.nol,4) ##bl indicates branch length
#bl22=0
#bl20=round(length(n18.16[which(n18.16!=0)])/total.unambi/bl20.nol,4)
#bl21=round(length(n18.12[which(n18.12!=0)])/total.unambi/bl21.nol,4)
#bl17=round(length(n23.22[which(n23.22!=0)])/total.unambi/bl17.nol,4)
#bl17=0
#bl15=round(length(n22.21[which(n22.21!=0)])/total.unambi/bl15.nol,4)
#bl16=round(length(n22.9[which(n22.9!=0)])/total.unambi/bl16.nol,4)
#bl14=round(length(n21.15[which(n21.15!=0)])/total.unambi/bl14.nol,4)
#bl11=round(length(n21.20[which(n21.20!=0)])/total.unambi/bl11.nol,4)
#bl9=round(length(n20.19[which(n20.19!=0)])/total.unambi/bl9.nol,4)
#bl10=round(length(n20.6[which(n20.6!=0)])/total.unambi/bl10.nol,4)
#bl8=round(length(n19.14[which(n19.14!=0)])/total.unambi/bl8.nol,4)
#bl5=round(length(n19.17[which(n19.17!=0)])/total.unambi/bl5.nol,4)
#bl3=round(length(n17.13[which(n17.13!=0)])/total.unambi/bl3.nol,4)
#bl4=round(length(n17.1[which(n17.1!=0)])/total.unambi/bl4.nol,4)
#bl1=round(length(n13.2[which(n13.2!=0)])/total.unambi/bl1.nol,4)
#bl2=round(length(n13.3[which(n13.3!=0)])/total.unambi/bl2.nol,4)
#bl6=round(length(n14.4[which(n14.4!=0)])/total.unambi/bl6.nol,4)
#bl7=round(length(n14.5[which(n14.5!=0)])/total.unambi/bl7.nol,4)
#bl12=round(length(n15.7[which(n15.7!=0)])/total.unambi/bl12.nol,4)
#bl13=round(length(n15.8[which(n15.8!=0)])/total.unambi/bl13.nol,4)
#bl18=round(length(n16.10[which(n16.10!=0)])/total.unambi/bl18.nol,4)
#bl19=round(length(n16.11[which(n16.11!=0)])/total.unambi/bl19.nol,4)




#ncol(leaf)

#total.len=sum(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9,bl10,bl11,bl12,bl13,bl14,bl15,bl16,bl17,bl18,bl19,bl20,bl21,bl22)
#total.len

## Print the newick tree #######################################################
newick=paste("(((((((dsim:",bl1,",dsec:",bl2,"):",bl3,",dmel:",bl4,"):",bl5,",(dyak:",bl6,",dere:",bl7,"):",bl8,"):",bl9,",dana:",bl10,"):",bl11,",(dpse:",bl12,",dper:",bl13,"):",bl14,"):",bl15,",dwil:",bl16,"):",bl17,",((dvir:",bl18,",dmoj:",bl19,"):",bl20,",dgri:",bl21,"):",bl22,");",sep="")
treeout=paste(args[1],".newick",sep="")
write.table(newick,file=treeout,row.names=FALSE,col.names=FALSE,quote=FALSE)


## Build data matrix for change step site spectrum ##############################
be23.18=as.numeric(paste(l.n23.1,l.n18.1,sep="")) #be indicated branch edge
be18.16=as.numeric(paste(l.n18.1,l.n16.1,sep=""))
be18.12=as.numeric(paste(l.n18.1,l.n12.1,sep=""))
be23.22=as.numeric(paste(l.n23.1,l.n22.1,sep=""))
be22.21=as.numeric(paste(l.n22.1,l.n21.1,sep=""))
be22.9=as.numeric(paste(l.n22.1,l.n9.1,sep=""))
be21.15=as.numeric(paste(l.n21.1,l.n15.1,sep=""))
be21.20=as.numeric(paste(l.n21.1,l.n20.1,sep=""))
be20.19=as.numeric(paste(l.n20.1,l.n19.1,sep=""))
be20.6=as.numeric(paste(l.n20.1,l.n6.1,sep=""))
be19.14=as.numeric(paste(l.n19.1,l.n14.1,sep=""))
be19.17=as.numeric(paste(l.n19.1,l.n17.1,sep=""))
be17.13=as.numeric(paste(l.n17.1,l.n13.1,sep=""))
be17.1=as.numeric(paste(l.n17.1,l.n1.1,sep=""))
be13.2=as.numeric(paste(l.n13.1,l.n2.1,sep=""))
be13.3=as.numeric(paste(l.n13.1,l.n3.1,sep=""))
be14.4=as.numeric(paste(l.n14.1,l.n4.1,sep=""))
be14.5=as.numeric(paste(l.n14.1,l.n5.1,sep=""))
be15.7=as.numeric(paste(l.n15.1,l.n7.1,sep=""))
be15.8=as.numeric(paste(l.n15.1,l.n8.1,sep=""))
be16.10=as.numeric(paste(l.n16.1,l.n10.1,sep=""))
be16.11=as.numeric(paste(l.n16.1,l.n11.1,sep=""))

branch.site=cbind(be23.18,be18.16,be18.12,be23.22,be22.21,be22.9,be21.15,be21.20,be20.19,be20.6,be19.14,be19.17,be17.13,be17.1,be13.2,be13.3,be14.4,be14.5,be15.7,be15.8,be16.10,be16.11)
change.test=matrix(as.numeric(branch.site%%11!=0),nrow=length(be23.18))
change.step=rep(0,length(be23.18))
for (i in 1:length(be23.18)) {
change.step[i]=as.vector(table(change.test[i,])[2])
}
for (i in 1:length(be23.18)) {
if (is.na(change.step[i])){
change.step[i]=0}
}

stepout=paste(args[1],"-change-step.txt",sep="")
write.table(change.step,file=stepout,row.names=FALSE,col.names=FALSE,quote=FALSE)



