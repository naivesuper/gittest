## This script takes ol type as command line arguments and produces the output from sankoff algorithm (including ancestor states and newick tree based on state change) and can be called using $ Rscript sankoff.R oltype
args <- commandArgs(trailingOnly = TRUE)
library("scales")
library("matrixStats")
library("plotrix")
source("~/research/phylogeny/bin/sankoff.R")

########## output direction
out=paste(args,".Rout",sep="")
zz <- file(out, open = "wt")
sink(zz)
sink(zz, type = "message")


## slice the leaf nodes and make them column vector, the leaf nodes are numbered from 1-12
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
## retrieve the states for internal nodes, the internal nodes are numbered from 13-23
n23=sapply(ances.state,"[",i=1) ## to slice the list into seperate nodes
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

########## discard ambiguous sites: unambiguous sites meaning all the internal nodes INCLUDING ancestor states should have only 1 state
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

total.unambi=length(sure); # print out total number of unambiguous sites:
print("total unambiguous site");total.unambi

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

##--pending when it comes to do the stack plot, since all the unambiguous sites would have no change on the ancestral branch without an outgroup, in order to refer changes on the ancestral branch we have to instead use all sites. in other words:
# 1-to calculate the changes on the other branches, we use unambiguous sites including node 23, if not, there will be ambiguity in other nodes and we cannot calculate changes.
# 2-to calcualte the change on the ancestral branch, we compare the states for nodes 22 and 18, if they are the same there is no change, if they are different then we infer the change as 1 based on parsimony principle.

########### build data matrix for the following analysis. The format is SITE(row) X BRANCH(column)
#be23.18=as.numeric(paste(l.n23.1,l.n18.1,sep="")) ##ancestral state will be dealt seperately
be18.16=as.numeric(paste(l.n18.1,l.n16.1,sep=""))
be18.12=as.numeric(paste(l.n18.1,l.n12.1,sep=""))
#be23.22=as.numeric(paste(l.n23.1,l.n22.1,sep=""))
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
branch.site=cbind(be18.16,be18.12,be22.21,be22.9,be21.15,be21.20,be20.19,be20.6,be19.14,be19.17,be17.13,be17.1,be13.2,be13.3,be14.4,be14.5,be15.7,be15.8,be16.10,be16.11)
## the ancestor species is not included in the data matrix right now.

################ States composition in the common ancestor :
l.n23.1=as.vector(unlist(n23[sure])) ## tmp use of umbiguous states only.
print("ancestor with state 1")
ances.1=length(l.n23.1[which(l.n23.1==1)]); ances.1
print("ancestor with state 2")
ances.2=length(l.n23.1[which(l.n23.1==2)]);ances.2
print("ancestor with state 3")
ances.3=length(l.n23.1[which(l.n23.1==3)]);ances.3
print("ancestor with state 4")
ances.4=length(l.n23.1[which(l.n23.1==4)]);ances.4 
print("ancestor with state 5")
ances.5=length(l.n23.1[which(l.n23.1==5)]);ances.5


############################ Ancestral states inference ########################################################
## generate leaf states plot
tablenorm=function(x){ ## the input is of the class of table
   tmp=table(as.matrix(x))
  if (is.na(tmp["1"])){tmp["1"]=0}
  if (is.na(tmp["2"])){tmp["2"]=0}
  if (is.na(tmp["3"])){tmp["3"]=0}
  if (is.na(tmp["4"])){tmp["4"]=0}
  if (is.na(tmp["5"])){tmp["5"]=0}
   tmp=as.matrix(tmp)
   return(tmp)
 }

## data structure, leaf is data frame convert to matrix and then calcuate frequency by table and convert to matrix again.
leaf.plot=c()
for (i in 1: 12){
tmp=tablenorm(leaf[i,]) 
leaf.plot=cbind(leaf.plot,as.matrix(tmp))
}

## deal with ancestral node
ances.plot=c()
tmp=sapply(n23,length)
one.state=unlist(n23[which(tmp==1)])
two.state=unlist(n23[which(tmp==2)])
one.state=tablenorm(one.state)
two.state=tablenorm(two.state)/2

if (max(tmp)==3){
three.state=unlist(n23[which(tmp==3)])
three.state=tablenorm(three.state)/3
ances.plot=one.state+two.state+three.state
} else if (max(tmp)==2){
  ances.plot=one.state+two.state
} else if (max(tmp)==4){
three.state=unlist(n23[which(tmp==3)])
three.state=tablenorm(three.state)/3
four.state=unlist(n23[which(tmp==4)])
four.state=tablenorm(four.state)/4
ances.plot=one.state+two.state+three.state+four.state
} else if (max(tmp)==5){
three.state=unlist(n23[which(tmp==3)])
three.state=tablenorm(three.state)/3
four.state=unlist(n23[which(tmp==4)])
four.state=tablenorm(four.state)/4
five.state=unlist(n23[which(tmp==5)])
five.state=tablenorm(five.state)/5
ances.plot=one.state+two.state+three.state+four.state+five.state
}

state.spectrum=cbind(ances.plot,leaf.plot)/dim(leaf)[2]
colnames(state.spectrum)=c("ances", "dmel", "dsim","dsec","dyak","dere","dana","dpse","dper","dwil","dvir","dmoj","dgri")

file=paste(args,"_state_spectrum.svg",sep="")
svg(file,width=6,height=6)
test=barplot(state.spectrum,col=c("blue","red","orange","yellow","green"),xlab="Species",ylab="Percentage",cex.main=0.6,xaxt="n",yaxt="n")
axis(1,at=test,las=1,labels=c("ances", "dmel", "dsim","dsec","dyak","dere","dana","dpse","dper","dwil","dvir","dmoj","dgri"),cex.axis=0.6, tck=0.01)
axis(2,las=1,tck=0.01,cex.lab=0.7, cex.axis=0.6)
dev.off()

###########################Evolutionary events inference##############################################################
################### do row sums / column /color sums for the data matrix (row-pairs[when]; column-branches[where]; color-events[what])
## data matrix sorting:
site.change.raw=apply(branch.site,1,function(x){return(length(which(x%%11!=0)))}) ## calculate row-wise events
branch.site.count=cbind(branch.site,l.n23.1,site.change.raw) ## implement with ancestral branch and the row-wise events
## overlapping pairs sorting first by ancestor states, then by events number
branch.site.order=branch.site.count[order(branch.site.count[,21],branch.site.count[,22]),]



#### column stats calculation and visualization
## calcualte column sums
branch.events=apply(branch.site.order[,1:20],2,function(x){length(which(x%%11!=0))})
## count sites switching from state 2/3/4 to state 1
gain=matrix(as.numeric(branch.site%%10==1 & branch.site%%11!=0 & branch.site != 51),nrow=total.unambi)
## count sites switching from state 5 to state 1
rearrange.in=matrix(as.numeric(branch.site == 51),nrow=total.unambi)
## count sites switching from state 1 to state 2/3/4
loss=matrix(as.numeric(branch.site==12 | branch.site==13 | branch.site==14),nrow=total.unambi)
## count sites switching from state 1 to state 5
rearrange.out=matrix(as.numeric(branch.site == 15),nrow=total.unambi)
## sum events for each branch 
gain.branchsum=colSums(gain)
rearrange.in.branchsum=colSums(rearrange.in)
loss.branchsum=colSums(loss)
rearrange.out.branchsum=colSums(rearrange.out)
## calculate changes on the ancestral branch separately
tmp=sapply(n23,length)
consensus=unlist(n23[which(tmp==1)])
ances.conserve=length(consensus[which(consensus==1)])
ances.nonconserve=length(consensus[which(consensus!=1)])
##count events on ancestral branch ## eg. n18=2, n22=5,
ances.change=length(n23)-length(consensus)
print("ances changes:"); ances.change
## ouput the data table
change.total=c()
for (i in 1:20){
change.total=rbind(change.total,c(gain.branchsum[i],rearrange.in.branchsum[i],loss.branchsum[i],rearrange.out.branchsum[i]))   
}
rownames(change.total)=c("18.16","18.12","22.21","22.9","21.15","21.20","20.19","20.6","19.14","19.17","17.13","17.1","13.2","13.3","14.4","14.5","15.7","15.8","16.10","16.11")
write.table(t(change.total),file=paste(args,"-events.txt",sep=""))
## stacked bar plot
file=paste(args,"-branch-length.svg",sep="")
svg(file,width=7,height=3)
barplot(t(change.total), col=c("red","blue","green","yellow"),horiz=F,cex.axis=0.8,las=1)
dev.off()



