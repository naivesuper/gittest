#This script takes overlapping pairs in dmel and look for the orthologs in the additional speicies, calculates the order difference and maps it to state. The format is 12 species (rows) by overlapping pairs (columns).

# state class is
#1-conserved,2-both missing,3-left missing,4-right missing,5-order shuffle
args=commandArgs(trailingOnly = TRUE)

##########setwddata loading
setwd("~/research/phylogeny/results")
input=paste(args,".txt",sep="")
inpath="~/research/phylogeny/data/"
equiva=read.table(file=paste(inpath,"fbgn-equivalence-table.txt",sep=""), sep="\t", head=T)
dmel=read.table(file=paste(inpath,input,sep=""), sep=' ',head=F)

#species indexing
dsim=14;
dsec=15;
dere=16;
dyak=17;
dana=18;
dpse=19;
dper=20;
dwil=21;
dvir=22;
dmoj=23;
dgri=24;

#############ortholog order mapping and difference calculation
sankoffin.order=apply(dmel,
  FUN=function(x) {
states=array(12)
states[1]=1 #this is the dmel ordiff should be 1 for conserved gene pair

for (species in 14:24) {
  
left=equiva[match(as.numeric(x[15]), equiva[,13]),species] #lookup the order from equivalence table
right=equiva[match(as.numeric(x[32]), equiva[,13]),species]


if (left==0 & right ==0){state=2} #order difference calculation and states mapping
if (left ==0 & right !=0){state=3}
if (left !=0 & right ==0){state=4}
if (left !=0 & right !=0){
delta=abs(left-right)
if (delta > 1){state=5}
if (delta==1) {state=1}
}
states=c(states, state)
}
return(states)
}, 1)

output=paste("sankoffin-",args,".txt",sep="")

write.table(sankoffin.order,file=output,row.names = FALSE,
            col.names = FALSE)
