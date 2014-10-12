ances=function(site){
#initializing terminal nodes  
s1=site[1];s2=site[2];s3=site[3];s4=site[4];s5=site[5];s6=site[6];s7=site[7];s8=site[8];s9=site[9];s10=site[10];s11=site[11];s12=site[12];
#initializing internal nodes
s13=s14=s15=s16=s17=s18=s19=s20=s21=s22=s23=rep(0,5)
## 2nd node
for(i in 1:5) {s13[i]=cost[i,s2]+cost[i,s3]}
for(i in 1:5) {s14[i]=cost[i,s4]+cost[i,s5]}
for(i in 1:5) {s15[i]=cost[i,s7]+cost[i,s8]}
for(i in 1:5) {s16[i]=cost[i,s10]+cost[i,s11]}
S13=matrix(rep(t(s13),5),byrow=F,5,5)
S14=matrix(rep(t(s14),5),byrow=F,5,5)
S15=matrix(rep(t(s15),5),byrow=F,5,5)
S16=matrix(rep(t(s16),5),byrow=F,5,5)
## 3rd node
S17.13=cost+S13
S17.13.min=apply(S17.13,2,function(m){c=min(m);return(c)})
for(i in 1:5) {s17[i]=S17.13.min[i]+cost[i,s1]}
S18.16=cost+S16
S18.16.min=apply(S18.16,2,function(m){c=min(m);return(c)})
for(i in 1:5) {s18[i]=S18.16.min[i]+cost[i,s12]}
S17=matrix(rep(t(s17),5),byrow=F,5,5)
S18=matrix(rep(t(s18),5),byrow=F,5,5)
## 4th node
S19.17=cost+S17
S19.14=cost+S14
S19.17.min=apply(S19.17,2,function(m){c=min(m);return(c)})
S19.14.min=apply(S19.14,2,function(m){c=min(m);return(c)})
s19=S19.17.min+S19.14.min
S19=matrix(rep(t(s19),5),byrow=F,5,5)

S20.19=cost+S19
S20.19.min=apply(S20.19,2,function(m){c=min(m);return(c)})
for(i in 1:5) {s20[i]=S20.19.min[i]+cost[i,s6]}
S20=matrix(rep(t(s20),5),byrow=F,5,5)

S21.20=cost+S20
S21.15=cost+S15
S21.20.min=apply(S21.20,2,function(m){c=min(m);return(c)})
S21.15.min=apply(S21.15,2,function(m){c=min(m);return(c)})
s21=S21.20.min+S21.15.min
S21=matrix(rep(t(s21),5),byrow=F,5,5)

S22.21=cost+S21
S22.21.min=apply(S22.21,2,function(m){c=min(m);return(c)})
for(i in 1:5) {s22[i]=S22.21.min[i]+cost[i,s9]}
S22=matrix(rep(t(s22),5),byrow=F,5,5)

S23.22=cost+S22
S23.18=cost+S18
S23.22.min=apply(S23.22,2,function(m){c=min(m);return(c)})
S23.18.min=apply(S23.18,2,function(m){c=min(m);return(c)})
s23=S23.22.min+S23.18.min

## termination
s23.min=min(s23)
s23.state=which(s23==min(s23))

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

