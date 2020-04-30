# Q1) Implement KNN function

myknn<-function(train,test,cl,k,method){
  predicted.class<-c()
  for (i in 1:nrow(test)){
    newdf<-data.frame(Dist=double(),J=integer())
      for (j in 1:nrow(train)){
        sum<-0
        for (k in 1:ncol(train)){
          #standardization
          train.row.mean<-mean(train[,k])
          train.row.sd<-sd(train[,k])
          cur.train<-(train[j,k]-train.row.mean)/train.row.sd
          test.row.mean<-mean(test[,k])
          test.row.sd<-sd(test[,k])
          cur.test<-(train[i,k]-test.row.mean)/test.row.sd
          sum<-sum+(cur.test-cur.train)^2
        }
        sum<-sqrt(sum)
        newdf<-rbind(newdf, data.frame(Dist=sum,J=j))
      }
    newdf<-cbind(newdf,Class=cl)
    newdf_order<-newdf[order(newdf$Dist),]
    newdf_order<-newdf_order[1:k,]
    if(method==1) classify<-names(sort(table(newdf_order$Class),decreasing=TRUE)[1])
    else if(method==2){
      max<-0
      lev_num<-nlevels(newdf_order$Class)
      for(x in levels(newdf_order$Class)){
        subdf<-newdf_order[newdf_order$Class==x,]
        subdf$Dist<-1/subdf$Dist
        sum<-sum(subdf$Dist)
        if(max<sum){
          max<-sum
          classify<-x
        }
      }
    }
    predicted.class<-append(predicted.class,classify)
  }
  predicted.class<-as.factor(predicted.class)
  return(predicted.class)
}
   
# Q2) Implement naive Bayes function (Using Laplace for probability)

myNaiveBayes<-function(train,cl){
  temp<-table(cl)
  sum.temp<-sum(temp)
  for(i in 1:length(temp)){
    temp[i]<-temp[i]/sum.temp
  }
  priors<-temp
  sub<-cbind(Class=cl,train)
  rname<-c()
  cond.prob<-list()
  for(k in levels(cl)) rname<-append(rname,k)
  for (i in 2:ncol(sub)){
    vec<-c()
    cname<-c()
    for(x in levels(sub[,i])) cname<-append(cname,x)
    for (k in 1:length(rname)){
      sub2<-sub[sub[,1]==rname[k],]
      dnom<-nrow(sub2)+length(rname)-sum(is.na(sub2[,i]))
      for(x in 1:length(cname)){
        sub3<-sub2[sub2[,i]==cname[x],]
        numer<-nrow(sub3)+1-sum(is.na(sub3[,i]))
        vec<-append(vec, (numer/dnom))
      }
    }
    mymatrix<-matrix(vec,nlevels(cl),nlevels(sub[,i]),byrow=TRUE,dimnames=list(rname,cname))
    mymatrix_list<-list(mymatrix)
    cond.prob<-append(cond.prob,mymatrix_list)
  }
  return(list(priors,cond.prob))
}


# Q3-1) Implement CART - sum of entropy in a particular node

myInfoEntropy<-function(cl.1,cl.2){
  n.t1<-length(cl.1) 
  n.t2<-length(cl.2)
  n.s<-n.t1+n.t2
  sum1.t1<-0
  for(i in cl.1){
    if(i==levels(cl.1)[1]) sum1.t1<-sum1.t1+1
  }
  sum1.t2<-0
  for(j in cl.2){
    if(j==levels(cl.2)[1]) sum1.t2<-sum1.t2+1
  }
  t1.e1<-sum1.t1/n.t1
  t1.e2<-1-t1.e1
  t2.e1<-sum1.t2/n.t2
  t2.e2<-1-t2.e1
  if(t1.e1==0||t1.e2==0) e.t1<-0 
  else e.t1<-(-1)*(t1.e1*log2(t1.e1)+t1.e2*log2(t1.e2))
  if(t2.e1==0||t2.e2==0) e.t2<-0
  else e.t2<-(-1)*(t2.e1*log2(t2.e1)+t2.e2*log2(t2.e2))
  sum.ent<-n.t1/n.s*e.t1+n.t2/n.s*e.t2
  return(sum.ent)
}

# Q3-2) Implement CART - split nodes into two children (based on lowest entropy)

mySplit<-function(dt,cl){
  new.dt<-cbind(dt,cl)
  min.max<-rbind(apply(dt,2,min),apply(dt,2,max))
  min<-10
  flag<-0
  for(j in 1:2){
    dt_ord<-new.dt[order(new.dt[,j]),]
    i<-1
    k<-1
    while(i<nrow(dt_ord)){
      cur<-dt_ord[i,j]
      while(dt_ord[i,j]<=cur) i<-i+1
      sub1<-dt_ord[j:i-1,]
      sub2<-dt_ord[i:nrow(dt_ord),]
      k<-i
      cl.1<-sub1[,3]
      cl.2<-sub2[,3]
      ent<-myInfoEntropy(cl.1,cl.2)
      if(flag==0||(flag==1&min>ent)){
        flag<-1
        min<-ent
        dt.left<-sub1
        dt.right<-sub2
      }
      if(dt_ord[i,j]==min.max[2,j]) break
    }
  }
  dt.left<-sub1[,1:2]
  cl.left<-sub1[,3]
  dt.right<-sub2[,1:2]
  cl.right<-sub2[,3]
  t1<-list(dt.left,cl.left)
  t2<-list(dt.right,cl.right)
  return (list(t1,t2))
}