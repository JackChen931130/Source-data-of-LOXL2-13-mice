rm(list = ls())
Factorials<-function(n){
  init=1
  for (each in 1:n){
    init=init*each
  }
  return (init)
}#上面这个是阶乘的函数
Q_cal<-function(vec){#对任意由rank构成的样本排列向量的Q值函数
  #N<-length(vec)
  #vec<-na.omit(vec)#对没有rank的数据集进行剔除
  vec[is.na(vec)] <- 0#对没有rank的数据集进行0值填充
  n<-length(vec)
  V_vec<-rep(0,n)#先构造一个全0的长度为n的全零V向量，之后填充
  for (k in 1:n){
    sumnum=0#连加号的起点
    for (i in 1:k){#对连加号里的i在1:k上进行遍历
      first<-vec[i]*(-1)^(i-1)/Factorials(i)#首先是排除V(k-i)，这里的vec[i]就是r
      if (i!=k){
        first<-first*V_vec[k-i]#如果i不是k，那么V(k-i)不为V(0)
      }#如果是k则跳过，因为V(0)=1,直接用first即可
      sumnum<-sumnum+first#逐次累加first
    }
    V_vec[k]<-sumnum#把sumnum赋给V_vec第k个值，继续直至算到n
  }
  V<-V_vec[n]#V(n)
  Q<-Factorials(n)*V_vec[n]#最后用n!乘以V(n)即可
  return (Q)
  #return (V)
}
setwd("E:\\2021-10-排序算法")
try_csv1<-read.csv("test.csv")#测试样本
#先把第一列（sample）去掉，然后转置，这是因为dataframe的行无法
#变成向量（因为编写的Q_cal需要对每一个样本的行向量进行计算），
#所以变成列
rank_part<-t(try_csv1[-1])
sample<-ncol(rank_part)#样本数量
vec_Q<-rep(0,sample)
#vec_V<-rep(0,sample)
for (each in 1:sample){
  vec_Q[each]<-Q_cal(rank_part[,each])#累次计算每一个样本的Q值
}
print(vec_Q)#打印Q向量
try_csv1$Q1=vec_Q
###################################################################



#############生成rank矩阵##########################
rm(list = ls())
#setwd("E:\\2021-10-排序算法\\our data")
compare_name<-c("2021-10-Rank-mRNA.csv","2021-10-28-Rank-Protein-final.csv",
              "2021-10-28-Rank-Protein-phos-final.csv","2021-10-28-Rank-Protein-acety-final.csv",
              "2021-10-27-Rank-Pathway-genes-freq-final.csv")#所有的比对文件按顺序输入
#我看到你下面有个list的csv，要用它作为主文件进行其余五个文件的匹配，因此我单独将其拎出来
list_csv<-read.csv("2021-10-28-Rank-list.csv")
compare_ls<-list()
for (each in 1:5){
  compare_ls[[each]]<-read.csv(compare_name[each]) #list可以将每个dataframe储存
}
list_csv$index=1:nrow(list_csv)#确定原始list_csv的排序
for (each in 1:5){
  frame<-compare_ls[[each]]#去除比对frame里的每个dataframe
  #merge匹配，按照左边的csv为模板（即list_csv）,右边匹配到的y有的即有值，没有的就填充na
  x<-merge(list_csv,frame,by.x=names(list_csv)[1],by.y=names(frame)[1],all.x=TRUE)
  sort_index<-x$index#发现merge的时候顺序会乱，因此检查一下merge后的原始数据排序
  #x[,length(x)]即x的最后一列就是rank，那么each+1是list_csv的填充列，按照sort_index顺序填充
  list_csv[sort_index,each+1]=x[,length(x)]
}
total_csv<-list_csv[,-length(list_csv)]#去掉索引列
#这样得到的最终结果的symbol序列是正好按照原始list_csv排序的
#我改了一下row.names,这样保存csv可以不输出索引
write.table(total_csv,file="2021-10-28-Rank-all.csv",sep=",",row.names=FALSE)

#########################################排序算法################################

rm(list = ls())
Factorials<-function(n){
  init=1
  for (each in 1:n){
    init=init*each
  }
  return (init)
}#上面这个是阶乘的函数
Q_cal<-function(vec){#对任意由rank构成的样本排列向量的Q值函数
  N<-length(vec)
  vec<-na.omit(vec)#对没有rank的数据集进行剔除
  #vec[is.na(vec)] <- 0#对没有rank的数据集进行0值填充
  n<-length(vec)
  V_vec<-rep(0,n)#先构造一个全0的长度为n的全零V向量，之后填充
  for (k in 1:n){
    sumnum=0#连加号的起点
    for (i in 1:k){#对连加号里的i在1:k上进行遍历
      first<-vec[i]*(-1)^(i-1)/Factorials(i)#首先是排除V(k-i)，这里的vec[i]就是r
      if (i!=k){
        first<-first*V_vec[k-i]#如果i不是k，那么V(k-i)不为V(0)
      }#如果是k则跳过，因为V(0)=1,直接用first即可
      sumnum<-sumnum+first#逐次累加first
    }
    V_vec[k]<-sumnum#把sumnum赋给V_vec第k个值，继续直至算到n
  }
  V<-V_vec[n]#V(n)
  Q<-Factorials(N)*V_vec[n]#最后用n!乘以V(n)即可
  return (Q)
  #return (V)
}

#setwd("E:\\2021-10-排序算法\\our data")
data<-read.csv("2021-10-28-Rank-all.csv")#测试样本
data1<-data[rowSums(is.na(data))<4,]######只保留NA小于2个的蛋白(?是非NA大于两个吧)

#先把第一列（sample）去掉，然后转置，这是因为dataframe的行无法
#变成向量（因为编写的Q_cal需要对每一个样本的行向量进行计算），
#所以变成列
rank_part<-t(data1[,-1])#改了一下，便于我这边读取（因为第一列是symbol）
sample<-ncol(rank_part)#样本数量
vec_Q<-rep(NA,sample)
#vec_V<-rep(0,sample)
for (each in 1:sample){
  vec_Q[each]<-Q_cal(rank_part[,each])#累次计算每一个样本的Q值
}
#print(vec_Q)#打印Q向量
data1$Q1=vec_Q
write.table(data1,file="2021-10-28-Rank-all-final3.csv",sep=",",row.names=FALSE)
#vec_Q[1]<-Q_cal(124,0,0)

