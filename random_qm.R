require(stringr)

#uniform distributed questionmarks in the whole dataset
data = read.table("tb_data_10.txt",colClasses = c("character", "character"))
L = 10
threshold = 0.5   #approximately this amount of positions will be replaced by a ?

length = nrow(data)
for (i in 1:length){
  dp1 = data$V1[i]
  dp2 = data$V2[i]
  for (j in 1:L){
    rn = runif(2,0,1)
    if (rn[1] < threshold){
      str_sub(dp1,j,j) = "?"
      data$V1[i] = dp1
    }
    if (rn[2] < threshold){
      str_sub(dp2,j,j) = "?"
      data$V2[i] = dp2
    }
  }
}

write.table(data, "tb_data_10_qm50.txt", row.names=FALSE, quote=FALSE)

#uniform distributed questionmarks at a specific position
data.2 = read.table("full.txt", colClasses = c("character","character"))
L = 6
feature = 1   #specify the position in which the uncertainty markers will be inserted
threshold = 0.4 #in approximately this amount of strings the specified position will be replaced by a ?

length = nrow(data.2)
for (i in 1:length){
  dp1 = data.2$V1[i]
  dp2 = data.2$V2[i]
  rn = runif(2,0,1)
  if (rn[1] < threshold){
    str_sub(dp1,feature,feature) = "?"
    data.2$V1[i] = dp1
  }
  if (rn[2] < threshold){
    str_sub(dp2, feature, feature) = "?"
    data.2$V2[i] = dp2
  }
}

write.table(data.2, "feature1_40.txt", row.names=FALSE, quote=FALSE)
