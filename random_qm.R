require(stringr)

#uniform distributed questionmarks in the whole dataset
data = read.table("tb_data_9.txt",colClasses = c("character", "character"))
L = 9
threshold = 0.5

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

write.table(data, "tb_data_9_qm50.txt", row.names=FALSE, quote=FALSE)

#uniform distributed questionmarks at a specific position
data.2 = read.table("full.txt", colClasses = c("character","character"))
L = 6
feature = 6
threshold = 0.1

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

write.table(data.2, "feature6_10.txt", row.names=FALSE, quote=FALSE)