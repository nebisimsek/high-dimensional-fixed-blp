
#Get the RData Files

# all.files <- c("0.2_char10_ins10_rel5.RData", "0.2_char10_ins20_rel5.RData",
#                "0.2_char10_ins50_rel5.RData","0.2_char50_ins10_rel5.RData", 
#                "0.2_char50_ins20_rel5.RData", "0.2_char50_ins50_rel5.RData")

all.files <- c("endo0.8_0.2_char10_ins10_rel5_lambda0.1.RData",
               "endo0.8_0.2_char10_ins50_rel5_lambda0.1.RData", 
               "endo0.8_0.2_char10_ins50_rel5_lambda0.2.RData")


#Take specific value from the Rdata files
mylist <- lapply(all.files, function(x) {
  load(file = x)
  get(ls()[ls() == "p.values"])
})

#Name Rdata files
# names(mylist) <- c("endo0.8_0.2_char10_ins10_rel5_lambda0.1.RData",
#                     "endo0.8_0.2_char10_ins50_rel5_lambda0.1.RData", 
#                     "endo0.8_0.2_char10_ins50_rel5_lambda0.2.RData")

#Get the means of specified data
means <- lapply(mylist, colMeans)

#Get the bias 
lapply(means, function(x) abs(-5-x))
write.table(lapply(means, function(x) abs(-5-x)))


#Get the standard error
rapply(mylist, function(x) apply(x, 2, sd), how = "replace")
write.table(rapply(mylist, function(x) apply(x, 2, sd), how = "replace")
)


#Get the Root Means squared error
rapply(mylist, function(x) apply(x, 2, function(x) sqrt(mean((x+5)^2))), how = "replace")
write.table(rapply(mylist, function(x) apply(x, 2, function(x) sqrt(mean((x+5)^2))), how = "replace")
)


# mylist.z <- lapply(all.files, function(x) {
#   load(file = x)
#   get(ls()[ls() == "z.count"])
# })
# 
# lapply(mylist.z, colMeans)
# 
# 
# mylist.z.select <- lapply(all.files, function(x) {
#   load(file = x)
#   get(ls()[ls() == "z.selection"])
# })
# 
# lapply(mylist.z.select, colMeans)
