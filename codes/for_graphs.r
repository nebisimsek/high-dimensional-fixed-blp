#Get the Rdata files

all.files <- c("0.2_char10_ins10_rel5.RData",
               "0.2_char10_ins50_rel5.RData", 
               "0.2_char50_ins10_rel5.RData",
               "0.2_char50_ins50_rel5.RData")

#Get the specified data from the RData files
mylist <- lapply(all.files, function(x) {
  load(file = x)
  get(ls()[ls() == "p.values"])
})

mylist2 <- lapply(all.files, function(x) {
  load(file = x)
  get(ls()[ls() == "x.count"])
})

#some specifications for plot

par(mfrow = c(2,2))
op <- par(cex = 0.5)

p.values <- mylist[[1]][]

dens1 <- density(p.values[,'oracle'])
dens2 <- density(p.values[,'full_model'])
dens3 <- density(p.values[,'triple-lasso'])
dens4 <- density(p.values[,'only_first'])
dens5 <- density(p.values[,'wo_third'])


plot(density(p.values[,'oracle']), xlim = range(p.values[,'oracle'], 
                                                p.values[,'full_model'],
                                                p.values[,'triple-lasso'],
                                                p.values[,'only_first'],
                                                p.values[,'wo_third']),
     ylim = range(dens1$y,dens2$y,dens3$y, dens4$y, dens5$y), 
     col = 'ivory4', lwd = 2, lty = 2,
     xlab = "", main = "",sub = "K = 10, L = 10", ylab = "")
abline(v = -5, col = 'black')
lines(density(p.values[,'full_model']), col = 'red', lwd = 2, lty = 2)
lines(density(p.values[,'triple-lasso']), col = 'blue', lwd = 2)
lines(density(p.values[,'only_first']), col = 'lightgreen', lwd = 2, lty = 5)
#lines(density(p.values[,'wo_second']), col = 'darkmagenta', lwd = 2)
lines(density(p.values[,'wo_third']), col = 'goldenrod1', lwd = 2,lty = 1)

legend("topright",legend = c("Oracle", "Unpenalized", 
                             "Single-Lasso", "Double-Lasso","Triple-Lasso"),
       col = c('ivory4', 'red', 'lightgreen',' goldenrod1', 'blue'),
       lty = c(2, 2, 2, 1, 1), lwd = 2)


p.values <- mylist[[2]][]

dens1 <- density(p.values[,'oracle'])
dens2 <- density(p.values[,'full_model'])
dens3 <- density(p.values[,'triple-lasso'])
dens4 <- density(p.values[,'only_first'])
dens5 <- density(p.values[,'wo_third'])


plot(density(p.values[,'oracle']), xlim = range(p.values[,'oracle'], 
                                                p.values[,'full_model'],
                                                p.values[,'triple-lasso'],
                                                p.values[,'only_first'],
                                                p.values[,'wo_third']),
     ylim = range(dens1$y,dens2$y,dens3$y, dens4$y, dens5$y), 
     col = 'ivory4', lwd = 2, lty = 2,
     xlab = "", main = "",sub = "K = 10, L = 50", ylab = "")
abline(v = -5, col = 'black')
lines(density(p.values[,'full_model']), col = 'red', lwd = 2, lty = 2)
lines(density(p.values[,'triple-lasso']), col = 'blue', lwd = 2)
lines(density(p.values[,'only_first']), col = 'lightgreen', lwd = 2, lty = 5)
#lines(density(p.values[,'wo_second']), col = 'darkmagenta', lwd = 2)
lines(density(p.values[,'wo_third']), col = 'goldenrod1', lwd = 2,lty = 1)

legend("topright",legend = c("Oracle", "Unpenalized", 
                             "Single-Lasso", "Double-Lasso","Triple-Lasso"),
       col = c('ivory4', 'red', 'lightgreen',' goldenrod1', 'blue'),
       lty = c(2, 2, 2, 1, 1), lwd = 2)


p.values <- mylist[[3]][]
x.count <- mylist2[[3]][]

#Get the index of not identified cases
#notid <- which(x.count[,1]*4 + x.count[,2] > 19,)
#Drop them from the values
#p.values <- p.values[-notid,]


dens1 <- density(p.values[,'oracle'])
dens3 <- density(p.values[,'triple-lasso'])
dens4 <- density(p.values[,'only_first'])
dens5 <- density(p.values[,'wo_third'])


plot(density(p.values[,'oracle']), xlim = range(p.values[,'oracle'],
                                                p.values[,'triple-lasso'],
                                                p.values[,'only_first'],
                                                p.values[,'wo_third']),
     ylim = range(dens1$y,dens3$y, dens4$y, dens5$y), 
     col = 'ivory4', lwd = 2, lty = 2,
     xlab = "", main = "",sub = "K = 50, L = 10", ylab = "")
abline(v = -5, col = 'black')
#lines(density(p.values[,'full_model']), col = 'red', lwd = 2, lty = 2)
lines(density(p.values[,'triple-lasso']), col = 'blue', lwd = 2)
lines(density(p.values[,'only_first']), col = 'lightgreen', lwd = 2, lty = 5)
#lines(density(p.values[,'wo_second']), col = 'darkmagenta', lwd = 2)
lines(density(p.values[,'wo_third']), col = 'goldenrod1', lwd = 2,lty = 1)

legend("topright",legend = c("Oracle",
                             "Single-Lasso", "Double-Lasso","Triple-Lasso"),
       col = c('ivory4', 'lightgreen',' goldenrod1', 'blue'),
       lty = c(2, 2, 1, 1), lwd = 2)



p.values <- mylist[[4]][]
x.count <- mylist2[[4]][]

#Get the index of not identified cases
#notid <- which(x.count[,1]*4 + x.count[,2] > 19,)
#Drop them from the values
#p.values <- p.values[-notid,]

dens1 <- density(p.values[,'oracle'])
dens3 <- density(p.values[,'triple-lasso'])
dens4 <- density(p.values[,'only_first'])
dens5 <- density(p.values[,'wo_third'])


plot(density(p.values[,'oracle']), xlim = range(p.values[,'oracle'],
                                                p.values[,'triple-lasso'],
                                                p.values[,'only_first'],
                                                p.values[,'wo_third']),
     ylim = range(dens1$y,dens3$y, dens4$y, dens5$y),
     col = 'ivory4', lwd = 2, lty = 2,
     xlab = "", main = "",sub = "K = 50, L = 50", ylab = "")
abline(v = -5, col = 'black')
#lines(density(p.values[,'full_model']), col = 'red', lwd = 2, lty = 2)
lines(density(p.values[,'triple-lasso']), col = 'blue', lwd = 2)
lines(density(p.values[,'only_first']), col = 'lightgreen', lwd = 2, lty = 5)
#lines(density(p.values[,'wo_second']), col = 'darkmagenta', lwd = 2)
lines(density(p.values[,'wo_third']), col = 'goldenrod1', lwd = 2,lty = 1)

legend("topright",legend = c("Oracle",
                             "Single-Lasso", "Double-Lasso","Triple-Lasso"),
       col = c('ivory4', 'lightgreen',' goldenrod1', 'blue'),
       lty = c(2, 2, 1, 1), lwd = 2)





