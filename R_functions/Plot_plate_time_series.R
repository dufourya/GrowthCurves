fileName = 'magellan_sample.asc'
magellan <- read.table(fileName, header = T, sep = ',',strip.white = T, stringsAsFactors = FALSE, fileEncoding = "latin1")
timeStamps <- as.numeric(gsub("[^0-9]", "", magellan[1,-(1:which(colnames(magellan) == "Raw.data") - 1)]))
temperature <- as.numeric(gsub("[^0-9\\.]", "", magellan[2,-(1:which(colnames(magellan) == "Raw.data") - 1)]))
wellPos <- as.vector(magellan[-(1:2),1])
methodName <- as.vector(magellan[-(1:2),2])
dilutionFactor <- as.vector(magellan[-(1:2),3])
plateLayout <- as.vector(magellan[-(1:2),4])
sampleReplicate <- as.vector(magellan[-(1:2),5])
data <- data.matrix(magellan[-(1:2),-(1:which(colnames(magellan) == "Raw.data") - 1)], rownames.force = NA)
rownames(data) <- wellPos
colnames(data) <- timeStamps

##
par(mfrow = c(8,12))
par(mar = c(0,0,0,0))
for (i in 1:length(wellPos)) {
  plot(timeStamps, data[i,],xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(min(data,na.rm = T), max(data,na.rm = T)), type = 'l',log = 'y')
  text(max(timeStamps,na.rm = T)*0.1, max(data,na.rm = T)*0.8, plateLayout[i], adj = c(0,0))
}
##

growthRate <- t(diff(t(log2(data)),lag = 3))
par(mfrow = c(8,12))
par(mar = c(0,0,0,0))
for (i in 1:length(wellPos)) {
  plot(timeStamps[1:length(growthRate[i,])], growthRate[i,],xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(-0.1, max(growthRate,na.rm = T)), type = 'l')
  text(max(timeStamps,na.rm = T)*0.5, max(growthRate,na.rm = T)*0.8, plateLayout[i], adj = c(0,0))
}

##
growthRate <- t(diff(t(log2(data)),lag = 3))
par(mfrow = c(8,12))
par(mar = c(0,0,0,0))
for (i in 1:length(wellPos)) {
  hist(growthRate[i,], col = 'black',xlab = '', ylab = '', yaxt = 'n', main ='', xlim=c(-0.1,0.6),breaks = seq(from = min(growthRate,na.rm=T)-0.025, to = max(growthRate,na.rm=T)+0.025, by = 0.01))
}
##

par(op)
par(mfrow=c(2,4), bty="n", cex = 0.9, mgp = c(3, 1, 0))
plot(timeStamps/60, data[12+12,], ylim = c(min(data,na.rm = T), max(data,na.rm = T)), axes = F,lwd = 2, type = 'l', main='P. aeruginosa (ancestor)',log="y",xlab="Time (min)", ylab="OD590")
axis(side = 1, las = 0)
axis(side = 2, las = 2)

plot(timeStamps/60, data[12+3,], ylim = c(min(data,na.rm = T), max(data,na.rm = T)), axes = F,lwd = 2, type = 'l',col='red', main='P. aeruginosa pos. A',log="y",xlab="Time (min)", ylab="OD590")
points(timeStamps/60, data[12+4,], lwd = 2, type = 'l',col='blue')
points(timeStamps/60, data[12+5,], lwd = 2, type = 'l',col='green')
axis(side = 2, las = 2)
axis(side = 1, las = 0)

plot(timeStamps/60, data[12+6,], ylim = c(min(data,na.rm = T), max(data,na.rm = T)), axes = F,lwd = 2, type = 'l',col='red', main='P. aeruginosa pos. B',log="y",xlab="Time (min)", ylab="OD590")
points(timeStamps/60, data[12+7,], lwd = 2, type = 'l',col='blue')
points(timeStamps/60, data[12+8,], lwd = 2, type = 'l',col='green')
axis(side = 2, las = 2)
axis(side = 1, las = 0)

plot(timeStamps/60, data[12+9,], ylim = c(min(data,na.rm = T), max(data,na.rm = T)), axes = F,lwd = 2, type = 'l',col='red', main='P. aeruginosa pos. C',log="y",xlab="Time (min)", ylab="OD590")
points(timeStamps/60, data[12+10,], lwd = 2, type = 'l',col='blue')
points(timeStamps/60, data[12+11,], lwd = 2, type = 'l',col='green')
axis(side = 2, las = 2)
axis(side = 1, las = 0)

plot(timeStamps/60, data[36+12,], ylim = c(min(data,na.rm = T), max(data,na.rm = T)), axes = F,lwd = 2, type = 'l', main='V. cholerae (ancestor)',log="y",xlab="Time (min)", ylab="OD590")
axis(side = 2, las = 2)
axis(side = 1, las = 0)

plot(timeStamps/60, data[36+3,], ylim = c(min(data,na.rm = T), max(data,na.rm = T)), axes = F,lwd = 2, type = 'l',col='red', main='V. cholerae pos. A',log="y",xlab="Time (min)", ylab="OD590")
points(timeStamps/60, data[36+4,], lwd = 2, type = 'l',col='blue')
points(timeStamps/60, data[36+5,], lwd = 2, type = 'l',col='green')
axis(side = 2, las = 2)
axis(side = 1, las = 0)

plot(timeStamps/60, data[36+6,], ylim = c(min(data,na.rm = T), max(data,na.rm = T)), axes = F,lwd = 2, type = 'l',col='red', main='V. cholerae pos. B',log="y",xlab="Time (min)", ylab="OD590")
points(timeStamps/60, data[36+7,], lwd = 2, type = 'l',col='blue')
points(timeStamps/60, data[36+8,], lwd = 2, type = 'l',col='green')
axis(side = 2, las = 2)
axis(side = 1, las = 0)

plot(timeStamps/60, data[36+9,], ylim = c(min(data,na.rm = T), max(data,na.rm = T)), axes = F,lwd = 2, type = 'l',col='red', main='V. cholerae pos. C',log="y",xlab="Time (min)", ylab="OD590")
points(timeStamps/60, data[36+10,], lwd = 2, type = 'l',col='blue')
points(timeStamps/60, data[36+11,], lwd = 2, type = 'l',col='green')
axis(side = 2, las = 2)
axis(side = 1, las = 0)


resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}