# read data
path = "res/way1/"
power.cor = read.csv(paste0(path, "res_cor.csv"), row.names = 1)
power.g2m = read.csv(paste0(path, "res_g2m.csv"), row.names = 1)
power.g2t = read.csv(paste0(path, "res_g2t.csv"), row.names = 1)

png(paste0(path, "power_way1.png"), width = 1080, height = 1260, pointsize = 14)
par(mfrow = c(4,2), cex = 1)
xrange = (1:30)*0.5/30
titles = c("Linear", "Quadratic", "Cubic", "Radical","Sine: period 1/2", "Triangle", "Sine: period 1/8", "Step function")
for (i in 1:8)
{
  plot(xrange, power.cor[,i], ylim = c(0,1), main = titles[i], xlab = expression((1+sigma^2)^{-1}), ylab = "Power", pch = 1, col = "orange", type = 'o')
  lines(xrange, power.g2m[,i], pch = 2, col = "green", type = 'o')
  lines(xrange, power.g2t[,i], pch = 3, col = "blue", type = 'o')
  legend("topright",c("cor","g2m","g2t"), pch = c(1,2,3), col = c("orange","green","blue"))
}
dev.off()

plot(xrange, power.cor[,], ylim = c(0,1), main = "Linear", xlab = expression((1+sigma^2)^{-1}), ylab = "Power", pch = 1, col = "orange", type = 'o')
lines(xrange, power.g2m[,1], pch = 2, col = "green", type = 'o')
lines(xrange, power.g2t[,1], pch = 3, col = "blue", type = 'o')
legend("topright",c("cor","g2m","g2t"), pch = c(1,2,3), col = c("orange","green","blue"))


plot((1:30)/10, power.cor[2,], ylim = c(0,1), main = "Quadratic", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[2,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[2,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","MIC"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[3,], ylim = c(0,1), main = "Cubic", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[3,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[3,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","MIC"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[5,], ylim = c(0,1), main = "Sine: period 1/8", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[5,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[5,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","MIC"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[4,], ylim = c(0,1), main = "Sine: period 1/2", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[4,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[4,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","MIC"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[6,], ylim = c(0,1), main = "X^(1/4)", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[6,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[6,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","MIC"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[7,], ylim = c(0,1), main = "Circle", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[7,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[7,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","MIC"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[8,], ylim = c(0,1), main = "Step function", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[8,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[8,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","MIC"), pch = c(1,2,3), col = c("black","green","red"))
dev.off()