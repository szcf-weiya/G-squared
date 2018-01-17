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
