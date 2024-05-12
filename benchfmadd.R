library(ggplot2)
df <- read.csv("./benchfmadd.csv")
p <- ggplot(data = df, aes(x=Mem, y=Time, color=optim, group=optim))+geom_line()
p <- p + ggtitle("fmadd") + labs(y = "Time (ms)", x = "Memory (Mb)")
ggsave("bench.pdf", plot=p)

transform <- subset(df, optim == "transform")
mipp <- subset(df, optim == "mipp")
for(i in 1:length(mipp$type)) {
    mipp$TotSpeedUp[i] = transform$Total[i] / mipp$Total[i]
    mipp$SpeedUp[i] = transform$Time[i] / mipp$Time[i]
    print(mipp$SpeedUp[i])
}
p <- ggplot(data = mipp, aes(x=Mem, y=SpeedUp))+geom_line()
p <- p + ggtitle("fmadd") + labs(y = "SpeedUp", x = "Memory (Mb)")
ggsave("SpeedUp.pdf", plot=p)
#p <- p+facet_grid(Sthetac ~ Sic+Swhichc+l)+labs(y="Backward error")+labs(x="Solver iter.")