library(ggplot2)
df <- read.csv("./benchmatmul.csv")
p <- ggplot(data = df, aes(x=Mem, y=Time, color=optim, group=optim))+geom_line()
p <- p+facet_grid(device ~ optim)+labs(y="Backward error")+labs(x="Solver iter.")
p <- p + ggtitle("matmul") + labs(y = "Time (ms)", x = "Memory (Mb)")
ggsave("benchmatmul.pdf", plot=p)

#p <- p+facet_grid(Sthetac ~ Sic+Swhichc+l)+labs(y="Backward error")+labs(x="Solver iter.")