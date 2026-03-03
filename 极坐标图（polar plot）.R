library(ggplot2)
library(dplyr)
library(tidyr)

# 模拟数据
df <- data.frame(
  Time = rep(c(0, 152, 248, 344), 3),   # 时间点
  Reactor = rep(c("A","B","C"), each=4),
  Distance = c(0.1, 0.2, 0.3, 0.1, 0.05, 0.15, 0.25, 0.1, 0.08, 0.18, 0.28, 0.12),
  Method = rep("Canberra", 12)
)

# 极坐标图
ggplot(df, aes(x=Time, y=Distance, color=Reactor, group=Reactor)) +
  geom_line(size=1) +
  geom_point(aes(shape=Reactor), size=3) +
  coord_polar(theta="x") +
  labs(
    title="MDb-MM Stability Properties (Resilience / Deviation)",
    x="Time (h)",
    y="Deviation from reference state",
    color="Bioreactor",
    shape="Bioreactor"
  ) +
  theme_minimal(base_size=14)

