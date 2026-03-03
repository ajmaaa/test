library(tidyverse)
library(vegan)
library(patchwork)

# 读取数据
asv_table <- read.csv("ASV_abundance.csv", row.names = 1)
sample_info <- read.csv("sample_info.csv")

# 转换为长格式
asv_long <- asv_table %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(sample_info, by = "SampleID")

# 参考状态 (对照组平均)
ref_samples <- sample_info %>% 
  filter(Group == "SC", Time %in% c("DA0","DB3","DC7","DD14","DE21","DF28","DG60"))

ref_abund <- asv_long %>%
  filter(SampleID %in% ref_samples$SampleID) %>%
  group_by(ASV) %>%
  summarise(Ref = mean(Abundance), .groups = "drop")

# 计算每个样本与参考状态的距离
calc_dev <- function(sample_id, method = "euclidean"){
  abund <- asv_long %>% filter(SampleID == sample_id) %>% select(ASV, Abundance)
  df <- left_join(ref_abund, abund, by = "ASV") %>% replace_na(list(Abundance = 0))
  dist(rbind(df$Ref, df$Abundance), method = method)[1]
}

deviation_df <- sample_info %>%
  mutate(
    Euclidean = sapply(SampleID, calc_dev, method = "euclidean"),
    Canberra  = sapply(SampleID, calc_dev, method = "canberra")
  ) %>%
  filter(Group == "SH") %>%
  mutate(Time = factor(Time, levels = c("DA0","DB3","DC7","DD14","DE21","DF28","DG60")))

# 计算RS和DS
deviation_df <- deviation_df %>%
  arrange(Time) %>%
  mutate(
    RS_Euclidean = max(Euclidean),
    DS_Euclidean = Euclidean - lag(Euclidean, default = first(Euclidean)),
    RS_Canberra = max(Canberra),
    DS_Canberra = Canberra - lag(Canberra, default = first(Canberra))
  )

# 干扰事件时间点
disturb_time <- "DB3"
disturb_index <- which(levels(deviation_df$Time) == disturb_time)

# 绘制Euclidean图
euclidean_plot <- ggplot(deviation_df, aes(x=Time)) +
  geom_line(aes(y=Euclidean), size=1, colour="blue") +
  geom_point(aes(y=Euclidean), size=2, colour="blue") +
  geom_segment(
    data = deviation_df %>% filter(!is.na(lead(Euclidean))),
    aes(x=Time, xend=lead(Time), y=Euclidean, yend=lead(Euclidean)),
    color="blue", size=0.5
  ) +
  geom_vline(xintercept = disturb_index, linetype="dashed", color="red") +
  labs(x="Time", y="Euclidean deviation", title="Euclidean Distance") +
  theme_bw(base_size=14) +
  theme(legend.position="none")

# 绘制Canberra图
canberra_plot <- ggplot(deviation_df, aes(x=Time)) +
  geom_line(aes(y=Canberra), size=1, colour="red") +
  geom_point(aes(y=Canberra), size=2, colour="red") +
  geom_segment(
    data = deviation_df %>% filter(!is.na(lead(Canberra))),
    aes(x=Time, xend=lead(Time), y=Canberra, yend=lead(Canberra)),
    color="red", size=0.5
  ) +
  geom_vline(xintercept = disturb_index, linetype="dashed", color="red") +
  labs(x="Time", y="Canberra deviation", title="Canberra Distance") +
  theme_bw(base_size=14) +
  theme(legend.position="none")

# 组合图
combined_plot <- euclidean_plot + canberra_plot +
  plot_layout(ncol=1) +
  plot_annotation(
    title = "Bioreactor Stability Analysis",
    caption = paste0(
      "Euclidean: RS = ", round(max(deviation_df$Euclidean),3),
      " | DS = ", round(deviation_df$DS_Euclidean[which(deviation_df$Time==disturb_time)],3), "\n",
      "Canberra: RS = ", round(max(deviation_df$Canberra),3),
      " | DS = ", round(deviation_df$DS_Canberra[which(deviation_df$Time==disturb_time)],3)
    )
  )

print(combined_plot)
