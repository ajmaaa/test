library(tidyverse)
library(vegan)
library(patchwork)

# === 数据准备 ===
asv_table <- read.csv("ASV_abundance.csv", row.names = 1)
sample_info <- read.csv("sample_info.csv")

asv_long <- asv_table %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(sample_info, by = "SampleID")

ref_samples <- sample_info %>% 
  filter(Group == "SC", Time %in% c("DA0","DB3","DC7","DD14","DE21","DF28","DG60"))

ref_abund <- asv_long %>%
  filter(SampleID %in% ref_samples$SampleID) %>%
  group_by(ASV) %>%
  summarise(Ref = mean(Abundance), .groups = "drop")

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
  arrange(Time) %>%
  mutate(Time = factor(Time, levels = c("DA0","DB3","DC7","DD14","DE21","DF28","DG60")))

# === 干扰事件索引 ===
disturb_time <- "DB3"
disturb_index <- which(levels(deviation_df$Time) == disturb_time)

# 将偏离向量提取出来
eu_vec <- deviation_df$Euclidean
ca_vec <- deviation_df$Canberra

# === 计算稳定性指标 ===
# RS: 最大偏离（干扰点）
RS_Euclidean <- eu_vec[disturb_index]
RS_Canberra  <- ca_vec[disturb_index]

# DS: 干扰瞬时偏离速度
DS_Euclidean <- eu_vec[disturb_index] - eu_vec[disturb_index-1]
DS_Canberra  <- ca_vec[disturb_index] - ca_vec[disturb_index-1]

# RL: 干扰后恢复能力 = 干扰点 - 最小偏离
RL_Euclidean <- RS_Euclidean - min(eu_vec[disturb_index:length(eu_vec)])
RL_Canberra  <- RS_Canberra - min(ca_vec[disturb_index:length(ca_vec)])

# E: 弹性 = RL / 恢复所需时间
recovery_index_eu <- which.min(eu_vec[disturb_index:length(eu_vec)]) + disturb_index - 1
recovery_index_ca <- which.min(ca_vec[disturb_index:length(ca_vec)]) + disturb_index - 1

E_Euclidean <- RL_Euclidean / (recovery_index_eu - disturb_index + 1)
E_Canberra  <- RL_Canberra  / (recovery_index_ca - disturb_index + 1)

# 汇总
stability_metrics <- tibble(
  Metric = c("RS", "DS", "RL", "E"),
  Euclidean = c(RS_Euclidean, DS_Euclidean, RL_Euclidean, E_Euclidean),
  Canberra  = c(RS_Canberra, DS_Canberra, RL_Canberra, E_Canberra)
)

print(stability_metrics)

# === 绘图 ===
euclidean_plot <- ggplot(deviation_df, aes(x=Time)) +
  geom_line(aes(y=Euclidean), size=1, colour="blue") +
  geom_point(aes(y=Euclidean), size=2, colour="blue") +
  geom_vline(xintercept = disturb_index, linetype="dashed", color="red") +
  labs(x="Time", y="Euclidean deviation", title="Euclidean Distance") +
  theme_bw(base_size=14)

canberra_plot <- ggplot(deviation_df, aes(x=Time)) +
  geom_line(aes(y=Canberra), size=1, colour="red") +
  geom_point(aes(y=Canberra), size=2, colour="red") +
  geom_vline(xintercept = disturb_index, linetype="dashed", color="red") +
  labs(x="Time", y="Canberra deviation", title="Canberra Distance") +
  theme_bw(base_size=14)

combined_plot <- euclidean_plot + canberra_plot +
  plot_layout(ncol=1) +
  plot_annotation(
    title = "Bioreactor Stability Analysis",
    caption = paste0(
      "Euclidean: RS=", round(RS_Euclidean,3), 
      " | DS=", round(DS_Euclidean,3),
      " | RL=", round(RL_Euclidean,3),
      " | E=", round(E_Euclidean,3), "\n",
      "Canberra: RS=", round(RS_Canberra,3),
      " | DS=", round(DS_Canberra,3),
      " | RL=", round(RL_Canberra,3),
      " | E=", round(E_Canberra,3)
    )
  )

print(combined_plot)


# 先创建一个数据框保存指标
stability_metrics <- tibble(
  TimePoint = disturb_time,
  Distance = c("Euclidean", "Canberra"),
  RS = c(RS_Euclidean, RS_Canberra),
  DS = c(DS_Euclidean, DS_Canberra),
  RL = c(RL_Euclidean, RL_Canberra),
  E  = c(E_Euclidean,  E_Canberra)
)

print(stability_metrics)
