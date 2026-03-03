# ==============================
# 扩增子测序 稳定性分析 (Fig.4A + Fig.4B, SC 和 SH 都保留)
# ==============================

library(tidyverse)
library(vegan)
library(patchwork)

# === 1. 导入数据 ===
asv_table <- read.csv("ASV_abundance.csv", row.names = 1)
sample_info <- read.csv("sample_info.csv")

asv_long <- asv_table %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(sample_info, by = "SampleID")

# === 2. 建立 时间映射 (请根据实验设计修改小时数) ===
time_map <- c(
  "DA0"  = 0,
  "DB3"  = 72,
  "DC7"  = 168,
  "DD14" = 336,
  "DE21" = 504,
  "DF28" = 672,
  "DG60" = 1440
)

sample_info <- sample_info %>%
  mutate(Time_num = time_map[Time])

# === 3. 定义参考状态 (这里统一用 SC 组的平均值作为参考) ===
ref_samples <- sample_info %>% 
  filter(Group == "SC", Time %in% names(time_map))

ref_abund <- asv_long %>%
  filter(SampleID %in% ref_samples$SampleID) %>%
  group_by(ASV) %>%
  summarise(Ref = mean(Abundance), .groups = "drop")

# === 4. 偏离计算函数 ===
calc_dev <- function(sample_id, method = "euclidean"){
  abund <- asv_long %>% filter(SampleID == sample_id) %>% select(ASV, Abundance)
  df <- left_join(ref_abund, abund, by = "ASV") %>% replace_na(list(Abundance = 0))
  dist(rbind(df$Ref, df$Abundance), method = method)[1]
}

# === 5. 计算每个 Group 的偏离 (不去掉 SC) ===
deviation_df <- sample_info %>%
  mutate(
    Euclidean = sapply(SampleID, calc_dev, method = "euclidean"),
    Canberra  = sapply(SampleID, calc_dev, method = "canberra")
  ) %>%
  arrange(Group, Time_num)

# === 6. 聚合轨迹 (每个 Group 每个时间点取均值) ===
traj_df <- deviation_df %>%
  group_by(Group, Time, Time_num) %>%
  summarise(
    Euclidean = mean(Euclidean, na.rm = TRUE),
    Canberra  = mean(Canberra,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Group, Time_num)

# === 7. 标准化 (缩放到 0~0.5 区间) ===
scale_distance <- function(x) (x - min(x)) / (max(x) - min(x))
traj_df <- traj_df %>%
  group_by(Group) %>%
  mutate(
    Euclidean = scale_distance(Euclidean),
    Canberra  = scale_distance(Canberra)
  ) %>%
  ungroup()

# === 8. 干扰点 ===
disturb_time <- "DB3"

# === 9. 循环按 Group 绘图 & 指标计算 ===
groups <- unique(traj_df$Group)
all_metrics <- list()
plot_list <- list()

for (g in groups) {
  df_g <- traj_df %>% filter(Group == g)
  disturb_index <- match(disturb_time, df_g$Time)
  
  eu_vec <- df_g$Euclidean
  ca_vec <- df_g$Canberra
  
  # === 指标计算 ===
  RS_Euclidean <- eu_vec[disturb_index]
  RS_Canberra  <- ca_vec[disturb_index]
  
  DS_Euclidean <- if (disturb_index > 1) eu_vec[disturb_index] - eu_vec[disturb_index-1] else NA
  DS_Canberra  <- if (disturb_index > 1) ca_vec[disturb_index] - ca_vec[disturb_index-1] else NA
  
  RL_Euclidean <- RS_Euclidean - min(eu_vec[disturb_index:length(eu_vec)], na.rm = TRUE)
  RL_Canberra  <- RS_Canberra  - min(ca_vec[disturb_index:length(ca_vec)], na.rm = TRUE)
  
  rec_idx_eu <- which.min(eu_vec[disturb_index:length(eu_vec)]) + disturb_index - 1
  rec_idx_ca <- which.min(ca_vec[disturb_index:length(ca_vec)]) + disturb_index - 1
  
  E_Euclidean <- RL_Euclidean / max(1, rec_idx_eu - disturb_index + 1)
  E_Canberra  <- RL_Canberra  / max(1, rec_idx_ca - disturb_index + 1)
  
  metrics_g <- tibble(
    Group = g,
    Distance = c("Euclidean", "Canberra"),
    RS = c(RS_Euclidean, RS_Canberra),
    DS = c(DS_Euclidean, DS_Canberra),
    RL = c(RL_Euclidean, RL_Canberra),
    E  = c(E_Euclidean,  E_Canberra)
  )
  
  all_metrics[[g]] <- metrics_g
  
  # === 绘图 ===
  # Fig.4A 极坐标（Time 作为因子 → 等间距）
  df_g <- df_g %>%
    mutate(Time_factor = factor(Time, levels = names(time_map)))
  
  polar_g <- ggplot(df_g, aes(x = Time_factor)) +
    geom_path(aes(y = Euclidean, colour = "Euclidean", group = 1), size = 1) +
    geom_point(aes(y = Euclidean, colour = "Euclidean"), size = 2) +
    geom_path(aes(y = Canberra,  colour = "Canberra", group = 1),  size = 1) +
    geom_point(aes(y = Canberra,  colour = "Canberra"),  size = 2) +
    coord_polar(theta = "x") +
    scale_color_manual(values = c("Euclidean"="blue","Canberra"="red")) +
    labs(y = "Deviation (d)", x = "Time", 
         title = paste("Fig.4A Deviation Trajectory - Group", g)) +
    theme_minimal(base_size=14)
  
  # Fig.4B 恢复曲线（数值时间 → 保留真实间隔）
  baseline_eu <- df_g$Euclidean[df_g$Time == disturb_time]
  baseline_ca <- df_g$Canberra[df_g$Time == disturb_time]
  
  df_g <- df_g %>%
    mutate(
      Euclidean_res = Euclidean - baseline_eu,
      Canberra_res  = Canberra  - baseline_ca
    )
  
  eu_plot <- ggplot(df_g, aes(Time_num, Euclidean_res)) +
    geom_path(color="blue", size=1) + geom_point(color="blue", size=2) +
    geom_hline(yintercept = 0, linetype="dashed") +
    scale_x_continuous(breaks = time_map, labels = names(time_map)) +
    labs(x="Time", y="Resilience_Euclidean", 
         title=paste("Fig.4B Euclidean - Group", g)) +
    theme_bw(base_size=14)
  
  ca_plot <- ggplot(df_g, aes(Time_num, Canberra_res)) +
    geom_path(color="red", size=1) + geom_point(color="red", size=2) +
    geom_hline(yintercept = 0, linetype="dashed") +
    scale_x_continuous(breaks = time_map, labels = names(time_map)) +
    labs(x="Time", y="Resilience_Canberra", 
         title=paste("Fig.4B Canberra - Group", g)) +
    theme_bw(base_size=14)
  
  final_g <- polar_g / (eu_plot + ca_plot) +
    plot_annotation(title = paste("Stability Analysis - Group", g),
                    caption = paste0(
                      "Euclidean: RS=", round(RS_Euclidean,3),
                      " | DS=", round(DS_Euclidean,3),
                      " | RL=", round(RL_Euclidean,3),
                      " | E=", round(E_Euclidean,3), "\n",
                      "Canberra: RS=", round(RS_Canberra,3),
                      " | DS=", round(DS_Canberra,3),
                      " | RL=", round(RL_Canberra,3),
                      " | E=", round(E_Canberra,3)
                    ))
  
  plot_list[[g]] <- final_g
  
  # 保存每个 Group 的图
  ggsave(paste0("Stability_", g, ".png"), final_g, width=10, height=12)
}

# === 10. 汇总所有 Group 的指标表 ===
stability_metrics_all <- bind_rows(all_metrics)
write.csv(stability_metrics_all, "Stability_metrics_SC_SH.csv", row.names = FALSE)

# === 11. 拼接所有 Group 的图在一张大图中 ===
final_all <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(title = "Stability Analysis - All Groups")
ggsave("Stability_All.png", final_all, width=20, height=12)

# === 12. 查看某个 Group 的图 (比如 SC) ===
plot_list[["SC"]]
plot_list[["SH"]]
