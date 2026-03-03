library(tidyverse)

library(tidyverse)
library(vegan)   # 用于距离计算

# 读取数据
asv_table <- read.csv("ASV_abundance.csv", row.names = 1)
sample_info <- read.csv("sample_info.csv")

# 转换为长格式
asv_long <- asv_table %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to="SampleID", values_to="Abundance") %>%
  left_join(sample_info, by="SampleID")

# 计算参考状态 (对照组 T1-T2 的平均丰度分布)
ref_samples <- sample_info %>% filter(Group=="SC", Time %in% c("DA0","DB3","DC7","DD14","DE21","DF28","DG60"))
ref_abund <- asv_long %>%
  filter(SampleID %in% ref_samples$SampleID) %>%
  group_by(ASV) %>%
  summarise(Ref=mean(Abundance), .groups="drop")

# 计算每个样本与参考状态的欧式距离
calc_dev <- function(sample_id){
  abund <- asv_long %>%
    filter(SampleID==sample_id) %>%
    select(ASV, Abundance)
  df <- left_join(ref_abund, abund, by="ASV") %>%
    replace_na(list(Abundance=0))
  dist(rbind(df$Ref, df$Abundance), method="euclidean")[1]
}

deviation_df <- sample_info %>%
  mutate(Deviation = sapply(SampleID, calc_dev))

# 假设 Treatment_T4 是干扰点
disturb_time <- "DB3"
dev_at_disturb <- deviation_df %>%
  filter(Group=="SH", Time==disturb_time) %>%
  pull(Deviation)

# 绘图
p <- ggplot(deviation_df %>% filter(Group=="SH"),
            aes(x=Time, y=Deviation, group=Group)) +
  geom_line(color="orange", size=1) +
  geom_point(size=2, color="orange") +
  
  # 干扰事件线
  geom_vline(xintercept=disturb_time, linetype="dashed", color="red") +
  annotate("text", x=disturb_time, y=max(deviation_df$Deviation)*0.9,
           label=paste0("干扰事件 (", disturb_time, ")"), color="red", size=3) +
  
  # RS 抗性
  annotate("segment", x=disturb_time, xend="DB3", 
           y=dev_at_disturb, yend=max(deviation_df$Deviation)*0.9,
           arrow=arrow(length=unit(0.2,"cm")), color="black") +
  annotate("text", x="DB3", y=max(deviation_df$Deviation)*0.95,
           label="RS 抗性\n瞬时偏离大小", size=3, hjust=0.5) +
  
  labs(x="时间点", y="与参考状态的偏离度",
       title="群落稳定性指标示意图 (基于 ASV 丰度表)") +
  theme_bw(base_size=12)

print(p)
