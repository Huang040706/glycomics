#导入包
library("tidyverse")

#导入文件
abundance_info = read_csv(file = "abundance_combated.csv")
clinical_info = read_csv(file = "clinical.csv")
group_info = read_csv(file = "groups.csv")

#产生总的数据
combined_info = full_join(arrange(group_info,sample),
                      arrange(clinical_info,sample),
                      by = "sample") 
combined_info = full_join(combined_info,abundance_info)

#练习：在HCC组中探究ALBI（白蛋白-胆红素评分，重要的肝功能指标）与糖之间的关系
HCC_ALBI = combined_info |>
    filter(group == "HCC") |>
    select(!group : AAR)
#先尝试分析ALBI和H3N3间的关系。
cor.test(HCC_ALBI$ALBI_score,
         HCC_ALBI$H3N3,
         alternative = "two.side",
         method = "spearman")$p.value




