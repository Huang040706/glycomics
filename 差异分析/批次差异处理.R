#1：导入分组信息
#2：取log2
#3：批次差异处理
#4：取2次方，导出

#导入包
library(tidyverse)
library(sva)

#导入文件，这里把批次处理放在了前处理计算相对强度之后
plate = read_csv("plates.csv")
abundance = read_csv("abundance_processed.csv")

#1：导入分组信息
abundance_plate = plate |>
    select("plate","sample") |>
    filter(sample %in% abundance$sample)
combat_plate = data.frame(abundance_plate)
rownames(combat_plate) = abundance_plate$sample
combat_plate$plate = as.character(abundance_plate$plate)

#2：取log2
abundance_log2 = abundance |>
    pivot_longer(-sample,
                 names_to = "glycan",
                 values_to = "value") |>
    mutate(value = log2(value)) |>
    pivot_wider(names_from = "glycan",
                values_from = "value")
abundance_matrix = as.matrix(abundance_log2)
rownames(abundance_matrix) = abundance$sample
abundance_matrix = abundance_matrix[,2:64]
combat_glycan = t(abundance_matrix)
combat_glycan = apply(combat_glycan,
                      c(1,2),
                      as.numeric)

#3：批次差异处理
combat_consequence = ComBat(dat = combat_glycan,
                            batch = combat_plate$plate)

#4：取2次方，导出
result = as_tibble(t(combat_consequence)) |>
    mutate(sample = abundance$sample,
           .before = 0) |>
    pivot_longer(-sample,
                 names_to = "glycan",
                 values_to = "value") |>
    mutate(value = 2^(value)) |>
    pivot_wider(names_from = "glycan",
                values_from = "value")
write_csv(result,
          file = "abundance_combated.csv")




