## 前面把加载包的代码写上
library(tidyverse)  #

## 读取数据的代码也写上，注意相对路径问题，参考 https://r4ds.hadley.nz/data-import
abundance = read_csv("abundance.csv")  # 注意不是read.csv

## 以上两点保证代码重复性

## 将数据转化为tidy格式，详见 https://tidyr.tidyverse.org/articles/tidy-data.html#defining
abundance_long <- abundance |>
  pivot_longer(-sample, names_to = "glycan", values_to = "value")
abundance_long

#1:过滤离群样本，生成abundance_1
# NA_count_abundance = apply(abundance,1,function(x) sum(is.na(x)))
# min(boxplot.stats(NA_count_abundance)$out)#最小离群值为19
# abundance_1 =  abundance |>
#     mutate(NA_count = unlist(NA_count_abundance),
#            .before = 2) |>
#     filter(NA_count < 19)

## ===== 关于apply函数 =====
## 最好避免使用base的函数，例如apply，尽量使用tidyverse。
## base函数通常更精简，但是可读性较差，不易维护。
## 对数据分析的代码来说，可读性 >> 效率。
## 因为代码大概率是需要反复修改的，未来的自己能看懂才能保证代码的可维护性。

## ===== 关于硬编码 =====
## 19是一个硬编码。未来如果输入数据有变化，这个19可能就不对了。
## 建议避免一切硬编码，尽量用变量来记录中间结果。

## ===== 关于离群样本 =====
## 有两种方式，一种是根据鉴定量来判断，一种是根据一张谱图的质谱总信号强度来判断。
## 这里我使用第二种，第一种的代码类似。
samples_to_keep <- abundance_long |> 
  # 计算每个样本的总信号强度
  summarise(total_abund = sum(value, na.rm = TRUE), .by = sample) |> 
  # 过滤离群样本
  filter(total_abund > quantile(total_abund, 0.25) - 1.5 * IQR(total_abund))  
  # 提取样本名
  pull(sample)
  
abundance_long_1 <- abundance_long |>
  filter(sample %in% samples_to_keep$sample)

#2:过滤无效糖，生成abundance_2
# abundance_2_1 = abundance_1 |>
#     select(!NA_count) |>
#     pivot_longer(cols = starts_with("H"),
#                  names_to = "glycan_type",
#                  values_to = "glycan_value") |>
#     pivot_wider(names_from = sample,
#                 values_from = glycan_value) 
# NA_count_abundance_1 = apply(abundance_2,1,function(x) sum(is.na(x)))
# abundance_2 = abundance_2_1 |>
#     mutate(NA_count = unlist(NA_count_abundance_1),
#            .before = 2) |>
#     filter(NA_count <= 87) |>
#     select(!NA_count) |>
#     pivot_longer(cols = starts_with("S"),
#                  names_to = "sample",
#                  values_to = "glycan_value") |>
#     pivot_wider(names_from = glycan_type,
#                 values_from = glycan_value)

## ===== 关于pivot_longer和pivot_wider =====
## 慢慢习惯在长数据上进行各种分析，避免使用宽数据。

glycans_to_keep <- abundance_long_1 |> 
  # 计算每个糖的缺失值比例
  summarise(na_prop = mean(is.na(value)), .by = glycan) |>
  # 过滤缺失值比例大于0.5的糖
  filter(na_prop <= 0.5) |>
  # 提取糖名
  pull(glycan)


abundance_long_2 <- abundance_long_1 |>
  filter(glycan %in% glycans_to_keep)


#3:填空值(最低值1/2)，生成abundance_3
# abundance_3_1 = abundance_2 |>
#     mutate(H_min_value = apply(abundance_2 |> select(!sample),1,min,na.rm = T)/2,
#            .before = 2) |> #这里加个H是为了方便转换
#     pivot_longer(cols = starts_with("H"),
#                  names_to = "glycan_type",
#                  values_to = "glycan_value") |>
#     pivot_wider(names_from = sample,
#                 values_from = glycan_value) 
# for(i in seq(2,175,1)){
#     abundance_3_1[,i][is.na(abundance_3_1[,i])] = min(abundance_3_1[,i],na.rm = T)
# }
# abundance_3 = abundance_3_1 |>
#     pivot_longer(cols = starts_with("S"),
#                  names_to = "sample",
#                  values_to = "glycan_value") |>
#     pivot_wider(names_from = glycan_type,
#                 values_from = glycan_value) |>
#     select(!H_min_value)

abundance_long_3 <- abundance_long_2 |> 
  # 计算每个样本的最小值
  group_by(sample) |> 
  mutate(min_value = min(value, na.rm = TRUE)) |>
  # 用最小值的一半填充缺失值
  mutate(value = if_else(is.na(value), min_value / 2, value)) |> 
  select(-min_value) |> 
  ungroup()

## ===== group_by + mutate =====
## 注意mutate函数在group_by之后的行为会发生变化，
## 详见 https://dplyr.tidyverse.org/articles/grouping.html#mutate


#4:归一化处理，生成前处理完毕的abundance_4
# abundance_4_1 = abundance_3 |>
#     pivot_longer(cols = starts_with("H"),
#                  names_to = "glycan_type",
#                  values_to = "glycan_value") |>
#     count(sample,wt = glycan_value)
# abundance_4 = abundance_3
# for(sample in seq(1,174,1)){
#     for(glycan in seq(2,65,1)){
#         abundance_4[sample,glycan] = 
#             unlist(abundance_3[sample,glycan])/unlist(abundance_4_1[sample,2])
#     }
# }

abundance_long_4 <- abundance_long_3 |> 
  group_by(sample) |> 
  mutate(value = value / sum(value)) |> 
  ungroup()

## 最后，转化为宽数据储存（宽数据更节省磁盘空间）
abundance_processed <- abundance_long_4 |> 
  pivot_wider(names_from = glycan, values_from = value)

#导出表格
# write.csv(abundance_4,file = "abundance_processed.csv")

## ===== 关于write_csv和write.csv =====
## 建议使用tidyverse的write_csv函数，而不是base的write.csv函数。
## 详见 https://readr.tidyverse.org/reference/write_delim.html

write_csv(abundance_processed, "abundance_processed.csv")