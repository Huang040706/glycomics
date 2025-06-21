#任务一：四组数据的差异分析（HCC/CHB/HC/LC）
#1：导入疾病分组信息
#2：筛选四组间ANOVA差异分析后p显著的糖
#3：热图绘制
#3.1：z-score归一化
#3.2：导入热图分组信息
#3.3：使用Complex heatmap绘制热图

#导入所用包
library(tidyverse)
library(rstatix)
library(ComplexHeatmap)
library(circlize)

#导入所用文件，这里用的是批次差异处理完后的文件
abundance_combated = read_csv("abundance_combated.csv")
group = read_csv("groups.csv")
group_sample = read_csv("group_sample.csv") |> arrange(sample)

#1：导入疾病分组信息
merged_abundance = merge(group,abundance_combated,by = "sample")

#2：筛选四组间ANOVA差异分析后p显著的糖，生成一个p-value list
merged_abundance_long = merged_abundance |>
    pivot_longer(starts_with("H"),
    names_to = "glycan",
    values_to = "value") 
merged_abundance_long_log = merged_abundance_long |>#用log10处理偏态数据
    mutate(value = log10(value))
modelA <- as.formula(str_c("value"," ~ ","group"))
abundance_p_value = merged_abundance_long_log |> 
    group_by(glycan) |>
    nest() |>
    mutate(model = map(data,~aov(modelA,data = .))) |>
    mutate(res = map(model,~tidy(.))) |>
    unnest(res) |>
    filter(term == "group") 
p_value_adjusted_list = p.adjust(abundance_p_value$p.value,method = "BH") #BH法校正
abundance_p_value = cbind(abundance_p_value,p.value_adjusted = p_value_adjusted_list)
glycans_to_keep = abundance_p_value |>
    filter(p.value_adjusted < 0.05) |>
    pull(glycan)

#3：热图绘制
#3.1：z-score归一化
z_score_normarlization = function(df){
    data = as_tibble(data.frame(df[[1]]))$value
    data_mean = mean(data)
    data_std_error = sd(data)
    new_df = as_tibble(data.frame((df[[1]]))) |>
        mutate(z_score = (value - data_mean) / data_std_error)
    return(new_df)
}
merged_abundance_long_log_z_score = merged_abundance_long_log |>
    filter(glycan %in% glycans_to_keep) |>
    group_by(glycan) |>
    nest() |>
    mutate(z_score = map(data,~z_score_normarlization(data))) |>
    unnest(z_score) |>
    select(!data) |>
    ungroup()

#3.2：导入热图分组信息
heatmap_df = merged_abundance_long_log_z_score |>
    select(!value) |>
    select(!group) |>
    pivot_wider(names_from = "sample",
                values_from = "z_score")
heatmap_matrix = as.matrix(heatmap_df |> select(!glycan))
rownames(heatmap_matrix) = heatmap_df$glycan
colnames(heatmap_matrix) = group_sample$group_sample
heatmap_matrix_transformed = t(heatmap_matrix)
heatmap_matrix_transformed_df = as.data.frame(heatmap_matrix_transformed)
heatmap_matrix_transformed_df <- cbind(RowName = rownames(heatmap_matrix_transformed_df), heatmap_matrix_transformed_df)
heatmap_matrix_transformed_df_sorted <- heatmap_matrix_transformed_df[order(heatmap_matrix_transformed_df$RowName),]
heatmap_matrix = t(as.matrix(heatmap_matrix_transformed_df_sorted[, -1]))

#3.3：使用Complex heatmap绘制热图
ht_opt$TITLE_PADDING = unit(c(5, 5), "points")
draw(Heatmap(heatmap_matrix,
             col = colorRamp2(c(-2, 0, 2),c("blue","white","red")),#限定z值范围防止离群值影响
             name = "Z-score",
             row_title = "Glycan",
             show_column_names = FALSE,
             row_km = 5,
             column_split = factor(rep(c("HC","CHB","LC","HCC"),each = 50),
                                   levels = c("HC","CHB","LC","HCC")),#这里是硬编码
             column_title_gp = gpar(fill = c("#FFEBCD", "#FFDEAD", "#FFD700","#DAA520"),
                                    border = "white",
                                    fontsize = 10),
             cluster_column_slices = FALSE,
             cluster_row_slices = FALSE,
             cluster_columns = FALSE,
             row_names_side = "right",
             row_names_rot = 45,
             row_names_gp = gpar(fontsize = 7),
             row_title_gp = gpar(fontsize = 10)))

