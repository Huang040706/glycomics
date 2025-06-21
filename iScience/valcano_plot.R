#所用的包
#library(tidyverse)
#library(ggrepel)

#提取相关列用于后续分析
valcano_sheet = raw |>
    select(!xy_filename:ID_Sample) |>
    select(!Severity:Age_groups) |>
    filter(Group_letter %in% c("Visceral Leishmaniasis",
                               "Asymptomatic",
                               "Healthy from Endemic area")) 

#求三组每种糖的平均值
part_VL = valcano_sheet |>
    filter(Group_letter == "Visceral Leishmaniasis") |>
    pivot_longer(
        cols = starts_with("H"),
        names_to = "glycan_type",
        values_to = "value"
    ) |>
    group_by(glycan_type) |>
    summarise(average_value = mean(value))

part_ASYMP = valcano_sheet |>
    filter(Group_letter == "Asymptomatic") |>
    pivot_longer(
        cols = starts_with("H"),
        names_to = "glycan_type",
        values_to = "value"
    ) |>
    group_by(glycan_type) |>
    summarise(average_value = mean(value))

part_EC = valcano_sheet |>
    filter(Group_letter == "Healthy from Endemic area") |>
    pivot_longer(
        cols = starts_with("H"),
        names_to = "glycan_type",
        values_to = "value"
    ) |>
    group_by(glycan_type) |>
    summarise(average_value = mean(value))

#生成平均值的log2FC值作为x轴
FC_sheet = part_VL |>
    mutate(average_value_VL = average_value) |>
    mutate(average_value_ASYMP = part_ASYMP$average_value) |>
    mutate(average_value_EC = part_EC$average_value) |>
    select(!average_value) |>
    mutate(Log2_FC_VL_ASYMP = log2(average_value_VL/average_value_ASYMP)) |>
    mutate(Log2_FC_VL_EC = log2(average_value_VL/average_value_EC))

#Wilcoxon分析函数
wilcoxon_test_ASYMP = function(x){
    sheet1 = valcano_sheet |>
        filter(Group_letter == "Visceral Leishmaniasis") |>
        select(x)
    sheet2 = valcano_sheet |>
        filter(Group_letter == "Asymptomatic") |>
        select(x)
    wilcoxon_test_result = wilcox.test(x = unlist(sheet1),y = unlist(sheet2))
    return(wilcoxon_test_result$p.value)
}

wilcoxon_test_EC = function(x){
    sheet1 = valcano_sheet |>
        filter(Group_letter == "Visceral Leishmaniasis") |>
        select(x)
    sheet2 = valcano_sheet |>
        filter(Group_letter == "Healthy from Endemic area") |>
        select(x)
    wilcoxon_test_result = wilcox.test(x = unlist(sheet1),y = unlist(sheet2))
    return(wilcoxon_test_result$p.value)
}

#将Wilcoxon分析结果汇总
glycan_list = unlist(FC_sheet$glycan_type)
wilcoxon_result_ASYMP = list()
wilcoxon_result_EC = list()
for (name in glycan_list) {
    wilcoxon_result_ASYMP = c(wilcoxon_result_ASYMP,
                              wilcoxon_test_ASYMP(name))
    wilcoxon_result_EC = c(wilcoxon_result_EC,
                           wilcoxon_test_EC(name))
}
p_sheet = tibble(
    glycan_type = glycan_list,
    p_value_ASYMP = wilcoxon_result_ASYMP,
    p_value_EC = wilcoxon_result_EC
)

#创建含-log10p与log2FC的总表，用于绘图
ultimate_sheet = FC_sheet |>
    mutate(P_ASYMP = unlist(p_sheet$p_value_ASYMP),
           P_EC = unlist(p_sheet$p_value_EC)) |>
    mutate(Log10_P_ASYMP = -log10(P_ASYMP),
           Log10_P_EC = -log10(P_EC)) |>
    select(glycan_type,
           Log2_FC_VL_ASYMP,
           Log2_FC_VL_EC,
           Log10_P_ASYMP,
           Log10_P_EC)

#通过FC值与显著性分组
split_group_ASYMP = function(){
    LIST = list()
    for (i in 1:73) {
        if(ultimate_sheet$Log2_FC_VL_ASYMP[i] > 0.3 &&
           ultimate_sheet$Log10_P_ASYMP[i] > -log10(5.0e-5)){
            LIST = c(LIST,"up")#显著上调
        }else if(ultimate_sheet$Log2_FC_VL_ASYMP[i] < -0.3 &&
                 ultimate_sheet$Log10_P_ASYMP[i] > -log10(5.0e-5)){
            LIST = c(LIST,"down")#显著下调
        }else{
            LIST = c(LIST,"ns")#无差异    
        }
    }
    return(LIST)
}

split_group_EC = function(){
    LIST = list()
    for (i in 1:73) {
        if(ultimate_sheet$Log2_FC_VL_EC[i] > 0.3 &&
           ultimate_sheet$Log10_P_EC[i] > -log10(5.0e-5)){
            LIST = c(LIST,"up")#显著上调
        }else if(ultimate_sheet$Log2_FC_VL_EC[i] < -0.3 &&
                 ultimate_sheet$Log10_P_EC[i] > -log10(5.0e-5)){
            LIST = c(LIST,"down")#显著下调
        }else{
            LIST = c(LIST,"ns")#无差异    
        }
    }
    return(LIST)
}
    
ultimate_sheet = ultimate_sheet |>
    mutate(Group_ASYMP = split_group_ASYMP(),
           Group_EC = split_group_EC()) 
    
#画VL/EC火山图
ultimate_sheet |>
    ggplot(aes(x = Log2_FC_VL_EC,
               y = Log10_P_EC)) +
    geom_point(data = ultimate_sheet |> filter(Group_EC == "up"),
               color = "red",alpha = 0.5) +
    geom_point(data = ultimate_sheet |> filter(Group_EC == "down"),
               color = "green",alpha = 0.5) +
    geom_point(data = ultimate_sheet |> filter(Group_EC == "ns"),
               color = "#696969",alpha = 0.5) +
    labs(title = "B",
         x = "Log2 fold change(VL/Control)",
         y = "-Log10P") +
    xlim(-1.2,2) +
    ylim(0,40) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.grid.major = element_line(color = "#DCDCDC"),
          panel.grid.minor = element_line(color = "#DCDCDC")) +
    theme(axis.line = element_line(color = "black")) +
    geom_hline(yintercept = -log10(5.0e-5),
               linetype = 2,
               size = 0.5,
               alpha = 0.5) +
    geom_vline(xintercept = c(-0.3,0.3),
               linetype = 2,
               size = 0.5,
               alpha = 0.5) +
    geom_text_repel(data = ultimate_sheet |> filter(Group_EC != "ns"),
                    aes(label = glycan_type),
                    arrow = arrow(length = unit(0.05,"inches"),ends = "first"),
                    min.segment.length = 0) +
    coord_fixed(0.045)
    
