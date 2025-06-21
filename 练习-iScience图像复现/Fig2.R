#iScience-Total serum N-glycans mark visceral leishmaniasis in human infections with Leishmania infantum
#Fig.2复现

#需要导入的包：
#library(tidyverse)
#library(rstatix)
#library(ggpubr)

#创建函数，方便对不同组别病例进行排序
wrap_sequence = function(x){
    wrap_list = c()
    for(i in x){
        if(i == "Visceral Leishmaniasis"){
            wrap_list = c(wrap_list,"3")
        }else if(i == "Healthy from Endemic area"){
            wrap_list = c(wrap_list,"1")
        }else if(i == "Asymptomatic"){
            wrap_list = c(wrap_list,"2")
        }
    }
    return(wrap_list)
}

#raw data预处理，保留性别、组别、组学数据，产生new data可以用于画图、差异分析
new_data = raw |>
    filter(Group_letter %in% c("Healthy from Endemic area",
                               "Asymptomatic",
                               "Visceral Leishmaniasis")) |> 
    select(Sex,Group_letter,H7N6F0E2L2) |>
    filter(Sex %in% c("M","F")) |>
    mutate(wrap_seq = wrap_sequence(Group_letter))

#确定高度，控制差异线位置、图像长度
relative_height = max(new_data$H7N6F0E2L2)*1.25
signficance_height_list = c(relative_height*0.97,
                            relative_height*0.97,
                            relative_height*0.9)
    
#用new data做差异分析(Kruskal-Wallis rank sum test & post-hoc Donn's test)
new_data$Group_letter = factor(new_data$Group_letter)
H_test_result = kruskal.test(H7N6F0E2L2 ~ Group_letter,data = new_data)
D_test_result = new_data |>
    dunn_test(H7N6F0E2L2 ~ Group_letter) |>
    add_significance("p",
                     cutpoints = c(0,1.0e-7,1.0e-6,1.0e-5,5.0e-5,1),
                     symbols = c("****","***","**","*","ns")) |>
    add_xy_position()
D_test_result$y.position = signficance_height_list
D_test_result$xmin = c(2,2,1)
D_test_result$xmax = c(1,3,3)

#用new data画箱线图+散点图
new_data |>
    ggplot(aes(x = wrap_seq,y = H7N6F0E2L2)) +
    geom_boxplot(outlier.shape = NA,aes(shape = Sex)) +
    geom_point(aes(group = Sex,shape = Sex,colour = Sex),
               position = position_jitterdodge(jitter.width = 0.1)) +
    lims(y = c(0,relative_height)) +
    scale_x_discrete(labels =c ("EC", "ASYMP", "VL")) +
    theme(axis.text.x = element_text(angle = 45,size = 12,hjust = 1,vjust = 1)) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(panel.background = element_rect(fill = "#B0E0E6")) +
    labs(title = "H7N6F0E2L2",
         subtitle = paste("p =",sprintf("%.2e",H_test_result$p.value)),
         tag = "F",
         x = NULL,
         y = "Relative area") +
    theme(plot.title = element_text(size = 15,hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 1)) +
    theme(plot.tag.position ="topleft") +
    theme(plot.tag = element_text(size = 15)) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_text(size = 15)) +
    theme(axis.ticks = element_blank()) +
    scale_color_manual(values = c("F" = "#00BFFF", "M" = "#4682B4")) +

    #添加显著性符号
    stat_pvalue_manual(D_test_result,
                       label = "p.signif",
                       hide.ns = T,
                       tip.length = 0.02,
                       label.size = 3) +
    theme(aspect.ratio = 1.42)
