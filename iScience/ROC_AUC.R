#载入使用的包
#tidyverse
#plotROC
#ggrepel
#magick

#载入使用的表格
#raw
#glycan_type

#EC-VL与ASYMP-VL的二分类函数，可同时用于两者的分别检验
#为了解决AUC小于0.5的问题，新增了一个反转函数
classification = function(x){
    list = c()
    for (i in x) {
        if(i == "Visceral Leishmaniasis"){
            list = c(list,1)
        }else{
            list = c(list,0)
        }    
    }
    return(unlist(list))
}
classification_adverse = function(x){
    list = c()
    for (i in x) {
        if(i == "Visceral Leishmaniasis"){
            list = c(list,0)
        }else{
            list = c(list,1)
        }    
    }
    return(unlist(list))
}


#EC-VL的ROC分析
EC_VL_ROC = raw |>
    filter(Group_letter %in% c("Visceral Leishmaniasis",
                               "Healthy from Endemic area")) |>
    select(!xy_filename:ID_Sample) |>
    select(!Severity:SampleType) |>
    select(!Group:Age_groups) 

#以H3N3F0E0L0为例，画ROC曲线
ggroc_EC = EC_VL_ROC |>
    mutate(P_N = classification(Group_letter)) |>
    ggplot(aes(d = P_N,m = H4N3F0E0L0,color = SampleSet)) +
    scale_x_continuous(guide = guide_axis(minor.ticks = TRUE),sec.axis = dup_axis()) +
    scale_y_continuous(guide = guide_axis(minor.ticks = TRUE),sec.axis = dup_axis()) +
    geom_roc(n.cuts = 0) +
    labs(title = "A    Healthy from Endemic area - Visceral Leishmaniasis",
         subtitle = "H4N3F0L0E0",
         x = "False positive rate",
         y = "True positive rate") +
    theme(panel.background = element_blank()) +
    theme(axis.line = element_line(color = "black")) +
    theme(axis.title.x.top = element_blank()) +
    theme(axis.title.y.right = element_blank()) +
    theme(axis.ticks.x.top = element_blank()) +
    theme(axis.ticks.y.right = element_blank()) +
    theme(axis.text.x.top = element_blank()) +
    theme(axis.text.y.right = element_blank()) +
    coord_fixed(1) +
    geom_abline(slope = 1,intercept = 0,alpha = 0.7,color = "grey")
ggroc_EC

#计算并显示AUC值（EC-VL组）
auc_list = calc_auc(ggroc_EC)
ggroc_EC = ggroc_EC +
    scale_color_discrete(labels = c(paste("Discovery  AUC=",round(auc_list$AUC[1],3),sep = ""),
                                    paste("Validation  AUC=",round(auc_list$AUC[2],3),sep = ""))) +
    theme(legend.title = element_blank(),
          legend.position = c(0.05,0.95),
          legend.background = element_blank(),
          legend.justification = c("left","top"))

#遍历所有糖类型，统计两种分类下的AUC值
#这里设计成函数主要为了防止多次运行程序时ROC列表反复叠加
EC_VL_ROC_FUNCTION = function(){
    AUC_LIST_glycan_type = c()
    AUC_LIST_discovery_auc = c()
    AUC_LIST_validation_auc = c()
    for(glycan in glycan_type$glycan_type){
        ggroc_EC_all = NULL
        ggroc_EC_all = EC_VL_ROC |>
            mutate(P_N = classification(Group_letter)) |>
            mutate(P_N_adverse = classification_adverse(Group_letter)) |>
            ggplot(aes(d = P_N,m = get(glycan),color = SampleSet)) +
            geom_roc(n.cuts = 0)
        auc_list_all = calc_auc(ggroc_EC_all)
        auc_discovery = round(auc_list_all$AUC[1],3)
        auc_validation = round(auc_list_all$AUC[2],3)
        auc_discovery = max(auc_discovery,1 - auc_discovery)
        auc_validation = max(auc_validation,1 - auc_validation)
        AUC_LIST_glycan_type = c(AUC_LIST_glycan_type,glycan)
        AUC_LIST_discovery_auc = c(AUC_LIST_discovery_auc,auc_discovery)
        AUC_LIST_validation_auc = c(AUC_LIST_validation_auc,auc_validation)
    }
    return(tibble(AUC_LIST_glycan_type,
           AUC_LIST_discovery_auc,
           AUC_LIST_validation_auc))
}
AUC_tibble = EC_VL_ROC_FUNCTION()
AUC_LIST_glycan_type = c()
AUC_LIST_discovery_auc = c()
AUC_LIST_validation_auc = c()
EC_VL_ROC = EC_VL_ROC |>
    mutate(P_N = classification(Group_letter))
for(glycan in glycan_type$glycan_type){
    ggroc_EC_all = NULL
    ggroc_EC_all = EC_VL_ROC |>
        ggplot(mapping = aes(d = P_N,m = get(glycan),color = SampleSet)) +
        geom_roc(n.cuts = 0)
    auc_list_all = calc_auc(ggroc_EC_all)
    auc_discovery = round(auc_list_all$AUC[1],3)
    auc_validation = round(auc_list_all$AUC[2],3)
    AUC_LIST_glycan_type = c(AUC_LIST_glycan_type,glycan)
    AUC_LIST_discovery_auc = c(AUC_LIST_discovery_auc,auc_discovery)
    AUC_LIST_validation_auc = c(AUC_LIST_validation_auc,auc_validation)
}
AUC_tibble = tibble(AUC_LIST_glycan_type,
                    AUC_LIST_discovery_auc,
                    AUC_LIST_validation_auc
                    )

#利用该AUC表做可视化-柱形散点图
AUC_tibble_filter = AUC_tibble |>
    pivot_longer(
        cols = ends_with("auc"),
        names_to = "Group",
        values_to = "AUC_value"
    ) 
    
AUC_tibble_filter |>
    ggplot(aes(x = Group,y = AUC_value,color = AUC_value)) +
    geom_jitter(width = 0.1) +
    scale_colour_gradient(low = "#7B68EE", high = "red") +
    scale_x_discrete(labels = c("Discovery","Validation")) +
    labs(y = "AUC value",
         title = "AUC values between EC and VL",
         subtitle = "a") +
    theme(panel.background = element_rect(fill = "#FFF0F5"),
          legend.title = element_blank()) +
    coord_fixed(2)

#可视化-散点图
AUC_tibble |>
    ggplot(aes(x = AUC_LIST_discovery_auc,
               y = AUC_LIST_validation_auc,
               fill = "1")) +
    geom_point(color = "#7B68EE") +
    scale_x_continuous(limits = c(0.5, 1)) +
    scale_y_continuous(limits = c(0.5, 1)) +
    coord_fixed(1) +
    labs(x = "Discovery",
         y = "Validation",
         title = " ",
         subtitle = "b") +
    theme(panel.background = element_rect(fill = "#FFF0F5")) +
    geom_text_repel(data = AUC_tibble |> 
                        filter(AUC_LIST_discovery_auc > 0.7 & AUC_LIST_validation_auc > 0.7),
                    aes(label = AUC_LIST_glycan_type),
                    arrow = arrow(length = unit(0.05,"inches"),ends = "first"),
                    min.segment.length = 0,
                    max.overlaps = 100) +
    theme(legend.position = "none")

p1 = image_read("AUC分布图-柱状.png")
p2 = image_read("AUC分布图-散点.png")

pp1 = image_ggplot(p1,interpolate = F)
pp2 = image_ggplot(p2,interpolate = F)

ggarrange(pp1,pp2)
