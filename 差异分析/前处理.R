#前处理步骤：
#1. 过滤离群样本：根据糖的鉴定数量过滤离群样本（鉴定量过低的样本），可以用boxplot离群值检测法；
#2. 过滤无效糖：过滤掉空值比例大于50%的糖，即如果一个糖在超过50%的样本中都是空值，舍弃这个糖；
#3. 填空值：对空值进行填充，填充每个样本中所有糖最低强度的1/2（也就是说每个样本填的值都不一样）；
#你如果感兴趣可以查一下更高级的填空值方法，比如KNN或者missForest；
#4. 归一化：每个糖的强度除以该样本中所有糖强度的总和。

#所用的包与数据
#tidyverse
#abundance.csv

#1:过滤离群样本，生成abundance_1
NA_count_abundance = apply(abundance,1,function(x) sum(is.na(x)))
min(boxplot.stats(NA_count_abundance)$out)#最小离群值为19
abundance_1 =  abundance |>
    mutate(NA_count = unlist(NA_count_abundance),
           .before = 2) |>
    filter(NA_count < 19)

#2:过滤无效糖，生成abundance_2
abundance_2_1 = abundance_1 |>
    select(!NA_count) |>
    pivot_longer(cols = starts_with("H"),
                 names_to = "glycan_type",
                 values_to = "glycan_value") |>
    pivot_wider(names_from = sample,
                values_from = glycan_value) 
NA_count_abundance_1 = apply(abundance_2,1,function(x) sum(is.na(x)))
abundance_2 = abundance_2_1 |>
    mutate(NA_count = unlist(NA_count_abundance_1),
           .before = 2) |>
    filter(NA_count <= 87) |>
    select(!NA_count) |>
    pivot_longer(cols = starts_with("S"),
                 names_to = "sample",
                 values_to = "glycan_value") |>
    pivot_wider(names_from = glycan_type,
                values_from = glycan_value)

#3:填空值(最低值1/2)，生成abundance_3
abundance_3_1 = abundance_2 |>
    mutate(H_min_value = apply(abundance_2 |> select(!sample),1,min,na.rm = T)/2,
           .before = 2) |> #这里加个H是为了方便转换
    pivot_longer(cols = starts_with("H"),
                 names_to = "glycan_type",
                 values_to = "glycan_value") |>
    pivot_wider(names_from = sample,
                values_from = glycan_value) 
for(i in seq(2,175,1)){
    abundance_3_1[,i][is.na(abundance_3_1[,i])] = min(abundance_3_1[,i],na.rm = T)
}
abundance_3 = abundance_3_1 |>
    pivot_longer(cols = starts_with("S"),
                 names_to = "sample",
                 values_to = "glycan_value") |>
    pivot_wider(names_from = glycan_type,
                values_from = glycan_value) |>
    select(!H_min_value)

#4:归一化处理，生成前处理完毕的abundance_4
abundance_4_1 = abundance_3 |>
    pivot_longer(cols = starts_with("H"),
                 names_to = "glycan_type",
                 values_to = "glycan_value") |>
    count(sample,wt = glycan_value)
abundance_4 = abundance_3
for(sample in seq(1,174,1)){
    for(glycan in seq(2,65,1)){
        abundance_4[sample,glycan] = 
            unlist(abundance_3[sample,glycan])/unlist(abundance_4_1[sample,2])
    }
}

#导出表格
write.csv(abundance_4,file = "abundance_processed.csv")
