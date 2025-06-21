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
raw |>
    filter(Group_letter %in% c("Healthy from Endemic area",
                               "Asymptomatic",
                               "Visceral Leishmaniasis")) |> 
    select(Sex,Group_letter,H3N4F1E0L0) |>
    filter(Sex %in% c("M","F")) |>
    mutate(wrap_seq = wrap_sequence(Group_letter)) |>
    ggplot(aes(x = Sex,y = H3N4F1E0L0)) +
    geom_boxplot(outlier.size = 0) +
    geom_point(position = position_jitter(width = 0.1),aes(shape = Sex,color = Sex)) +
    facet_wrap(~wrap_seq,labeller = labeller(wrap_seq = c(
        "1" = "EC",
        "2" = "ASYMP",
        "3" = "VL"))) +
    lims(y = c(0.00,0.23)) +
    labs(title = "H3N4F1E0L0",x = NULL,y = "Relative area")
