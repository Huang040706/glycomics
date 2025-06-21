p1 = image_read("valcano_plot_A.png")
p2 = image_read("valcano_plot_B.png")

pp1 = image_ggplot(p1,interpolate = F)
pp2 = image_ggplot(p2,interpolate = F)

ggarrange(pp1,pp2,ncol = 1)
