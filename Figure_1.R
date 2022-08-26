#Authored by: Carlie Perricone

ggplot(len, aes(x=Age, y=mm, color=Salinity, shape=Salinity)) +
  geom_point() +
  geom_smooth(method="loess", size=.7) +
  scale_color_manual(name = "Salinity (ppt)", values=c("#CC79A7", "#56B4E9", "#E69F00"), labels = c("10", "20", "30")) +
  scale_shape_discrete(name = "Salinity (ppt)", labels = c("10", "20", "30")) +
  theme_bw() +
  xlab("Age (DPH)") +
  ylab("Total Length (mm)") +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(3,6, 9, 12, 15, 18, 20, 24)) +
  annotate("text", x=c(6,24), y=c(3.7, 13), label="*", size=12)