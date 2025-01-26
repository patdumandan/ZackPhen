col_plot2=ggplot(col_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Collembola")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

aca_plot2=ggplot(aca_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Acari")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

bom_plot2=ggplot(bom_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Bombus")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

ich_plot2=ggplot(ich_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Ichneumonidae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

lyc_plot2=ggplot(lyc_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Lycosidae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

sci_plot2=ggplot(sci_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Sciaridae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

pho_plot2=ggplot(pho_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Phoridae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

nym_plot2=ggplot(nym_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Nymphalidae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

mus_plot2=ggplot(mus_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Muscidae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

lin_plot2=ggplot(lin_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Linyphiidae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

cul_plot2=ggplot(cul_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Culicidae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

coc_plot2=ggplot(coc_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Coccoidea")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

chi_plot2=ggplot(chi_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Chironomidae")+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=TotalCatch1, col=as.factor(Year)))+
  facet_wrap(~Plot.ID)+
  scale_y_log10()+
  ylab("log(abundances)")

ggarrange(bom_plot2, col_plot2, ich_plot2, aca_plot2,
          lyc_plot2, sci_plot2,
          common.legend = T, legend="right")

ggarrange(pho_plot2, nym_plot2,
          mus_plot2, lin_plot2, cul_plot2, coc_plot2, chi_plot2,
          common.legend = T, legend="right")
