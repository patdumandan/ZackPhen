dry_datA$doys_sqs=(dry_datA$DOYs^2)

require(lme4)
s1=glmer(cbind(tot_flwr, tot_all)~ DOYs+ doy_sq + years + # base linear terms
           DOYs*years + # timing and long-term trend interaction
           doy_sq * years+ (1 |Plot),
         family = "binomial", data=dry_datA)

summary(s1)

s1_preds=as.vector(predict(s1, type="response"))%>%cbind(dry_datA)
colnames(s1_preds)[1]="preds"

ggplot(s1_preds, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dryas (flowers)")+scale_color_viridis_d()+facet_wrap(~Plot)
#  geom_point(aes(x=DOY, y=prop_flwr, col=as.factor(year)))+
#scale_y_log10()


yrs=as.data.frame(rep(1996:2023, 6))
colnames(yrs)=c("year")

coef_names=as.data.frame(coef(s1)$Plot)%>%
tibble::rownames_to_column()%>%
  rename("Plot"="rowname", "Intercept"="(Intercept)", "DOY_yr"="DOYs:years", "doysq_yr"="doy_sq:years")%>%
relocate(c(DOY_yr, years), .after=DOYs)%>%
  mutate(beta_lin=rowSums(.[3:4]), beta_poly=rowSums(.[6:7]))%>%
select(Plot, Intercept, beta_lin, beta_poly)%>%
cross_join(yrs)%>%
  mutate(numerator=beta_lin*year, denominator=(2*(beta_poly*year)),
                                               x_val=numerator/denominator)

#summary table
coef_tab=as.data.frame(summary(s1)$coefficients)%>%
  tibble::rownames_to_column()
