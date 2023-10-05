

mod=gam(number~s(biog,k=3,bs='cr') + s(depth,k=3,bs='cr') + s(tpi,k=3,bs='cr'), family=tw,data=dat.species)

## This example is predicting for depth
testdata <- expand.grid(depth=seq(min(dat$depth),max(dat$depth),length.out = 20),
                        biog=mean(mod$model$biog),
                        tpi=mean(mod$model$tpi),
                        method = (mod$model$method)) %>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)

predicts.species.depth = testdata%>%
  data.frame(fits)%>%
  group_by(depth)%>% #only change here
  summarise(number=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()

## You would then go back and change this each time for the different variables 

## This is how you would then plot 
# depth ----
ggmod.species.depth<- ggplot() +
  ylab("")+
  xlab("Depth")+
  geom_point(data=dat.species,aes(x=depth,y=number),  alpha=0.2, size=1,show.legend=F)+
  geom_line(data=predicts.species.depth,aes(x=depth,y=number),alpha=0.5)+
  geom_line(data=predicts.species.depth,aes(x=depth,y=number - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.species.depth,aes(x=depth,y=number + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  Theme1
ggmod.species.depth



