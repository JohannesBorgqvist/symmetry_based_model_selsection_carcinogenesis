
library(nonnest2)

setwd("./")

datacsv=read.csv("./Data/Colon_cancer_all_fields.csv")
rates=datacsv[,"AgeAdjusted_Rate"][1:86]


data <- structure(list(R = rates, 
                       t = 1:86), 
                  .Names = c("R", "t"), 
                  class = "data.frame", 
                  row.names = c(NA, -86L))
data

library(drc)


mm.model.nls <- nls(R~A*t^B, data=data, 
                    start = list(A=max(data$R)/2, 
                                 B=4))
summary(mm.model.nls)

mm.model.nls2 <- nls(R~A/(exp(exp(-0.044*(t-tau)))-1), data=data, 
                    start = list(A=2,tau=58))

summary(mm.model.nls2)

vuongtest(mm.model.nls,mm.model.nls2)

# Hmm it seems they are distinguishable. Let's try taking log Risk

data <- structure(list(R = log(rates[14:86]), 
                       t = 14:86), 
                  .Names = c("R", "t"), 
                  class = "data.frame", 
                  row.names = c(NA, -73L))

mm.model.nls <- nls(R~A+B*log(t), data=data, 
                    start = list(A=max(data$R)/2, 
                                 B=4))
summary(mm.model.nls)

mm.model.nls2 <- nls(R~A-log(exp(exp(-0.044*(t-tau)))-1), data=data, 
                     start = list(A=5,tau=58))

summary(mm.model.nls2)

vuongtest(mm.model.nls,mm.model.nls2)

library(ggplot2)
MM.model <- drm(R~t, data=data, fct=MM.2())
mmdf <- data.frame(t=seq(14,max(data$t),length.out=100))
mmdf$R <- predict(mm.model.nls, newdata=mmdf)
mmdf2 <- data.frame(t=seq(14,max(data$t),length.out=100))
mmdf2$R <- predict(mm.model.nls2, newdata=mmdf)
ggplot(data, aes(x = t, y = R)) +
  theme_bw() +
  xlab("Age") +
  ylab("Risk") +
  ggtitle("Fitting both models with PLM in green and IM-II in orange") +
  geom_point(alpha = 0.5) +
  geom_line(data = mmdf, 
            aes(x = t, y = R), 
            colour = "green")+
  geom_line(data = mmdf2, 
            aes(x = t, y = R), 
            colour = "orange")
