## 1) Load libraries
library(tidyverse)
library(ggplot2)
library(plantecophys)
library(car)
library(MetBrewer)

setwd("~/Library/CloudStorage/GoogleDrive-kellykerr@ucsb.edu/My Drive/Projects/Walnut Creek AZ Garden/RStudio/Walnut Creek Garden AZ")

## 2)  Load a-ci curves
a.ci <- read.csv("aci-R.csv") 
str(a.ci)
head(a.ci)

## 3) Use "Curve" as a unique identifier for each sample and set it as a separate vector
# Check out each curve for fit
# a.ci$Curve <- paste(a.ci$time, a.ci$population, a.ci$clone, a.ci$bed, a.ci$plant, sep="-")
# a.ci$sample <- paste(a.ci$bed, a.ci$plant)
trees <- unique(a.ci$id)
trees <- sort(trees)
trees

for (i in 1:length(trees)){
  test <- a.ci[a.ci$id == trees[i],]
  t <- fitaci(test, fitmethod="default", Tcorrect=FALSE)
  plot(t, main=trees[i])
}

# Remove bad curves and rerun
a.ci2 <- subset(a.ci,  id !="H8-3" & id !="N 10-1") #remove branches that are not good
trees2 <- unique(a.ci2$id)
trees2 <- sort(trees2)
trees2

for (i in 1:length(trees2)){
  test <- a.ci2[a.ci2$id == trees2[i],]
  t <- fitaci(test, fitmethod="default", Tcorrect=FALSE)
  plot(t, main=trees2[i])
}

# Rerun bad curves with bilinear fit
a.ci.bilinear <- a.ci %>%
  filter(id %in% c("H8-3", "N 10-1"))
trees3 <- unique(a.ci.bilinear$id)
trees3 <- sort(trees3)
trees3

for (i in 1:length(trees3)){
  test <- a.ci.bilinear[a.ci.bilinear$id == trees3[i],]
  t <- fitaci(test, fitmethod="bilinear", Tcorrect=FALSE)
  plot(t, main=trees3[i])
}


## 4) Store A/Ci curve as dataframe
vcmax.default <- matrix(data=NA, nrow=length(trees2), ncol=6)

for (i in 1:length(trees2)) {
  test <- a.ci2[a.ci2$id == trees2[i],]
  t <- fitaci(test, fitmethod="default", Tcorrect=FALSE)
  vcmax.default[i,1]<-t$pars[1,1]
  vcmax.default[i,2]<-t$pars[1,2]
  vcmax.default[i,3]<-t$pars[2,1]
  vcmax.default[i,4]<-t$pars[2,2]
  vcmax.default[i,5]<-t$pars[3,1]
  vcmax.default[i,6]<-t$pars[3,2]
}

vcmax.default <- as.data.frame(vcmax.default)
vcmax.default <- cbind(trees2, vcmax.default)
colnames(vcmax.default) <- c("id", "Vcmax", "vcmax.error", "Jmax", "Jmax.error",
                     "Rd", "Rd.error")


vcmax.bilinear <- matrix(data=NA, nrow=length(trees3), ncol=6)

for (i in 1:length(trees3)) {
  test <- a.ci.bilinear[a.ci.bilinear$id == trees3[i],]
  t <- fitaci(test, fitmethod="bilinear", Tcorrect=FALSE)
  vcmax.bilinear[i,1]<-t$pars[1,1]
  vcmax.bilinear[i,2]<-t$pars[1,2]
  vcmax.bilinear[i,3]<-t$pars[2,1]
  vcmax.bilinear[i,4]<-t$pars[2,2]
  vcmax.bilinear[i,5]<-t$pars[3,1]
  vcmax.bilinear[i,6]<-t$pars[3,2]
}

vcmax.bilinear <- as.data.frame(vcmax.bilinear)
vcmax.bilinear <- cbind(trees3, vcmax.bilinear)
colnames(vcmax.bilinear) <- c("id", "Vcmax", "vcmax.error", "Jmax", "Jmax.error",
                     "Rd", "Rd.error")

vcmax <- bind_rows(vcmax.default, vcmax.bilinear)

# Create crosstype column
vcmax <- vcmax %>%
  mutate(crosstype = str_extract(id, "[A-Za-z]+"))

## 5) Summarize and graph Vcmax
vcmax.sum <- vcmax %>% 
  group_by(crosstype) %>%
  dplyr::summarise(mean.Vcmax = mean(Vcmax, na.rm = TRUE),
                   median.Vcmax = median(Vcmax, na.rm=TRUE),
                   sd.Vcmax = sd(Vcmax, na.rm = TRUE),
                   n.Vcmax = length(Vcmax[!is.na(Vcmax)])) %>%
  mutate(se.Vcmax = sd.Vcmax / sqrt(n.Vcmax),
         lower.ci.Vcmax = mean.Vcmax - qt(1 - (0.05 / 2), n.Vcmax - 1) * se.Vcmax,
         upper.ci.Vcmax = mean.Vcmax + qt(1 - (0.05 / 2), n.Vcmax - 1) * se.Vcmax)

vcmax.sum$crosstype <- factor(vcmax.sum$crosstype,levels = c("F","BCF","H","BCN","N"))


ggplot(data=vcmax.sum, aes(x=crosstype, y=mean.Vcmax, color=crosstype)) +
  geom_pointrange(aes(ymin = mean.Vcmax-se.Vcmax, ymax = mean.Vcmax+se.Vcmax), size=2,
                  position=position_dodge(width=2)) +
  theme_bw() 

ggplot(data=vcmax.sum, aes(x=crosstype, y=mean.Vcmax, color=crosstype)) +
  geom_point(size=10) +
  geom_errorbar(aes(ymin=mean.Vcmax-se.Vcmax, ymax=mean.Vcmax+se.Vcmax), size=1.5, width=0.5) +
  theme(plot.margin=unit(c(1,0.5,1,1),"cm")) +
  theme(plot.title = element_text(margin=margin(0,0,30,0))) +
  labs(x="", y=expression(italic(V[cmax])~(mu*mol~m^-2~s^-1))) +
  theme(axis.text=element_text(size=26, color="black"),
        axis.title=element_text(size=28), 
        panel.background = element_rect(fill = "white",
                                        colour = "white"),
        axis.line = element_line(color = 'black'),
        legend.text=element_text(size=22),
        legend.title=element_text(size=24)) +
  scale_color_manual(values=met.brewer("Johnson", 5)) +
  guides(color=guide_legend(title="Crosstype")) 
