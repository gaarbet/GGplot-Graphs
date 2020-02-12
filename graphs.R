df <- read.csv("SL_results.csv")
library(gridExtra)
library(tidyverse)
library(grid)
#dmaxp and dminp are for screen
#dmax and dmin are for non-screen
#thickness of points
#plotly for the other graphs to change timepoint? Seems useful to figure out, if possible

#Graph 1
ggplot(data = df, aes (x = factor(dminp,c("-3 to -1", "-1 to 1")), y = sum_diff,col = factor(dmin), shape = factor(dmin), group = factor(dmin)))+
 xlab("Range of Difficulty Parameters for Screen Items") + ylab("Bias Due to Skip-Logic") +
 stat_summary(geom = "point",fun.y = "mean",size = 8, aes(shape = factor(dmin), stroke = 1))+
 stat_summary(geom = "point",fun.y = "mean",size = 8, aes(col = factor(dmin), stroke = 1))+
 scale_shape_manual(name = "Range of Difficulty
Parameters for 
Non-Screen Items",values = c(8,5)) + theme_bw() + 
 scale_colour_manual(name = "Range of Difficulty
Parameters for 
Non-Screen Items",values = c("orange","purple"))+
 theme(panel.grid.major.x  =   element_blank(),
 panel.grid.major =   element_line(colour = "black",size=0.12))+ theme(legend.justification = "top") +
  coord_cartesian(ylim=c(-0.01,0.06))+ scale_y_continuous(breaks = c(-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06))+
  geom_hline(yintercept=0.00)
grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 2))


                       
#################graph 2
ggplot(data = df, aes (x = factor(amin,c("0.3 to 1.5", "1.5 to 3.5")), y = sum_diff, col = factor(probes),shape = factor(probes),group = factor(probes)))+
 xlab("Range of Discrimination Parameters") + ylab("Bias Due to Skip-Logic") +
 stat_summary(geom = "point",fun.y = "mean",size = 8, aes(shape = factor(probes)))+
 stat_summary(geom = "point",fun.y = "mean",size = 8, aes(col = factor(probes)))+
# scale_y_continuous(breaks=scales::pretty_breaks(n = 5))+
  scale_shape_manual(name = "Number
of Screen
Items",values = c(8,6,1)) + theme_bw() + 
  scale_colour_manual(name = "Number
of Screen
Items",values = c("blue","green","red")) + theme_bw()+
 theme(panel.grid.major.x  =   element_blank(),
 panel.grid.major =   element_line(colour = "black",size=0.12))+ theme(legend.justification = "top") +
  coord_cartesian(ylim=c(-0.01,0.06))+ scale_y_continuous(breaks = c(-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06))+
  geom_hline(yintercept=0.00)
  

################# side orientation for graph 3
p <- ggplot(data = df, aes (x = factor(dminp,c("-3 to -1", "-1 to 1")), y = sum_diff,col = factor(probes), shape = factor(probes),group = factor(probes)))+
   xlab("Range of Difficulty Parameters for Screen Items") + ylab("Bias Due to Skip-Logic") +
   stat_summary(geom = "point",fun.y = "mean",size = 8, aes(shape= factor(probes), stroke = 1))+
  stat_summary(geom = "point",fun.y = "mean",size = 8, aes(col= factor(probes)))+
  scale_shape_manual(name = "Number
of Screen
Items",values = c(8,6,1)) +
  scale_colour_manual(name = "Number
of Screen
Items",values = c("blue","green","red"))
  

p + facet_wrap(.~factor(dmin), strip.position = "right")+ theme(plot.title = element_text(size = rel(1))) +
  theme(panel.spacing = unit(1, "lines")) +scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +theme_bw() +
  labs(tag = "Range of Difficulty Parameters for Non-Screen Items") + theme(legend.box.margin=margin(l=20),
                                                                            plot.tag=element_text(angle=-90),
                                                                            plot.tag.position=c(.87, 0.5))+ theme(legend.justification = "top")+
  theme(panel.grid.major.x  =   element_blank(),
        panel.grid.major =   element_line(colour = "black",size=0.12))

  


############## Top orientation for graph 3

p<- ggplot(data = df, aes (x = factor(dminp,c("-3 to -1", "-1 to 1")), y = sum_diff,col = factor(probes), shape = factor(probes),group = factor(probes)))+
  xlab("Range of Difficulty Parameters for Screen Items") + ylab("Bias Due to Skip-Logic") +
  stat_summary(geom = "point",fun.y = "mean",size = 8, aes(shape= factor(probes), stroke = 1))+
  stat_summary(geom = "point",fun.y = "mean",size = 8, aes(col= factor(probes)))+
  ggtitle("Range of Difficulty Parameters for Non-Screen Items")+ 
  scale_shape_manual(name = "Number
of Screen
Items",values = c(8,6,1)) +
  scale_colour_manual(name = "Number
of Screen
Items",values = c("blue","green","red"))

p + facet_wrap(.~factor(dmin), strip.position = "top")+ theme(plot.title = element_text(size = rel(1))) +
  theme(panel.spacing = unit(1, "lines")) +scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +theme_bw() + theme(legend.justification = "top")+
  theme(panel.grid.major.x  =   element_blank(),
        panel.grid.major =   element_line(colour = "black",size=0.12))

################### Graph 4
p <- ggplot(data = df, aes (x = factor(dminp,c("-3 to -1", "-1 to 1")), y = sum_diff,col = factor(probes), shape = factor(probes), group = factor(probes)))+
  xlab("Range of Difficulty Parameters for Screen Items") + ylab("Bias Due to Skip-Logic") +
  stat_summary(geom = "point",fun.y = "mean",size = 8, aes(shape= factor(probes), stroke = 1))+
  ggtitle(
  "Range of Discrimination Parameters")+
  scale_shape_manual(name = "Number
of Screen
Items",values = c(8,6,1)) +
  scale_colour_manual(name = "Number
of Screen
Items",values = c("blue","green","red"))+ theme_bw()
  

p + facet_grid(dmin ~ factor(amin)) + labs(tag = "Range of Difficulty Parameters for Non-Screen Items") + theme(legend.box.margin=margin(l=20),
                                                                   plot.tag=element_text(angle=-90),
                                                                   plot.tag.position=c(.88, 0.5))+ theme(legend.justification = "top")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 14))+
  theme(panel.grid.major.x  =   element_blank(),
        panel.grid.major =   element_line(colour = "black",size=0.12)) 
 
################### Graph 5 
df2 <- read.csv("BlandAltman.csv")
df2$Bias <- factor(df2$Bias, levels = c("Bias = -40","Bias = -20","Bias = -10","Bias = -5","Bias = +5","Bias = +10","Bias = +20","Bias = +40"))

p <- ggplot(df2, aes(x = Speed_Mean, y = Speed_diff, group = factor(Bias), col = factor(Bias))) + stat_smooth(method = "loess", se = FALSE)+
  xlab("Mean of Scores (speed)")+ylab("Difference Between Scores")+labs(col = "Bias")+
  geom_point(size = 1.5, shape = 1) + scale_y_continuous(breaks = round(seq(min(df2$Speed_diff), max(df2$Speed_diff), by = .2),1))+
  theme_bw()

p + facet_wrap(.~factor(Bias_SD))+ theme(legend.justification = "top") +
  scale_colour_manual(name = "Bias",values = c("blue","green","red","purple","yellow","orange","cyan","black"))

################## Graph 6

p <- ggplot(df2, aes(x = Lapses_Mean, y = Lapses_diff, group = factor(Bias), col = factor(Bias))) + stat_smooth(method = "lm", se = FALSE)+
  geom_point(size = 1.5, shape = 1)+xlab("Mean of Scores")+ylab("Difference Between Scores")+labs(col = "Bias")+
  scale_y_continuous(breaks = round(seq(min(df2$Lapses_diff), max(df2$Lapses_diff), by = 1),1)) + theme_bw()

p + facet_wrap(.~factor(Bias_SD))+ theme(legend.justification = "top") +
  scale_colour_manual(name = "Bias",values = c("blue","green","red","purple","yellow","orange","cyan","black"))
########
#ggplot(data = temp, aes(x = test, y = value, shape = factor(CASE.CONTROL))) + 
#  stat_summary(geom = "point", fun.y = "mean", size = 5, aes(col = factor(gender)))+
#  scale_shape_manual(name = "Case/Control", values = c(2,3))+
#  scale_colour_manual(name = "Gender", values = c("red","blue"))+
#  theme(legend.justification = "top")+
##  #scale_x_discrete(labels=c("CPF", "ER40", "SNFB","PMAT24", "SVOLT","VSPLOT24","PCET","CPT"))+ggtitle("Accuracy Z Score by Case Type/Gender")+
#  xlab("Test Name") + ylab("Z score (Accuracy)") +
#  geom_hline(yintercept=0.00, size = .3)


#line plot
#ggplot(data = temp, aes(x = test, y = value, group = interaction(gender,CASE.CONTROL))) + 
#  stat_summary(geom = "line", fun.y = "mean", size = 1, aes(col = factor(gender), linetype = CASE.CONTROL))+
#  scale_x_discrete(labels=c("CPF", "CPT", "ER40","PCET", "PMAT24","SFNB2","SVOLT","VSPLOT"))+ggtitle("Accuracy Z Score by Case Type/Gender")+
#  xlab("Test Name") + ylab("Z score (Accuracy)")+scale_linetype_manual(name = "Case/Control", values = c(8,1))+
#  scale_colour_manual(name = "Gender", values = c("red","blue"))+
#  theme(legend.justification = "top")

#lower is better
#temp$line_group = with(temp,paste(gender,CASE.CONTROL))
#ggplot(data = temp, aes(x = test, y = value, shape =line_group,col=line_group)) + 
#  stat_summary(geom = "point", fun.y = "mean", size = 5)+
#  scale_shape_manual(name = "Gender with Case/Control",values = c(2,1,2,1))+
#  scale_colour_manual(name = "Gender with Case/Control",values = c("red","red","blue","blue")) +
#  stat_summary(geom = "line", fun.y = "mean", size = .25,aes(col=line_group, group = line_group))+ggtitle("Accuracy Z Score by Case Type/Gender")+
#  xlab("Test Name") + ylab("Z score (Accuracy)")+
#  scale_x_discrete(labels=c("CPF", "CPT", "ER40","PCET", "PMAT24","SFNB2","SVOLT","VSPLOT"))+
#  theme(legend.justification = "top")
######

######
#ggplot(data = temp, aes(x = factor(test), y = value, shape = factor(CASE.CONTROL))) + 
#  stat_summary(geom = "point", fun.y = "mean", size = 5, aes(col = factor(gender)))+
#  scale_shape_manual(name = "Case/Control", values = c(2,3))+
#  scale_colour_manual(name = "Gender", values = c("red","blue"))+
#  theme(legend.justification = "top")+
#  scale_x_discrete(labels=c("CPF", "ER40", "SNFB","SCTAP","PMAT24", "SVOLT","VSPLOT24","PCET","CPT"))+ggtitle("Speed Z Score by Case Type/Gender (Lower is Better)")+
#  xlab("Test Name") + ylab("Z score (Speed)") +
#  geom_hline(yintercept=0.00, size = .3)

#line
#ggplot(data = temp, aes(x = test, y = value, group = interaction(gender,CASE.CONTROL))) + 
#  stat_summary(geom = "line", fun.y = "mean", size = 1, aes(col = factor(gender), linetype = CASE.CONTROL))+
#  scale_x_discrete(labels=c("CPF", "ER40", "SNFB","SCTAP","PMAT24", "SVOLT","VSPLOT24","PCET","CPT"))+ggtitle("Accuracy Z Score by Case Type/Gender")+
#  xlab("Test Name") + ylab("Z score (Accuracy)")+scale_linetype_manual(name = "Case/Control", values = c(8,1))+
#  scale_colour_manual(name = "Gender", values = c("red","blue"))+
#  theme(legend.justification = "top")

#temp$line_group = with(temp,paste(gender,CASE.CONTROL))
#ggplot(data = temp, aes(x = test, y = value, shape =line_group,col=line_group)) + 
#  stat_summary(geom = "point", fun.y = "mean", size = 5)+
#  scale_shape_manual(name = "Gender with Case/Control",values = c(2,1,2,1))+
#  scale_colour_manual(name = "Gender with Case/Control",values = c("red","red","blue","blue")) +
##  stat_summary(geom = "line", fun.y = "mean", size = .25,aes(col=line_group, group = line_group))+ggtitle("Speed Z Score by Case Type/Gender")+
#  xlab("Test Name") + ylab("Z score (Accuracy) Lower is better")+
#  scale_x_discrete(labels=c("CPF", "CPT", "ER40","MPRACT","PCET","PMAT24","SCTAP","SFNB2","SVOLT","VSPLOT"))+
#  theme(legend.justification = "top")
