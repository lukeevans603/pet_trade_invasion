##############################
##### EVANS ET AL. 2022. #####
##############################

d_overlap_combo<-read.csv("D_combo_final.csv")


p<-ggplot(data = d_overlap_combo, aes(x = order)) +
  geom_line(aes(y = delta_d), size=1)+#, colour = "Current Climate"), size =1.25) #+
 # geom_line(aes(y = future_d, colour = "Climate 2050"), size = 1.25) +
  #scale_colour_manual("", 
                     # breaks = c("Current Climate", "Climate 2050"),
                      #values = c("Current Climate"="red", "Climate 2050"="blue")) +
  xlab("Species") +
  ylab("Δ Niche overlap")+
  geom_hline(yintercept = 0)+
  #geom_ribbon(aes(ymin= current_d,
               #   ymax= future_d),
             # fill = "grey")+
  theme(aspect.ratio = 1)
p

p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))

fit <- lm(Future_d ~ Current_d - 1, data = d_overlap_combo)
q<-ggplot(d_overlap_combo, aes(x=Current_d, y= Future_d)) + geom_point(stat='identity')+#, colour = "Current Climate"), size =1.25) #+
  # geom_line(aes(y = future_d, colour = "Climate 2050"), size = 1.25) +
  #scale_colour_manual("", 
  # breaks = c("Current Climate", "Climate 2050"),
  #values = c("Current Climate"="red", "Climate 2050"="blue")) +
  xlab("Current Niche Overlap") +
  ylab("Future Niche overlap")+
  geom_smooth(method = lm, se = FALSE, col = 'red', size=2)+
  # annotate(label = sprintf("y = %.3f x\nR² = %.2f", coef(fit), summary(fit)$r.squared),
  #          geom = "text", x = 0.1, y = 0.65, size = 8) +
  geom_abline(col = "grey 70", size =2)+
  #geom_ribbon(aes(ymin= current_d,
  #   ymax= future_d),
  # fill = "grey")+
  theme(aspect.ratio = 1)

q

q<-q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))

??binomial_smooth

c <- ggplot(d_overlap_combo, aes(order)); c2 <- ggplot(d_overlap_combo, aes(delta_d))
q<-ggplot(d_overlap_combo, aes(x = Current_d, y = Future_d))+ geom_bar()
q
p<-ggplot(d_overlap_combo, aes(delta_d)) +
  geom_freqpoly()
p


# Histogram with density plot
ggplot(d_overlap_combo, aes(x=order)) + 
  geom_histogram(aes(y=delta_d), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 


