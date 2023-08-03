
require(tidyverse)
require(ggfortify)
require(vegan)

#color palettes - colorblind friendly
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#E69F00", "#D55E00","56B4E9" , "#CC79A7", "#009E73", "#0072B2","#F0E442" , "#000000")

#Import data and filter dataset
geno_name <- read.delim("~/GitHub/HostMutants2021_C86meliloti/data/genotype_name.txt", sep="\t", header=TRUE)
freq_all<-read.table("~/GitHub/HostMutants2021_C86meliloti/data/freq.tsv", header=TRUE)
#C:\Users\Lab Manager\Documents\GitHub\HostMutants2021_C86meliloti\data
#initial = M_C86_NLA-D 
#culture gow-ups = M_C86_TYA-D
#exudate experiment = X_
#alfalfa S&R = M_888-890
#Wild types = M_860, M_869

#Remove initial inoculum community, TY grow-ups, exudate experiment, and wonky alfalfa sample
freq_filtered <- freq_all %>% filter(!str_detect(sample_id, "M_C86_")) %>% filter(!str_detect(sample_id, "X_")) %>% filter(sample_id != "M_888F_S168")




#Calculate Fitness
#Subset initial inoculum community and calculate mean
freq_initial <- freq_all %>% filter(str_detect(sample_id, "M_C86_NL"))
freq_initial_mean <- as_tibble(apply(data.frame(freq_initial[,c(-1:-1)]), MARGIN = 2, mean))

#Transform raw frequency data to fitness data. log2(treatment_freq/initial_freq)
fit<-as_tibble(cbind(freq_filtered[,c(1:1)],log2(mapply("/", freq_filtered[,c(-1:-1)], freq_initial_mean))))

#After running the above line it converts values from num to chr... I think this is because the cbind is joing a chr col and a num col ? Whatever the reason the below 3 lines converts matrix values to num
fit <- as.data.frame(apply(fit, 2, as.numeric))
sample_id <- freq_filtered[,c(1:1)]
fit<- cbind(sample_id, fit[,c(2:87)])

fit <- fit %>% mutate_if(is.numeric, ~ifelse(abs(.) == Inf,-8,.)) 

#Separate sample_id 
fit <- fit %>% mutate(genotype_id = substr(sample_id, 3,5), rep = substr(sample_id, 6,6))
fit$genotype_id <- as.factor(fit$genotype_id)

#Calculate median fitness and plot PCA
#Take Median of fitness
fit.med <- fit %>% group_by(genotype_id) %>% summarise_if(is.numeric,median)

#PCA to localize host genotypes based on strain community
pfit_h<-prcomp(fit.med[,-1],center=TRUE, scale. = FALSE)
summary(pfit_h)
hosts<-fit.med[,1]
pfit_h_results<-data.frame(cbind(hosts,pfit_h$x))
pfit_h_results

#Plot PCA results, uses ggfortify
pfit_h_results_geno <- inner_join(geno_name, pfit_h_results, by = "genotype_id")

p <- autoplot(pfit_h, data=pfit_h_results_geno, colour='background', label=TRUE, label.label="name", frame=TRUE, frame.type='norm') +
  theme_classic()
p

#A17
#Calculate median fitness and plot PCA
#Take Median of fitness

fit.A17 <- inner_join(fit, geno_name, by="genotype_id")
  
fit.med.A17 <- filter(fit.A17, background=="A17") %>% group_by(genotype_id) %>% summarise_if(is.numeric,median)
#PCA to localize host genotypes based on strain community
pfit_h<-prcomp(fit.med.A17[,c(-1,-88)],center=TRUE, scale. = FALSE)
summary(pfit_h)
hosts<-fit.med.A17[,1]
pfit_h_results<-data.frame(cbind(hosts,pfit_h$x))
pfit_h_results

#Plot PCA results, uses ggfortify
pfit_h_results_geno <- inner_join(geno_name, pfit_h_results, by = "genotype_id")

p <- autoplot(pfit_h, data=pfit_h_results_geno, colour='type', label=TRUE, label.label="name") +
  theme_classic()
p

#R108
#Calculate median fitness and plot PCA
#Take Median of fitness

fit.R108 <- inner_join(fit, geno_name, by="genotype_id")

fit.med.R108 <- filter(fit.R108, background=="R108") %>% group_by(genotype_id) %>% summarise_if(is.numeric,median)
#PCA to localize host genotypes based on strain community
pfit_h<-prcomp(fit.med.R108[,c(-1,-88)],center=TRUE, scale. = FALSE)
summary(pfit_h)
hosts<-fit.med.R108[,1]
pfit_h_results<-data.frame(cbind(hosts,pfit_h$x))
pfit_h_results

#Plot PCA results, uses ggfortify
pfit_h_results_geno <- inner_join(geno_name, pfit_h_results, by = "genotype_id")

p <- autoplot(pfit_h, data=pfit_h_results_geno, colour='type', label=TRUE, label.label="name") +
  theme_classic()
p

#MS
#Calculate median fitness and plot PCA
#Take Median of fitness

fit.MS <- inner_join(fit, geno_name, by="genotype_id")

fit.med.MS <- filter(fit.MS, background=="MS") %>% group_by(genotype_id) %>% summarise_if(is.numeric,median)
#PCA to localize host genotypes based on strain community
pfit_h<-prcomp(fit.med.MS[,c(-1,-88)],center=TRUE, scale. = FALSE)
summary(pfit_h)
hosts<-fit.med.MS[,1]
pfit_h_results<-data.frame(cbind(hosts,pfit_h$x))
pfit_h_results

#Plot PCA results, uses ggfortify
pfit_h_results_geno <- inner_join(geno_name, pfit_h_results, by = "genotype_id")

p <- autoplot(pfit_h, data=pfit_h_results_geno, colour='type', label=TRUE, label.label="name") +
  theme_classic()
p

#geom_boxplot(outlier.shape=NA) +
#  geom_jitter(aes(x=name, y=cfu_plant, color=type)) +
#  scale_y_continuous(trans='log10') +
#  labs(x="Plant Genotype",y="log(cfus per plant)")+
#  coord_flip() +
#  scale_color_manual(values = cbPalette) +
#  theme_classic() + theme(legend.position="top", text = element_text(size=18))


#freq table with sample info
#Separate sample_id 
freq <- freq_filtered %>% mutate(genotype_id = substr(sample_id, 3,5), rep = substr(sample_id, 6,6))
freq <- inner_join(freq, geno_name, by="genotype_id")

freq_div <-data.frame(shannon=diversity(as.matrix(freq[,c(-1:-1,-87:-93)]), index = "shannon"),Trt=factor(freq$background),Type=factor(freq$type), Genotype=factor(freq$genotype_id))

invsimpson = diversity(as.matrix(freq[,c(-1:-1,-87:-93)]), index = "invsimpson")
freq_div <- freq_div %>% add_column(invsimp = invsimpson)

anova(lm(shannon~Trt*Type,data=freq_div)) # First pass anova shows 
summary(lm(shannon~Trt*Type,data=freq_div)) # Examination of coeffecients suggests T

anova(lm(invsimpson~Trt*Type,data=freq_div)) # First pass anova shows 
summary(lm(invsimpson~Trt*Type,data=freq_div)) # Examination of coeffecients suggests T

p<-ggplot(data=freq_div, aes(x=Trt, y=shannon, fill=Type)) +
  geom_boxplot() +
  geom_point() + 
  theme_classic()
p

p<-ggplot(data=freq_div, aes(x=Trt, y=invsimpson, fill=Type)) +
  geom_boxplot() +
  geom_point() + 
  theme_classic()
p



#freq_zero <-df_long %>% mutate(value = ifelse(value < -1, 0, value))

#HEAT MAP
#df <- melt(df)
df_long <- pivot_longer(fit,"MAG5":"MAG761A")
df_long <- inner_join(df_long, geno_name, by="genotype_id")
#colnames(df) <- c("x", "y", "value")

p<-ggplot(df_long, aes(x=value)) + 
  geom_histogram(color="black", fill="white")
p

a17 <- df_long %>% filter(background=="A17")
p<- ggplot(a17, aes(x = genotype_id, y = name.x, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "blue",
                      high = "red",
                      guide = "colorbar") 

#ggsave(plot = p , width = 12, height = 12, dpi = 600, filename = "wild.png")
  
#scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) 

#coord_fixed()



###RDA - borrowing code from mesocosom_C68meliloti
# RDA analysis for fitness and plot RDA1 vs RDA2 and RDA3 vs RDA4
fit_details <- inner_join(fit, geno_name, by="genotype_id")
# Run a few different RDA Models
rda1<-rda(fit_details[,c(2:87)]~Trt*Time,meso_fit, scale=TRUE) #Full Model w/ interactions
rda2<-rda(fit_details[,c(2:87)]~Trt+Time,meso_fit, scale=TRUE) # Additive Model
rda3<-rda(fit_details[,c(2:87)]~Trt,meso_fit, scale=TRUE) # Treatment only
rda4<-rda(fit_details[,c(2:87)]~Time,meso_fit, scale=TRUE) # Time only

# Rsquared gives you an overall idea of well the overall model adjusted for the number of predictors you have
mymodels<-data.frame(data=c("all"),model=c("Trt*Time","Trt+Time","Trt","Time"),
                     Radj=c(RsquareAdj(rda1)$adj.r.squared,RsquareAdj(rda2)$adj.r.squared,RsquareAdj(rda3)$adj.r.squared,RsquareAdj(rda4)$adj.r.squared),
                     Prop.Exp=c(summary(rda1)$constr.chi/summary(rda1)$tot.chi,summary(rda2)$constr.chi/summary(rda2)$tot.chi,summary(rda3)$constr.chi/summary(rda3)$tot.chi,summary(rda4)$constr.chi/summary(rda4)$tot.chi))
mymodels
# The adjusted r-squared for the model allowing for an interaction between Time and Treatment has the highest rsquared even when penalized for the additional predictors
# Run a permutational significance test for the explanatory variables in the full model
# In this case the inference from the anova aligns with the inference from the adjusted R squared 

rda.model.all<-anova(rda1, step=1000, perm.max=1000, by= "terms") # The interaction is significant, but temperature has the largest effect
write.table(rda.model.all,file = "../tables/RDA_All_TimexTrt.tsv",sep="\t",col.names = NA,row.names=TRUE)

#Plot the FULL model with the interaction #

#Make the window view wider
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))

# Initiate the plot (Save format:FigX_Analysis_Data_Model)
pdf(file = "../figures/Fig2_RDA_All_TimexTrt.pdf",useDingbats = FALSE, height=4, width=4)
rdaplot <- ordiplot(rda1, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-2,2),ylim=c(-3,6),scaling=1, xlab=paste("RDA 1 (",round(summary(rda1)$cont$importance[2,1],2)*100,"% var.)",sep=""), 
                    ylab=paste("RDA 2 (",round(summary(rda1)$cont$importance[2,2],2)*100,"% var.)",sep=""))
points(rda1,"wa", cex=0.8,col=paste(mycols[meso_fit$Trt,]$cols),pch=myshapes[meso_fit$Time,]$shapes)
legend("topright", legend = levels(meso_fit$Trt), bty = "n",
       col = paste(mycols$cols),pt.bg = paste(mycols$cols),pch=15, title= "Trt")
legend("topleft", legend = levels(meso_fit$Time), bty = "n",
       pch = myshapes$shapes, title="Time")

rdaplot <- ordiplot(rda1,display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",
                    xlim=c(-2,2),ylim=c(-3,6),scaling=1, xlab=paste("RDA 3 (",round(summary(rda1)$cont$importance[2,3],2)*100,"% var.)",sep=""), 
                    ylab=paste("RDA 4 (",round(summary(rda1)$cont$importance[2,4],2)*100,"% var.)",sep=""))
points(rda1,"wa",choices = c(3,4), cex=0.8,col=paste(mycols[meso_fit$Trt,]$cols),pch=myshapes[meso_fit$Time,]$shapes)

dev.off()




