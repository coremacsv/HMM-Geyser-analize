require(tidyverse)
require(dplyr)
require(magrittr)
require(plotly)
require(mhsmm)
require(ggplot2)
require("viridis")
require(patchwork)
require(ggforce)
require(NbClust)
require(FactoMineR)
require(factoextra)
require(smotefamily)



geyser <- read.table("enter your path", header = FALSE)
colnames(geyser) <- c("WaitingTime", "EruptionDuration", "Code")
#print(geyser)

geyser <- na.omit(geyser)
geyser$Code <- as.factor(geyser$Code)
table(geyser[,3])

print("The codes 1, 2, and 3 (respectively) indicate whether the durations were observed only as short, medium or long, and are based on Table 1 of Azzalini and Bowman (1990). In this analysis we'll assume that code 0 stands for *Not Classified*")



#ERUPTION DURATION AND WAITING TIME DISTRIBUTION WITH HEATMAP

ggplot(geyser, aes(x = EruptionDuration, y = WaitingTime)) + 
  geom_density_2d_filled() +
  geom_density_2d(colour = "black") +
  guides(fill = guide_legend(title = "Level")) +
  geom_point(aes(color = Code), size = 3) +
  scale_color_manual(values = c("#ef4523", "#00aaff","#01ff00", "#c493ff"))


#CORRELATION BETWEEN WAITING TIME AND ERUPTION DURATION

cor(geyser$WaitingTime, geyser$EruptionDuration)


#ERUPTION DURATIONS SORTED

EruptionDurationSorted <- arrange(geyser, desc(EruptionDuration))
EruptionDurationSorted

Top10LongestEruptions <- head(EruptionDurationSorted, 10)
Top10LongestEruptions

max(geyser$WaitingTime)
min(geyser$WaitingTime)


#OBSERVATIONS' CODES DISTRIBUTION WITH RADAR PLOT

ggplot(geyser) + 
  geom_bar(aes(x = Code, fill = Code), stat = "count", width = 0.8) +
  coord_polar() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  labs(title = "Distribuzione dei Codici delle Eruzioni") +
  scale_fill_viridis(discrete = TRUE)


#GENERAL INFO DASHBOARD

WTBoxPlot <- ggplot(geyser) + 
  geom_boxplot(aes(x = WaitingTime, y = Code, fill = Code)) +
  scale_color_manual(values = c("#ef4523", "#00aaff","#01ff00", "#c493ff"))
WTBoxPlot

EDBarPlot <- ggplot(geyser) + 
  geom_bar(aes(x = WaitingTime, fill = Code)) +
  scale_color_manual(values = c("#ef4523", "#00aaff","#01ff00", "#c493ff"))
EDBarPlot

WTEDScatterPlot <- ggplot(geyser) + 
  geom_point(aes(x = WaitingTime, y = EruptionDuration, colour = Code, size = EruptionDuration)) +
  scale_color_manual(values = c("#ef4523", "#00aaff","#01ff00", "#c493ff"))
WTEDScatterPlot

(WTBoxPlot | EDBarPlot) / WTEDScatterPlot



#ERUPTION SCORE BOXPLOT

EruptionScore <- scale(geyser$WaitingTime)/geyser$EruptionDuration
geyser$ES <- EruptionScore


EruptionScorePlot <- plot_ly(geyser, x = ~Code, y = ~ES, text = ~Code, marker = list(color = ~ES, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE), type = "box", color = ~Code, colors = c("#ef4523", "#00aaff","#01ff00", "#c493ff"))


EruptionScorePlot <- EruptionScorePlot %>% layout(scene = list(xaxis = list(title = 'EruptionScore'), yaxis = list(title = 'Code')), annotations = list(
  x = 1.02,
  y = 1.0,
  text = 'Code',
  xref = 'paper',
  yref = 'paper',
  font = list(family="Arial", size=20),
  showarrow = FALSE),
  title = list(
    text = "Eruption Score Boxplot",
    font = list(family = "Arial", size = 20), y=0.98))

EruptionScorePlot

geyser <- geyser[,-4]

#ERUPTION SCALES COMPARISON

EP1 <- ggplot(geyser, aes(x = EruptionDuration, y = WaitingTime)) +
  geom_point() + 
  geom_mark_ellipse(aes(filter = Code == 1,
                        label = 'Bubble Eruptions'))


EP2 <- ggplot(geyser, aes(x = EruptionDuration, y = WaitingTime)) +
  geom_point() + 
  geom_mark_ellipse(aes(filter = Code == 2,
                        label = 'Medium Power Eruptions'))


EP3 <- ggplot(geyser, aes(x = EruptionDuration, y = WaitingTime)) +
  geom_point() + 
  geom_mark_ellipse(aes(filter = Code == 3,
                        label = 'Severe Eruptions'))


(EP1 | EP2 | EP3)



#-------------------------------------------------------------------------------------------

#MARKOV MODELS

geyser <- geyser[,-3]

WaitingTimeMean <- mean(geyser$WaitingTime)
EruptionDurationMean <- mean(geyser$EruptionDuration)


#Two States

K <- 2
EruptionMarkovModelInit <- hmmspec(init = rep(1/K, K),
                                   trans = matrix(1/K, nrow = K, ncol = K),
                                   parms.emis = list(mu = c(50,80), sigma=c(1,1)),
                                   dens.emis = dnorm.hsmm)


EruptionMarkovModel <- hmmfit(geyser[,1], EruptionMarkovModelInit, mstep = mstep.norm)

plot(EruptionMarkovModel$loglik, type = "b", ylab = "Log-likelihood", xlab = "Iteration")

states <- EruptionMarkovModel$yhat

plot(geyser[,1],col=states, main = "WaitingTime")
abline(h=EruptionMarkovModel$model$parms.emission$mu[1])
abline(h=EruptionMarkovModel$model$parms.emission$mu[2],col=1)



#Three States

K <- 3
EruptionMarkovModelInit3S <- hmmspec(init = rep(1/K, K),
                                   trans = matrix(1/K, nrow = K, ncol = K),
                                   parms.emis = list(mu = c(50,80, 110), sigma=c(1,1,1)),
                                   dens.emis = dnorm.hsmm)


EruptionMarkovModel3S <- hmmfit(geyser[,1], EruptionMarkovModelInit3S, mstep = mstep.norm)

plot(EruptionMarkovModel3S$loglik, type = "b", ylab = "Log-likelihood", xlab = "Iteration")

states <- EruptionMarkovModel3S$yhat

plot(geyser[,1],col=states, main = "WaitingTime")
abline(h=EruptionMarkovModel3S$model$parms.emission$mu[1])
abline(h=EruptionMarkovModel3S$model$parms.emission$mu[2],col=2)
abline(h=EruptionMarkovModel3S$model$parms.emission$mu[3],col=3)

#-------------------------------------------------------------------------------------------

remove(geyser)

geyser <- read.table("enter your path ", header = FALSE)
colnames(geyser) <- c("WaitingTime", "EruptionDuration", "Code")

geyser <- na.omit(geyser)
#geyser$Code <- as.numeric(geyser$Code)

#-------------------------------------------------------------------------------------------



# CLUSTERING, PCA AND SMOTE

set.seed(100)

NbClust(data = geyser[, -3], distance = "euclidean", method = "kmeans")
NbClust(data = geyser[, -3], distance = "euclidean", method = "ward.D")
NbClust(data = geyser[, -3], distance = "euclidean", method = "ward.D2")


geyserHCPCClustering <- HCPC(geyser[, -3], graph = FALSE, metric="manhattan", method="ward")

fviz_cluster(geyserHCPCClustering,
             show.clust.cent = TRUE,
             palette = c("darkorange", "darkgreen", "purple"),
             ggtheme = theme_minimal(),
             main = "Geyser Clustering")


geyserOversampled <- SMOTE(geyser, target=as.numeric(geyser$Code==2), K=1, dup_size=0)
geyserOversampled$data[,-4]


afterSMOTEGeyserClustering <- HCPC(geyserOversampled$data[,-c(3,4)], graph = FALSE, metric="manhattan", method="ward")

fviz_cluster(afterSMOTEGeyserClustering,
             show.clust.cent = TRUE,
             palette = c("darkorange", "darkgreen", "purple"),
             ggtheme = theme_minimal(),
             main = "Geyser Clustering")










