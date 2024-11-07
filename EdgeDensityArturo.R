# Initial Housekeeping
library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(igraph)

if (!dir.exists("~/TALL/MetricstoVariables")) {
  dir.create("~/TALL/MetricstoVariables")
}
MetricstoVariables <- "~/TALL/MetricstoVariables"



# 1. Data Import

##Single Sample Network Estimates 
SSN_files <- "/datos/ot/anakamura/Nuri/TALL_Ranked_NoNormal"
SSN_list <- list.files(SSN_files)

##Normalized Expression Matrix 
NorMATX <- readRDS("~/TALL/TALL_TARGET_norm.RDS")

##colDATA preguntar a arturo 
DATA  <- read.csv("~/TALL/ClinicalFiles/colData.csv")

##Clinical files 
CI_1 <- read_excel("~/TALL/ClinicalFiles/ALL_CI_1.xlsx")
CI_2 <- read_excel("~/TALL/ClinicalFiles/ALL_CI_2.xlsx")
CI_3 <- read_excel("~/TALL/ClinicalFiles/ALL_CI_3.xlsx")
CI_4 <- read_excel("~/TALL/ClinicalFiles/ALL_CI_4.xlsx")
CI_5 <- read_excel("~/TALL/ClinicalFiles/ALL_CI_5.xlsx")
CI_ALL <- rbind(CI_1,CI_2,CI_3,CI_4,CI_5)

##Sample filtering  (T-cell Acute Lymphoblastic Leukemia).
CI_ALL <- CI_ALL %>% filter(`Cell of Origin` == "T Cell ALL")


# 2. Clinical Files filtering and management 

## Single Sample Network ID
SSN_ORDER <- data.frame(patient = colnames(NorMATX), SSN_ORDER = 1:ncol(NorMATX))

## Patients present in both Clinical Files and Normalized Expression Matrix.
common_patients <- Reduce(intersect, list(colnames(NorMATX), CI_ALL$`TARGET USI`, DATA$patient))
        
## Common patients filtering and Binary Fist Event, no Censored mutate 
CI_ALL <- CI_ALL %>% filter(`TARGET USI` %in% common_patients) %>%
  mutate(`TARGET USI` = factor(`TARGET USI`, levels = common_patients)) %>%
  mutate(BinaryFirstEvent_NoCensored = ifelse(`First Event` == "None" | `First Event` == "Censored", 0, 1)) %>% 
  arrange(`TARGET USI`) 


#########se podría eliminar todo colData####### arturo 
## Filtering clinical files for common patients and Single Sample Network ID assignation   
DATA <- DATA %>% filter( patient %in% common_patients) %>%
  mutate(patient = factor(patient, levels = common_patients)) %>%
  arrange(patient) %>% left_join(SSN_ORDER, by = "patient") #arturo 


## Class change needed 
DATA$patient <- as.character(DATA$patient)
CI_ALL$`TARGET USI` <- as.character(CI_ALL$`TARGET USI`)

## Corroboration before bind 
all.equal(DATA$patient, CI_ALL$`TARGET USI`, common_patients)

## Complete clinical data 
Total_data <- bind_cols(DATA, CI_ALL) %>% mutate(SSN_ORDER = as.numeric(SSN_ORDER)) %>% arrange(SSN_ORDER)
write.csv(Total_data, file.path(MetricstoVariables,"Total_data.csv"))




# 3. Metrics computation 

metrics_summary <- data.frame()


for ( x in SSN_list) {
  Graph <- read.csv(file.path(SSN_files, x)) 
  Graph <- graph_from_data_frame(Graph[, c("Source", "Target")], directed = FALSE)
  
  patient <- sub("_.*", "", x)
  
  diameter <- diameter(Graph)
  edge_density <- edge_density(Graph)
  num_nodes <- vcount(Graph)
  mean_degree <- mean(degree(Graph))
  
  louvain <- cluster_louvain(Graph)
  modularity_score <- modularity(louvain)
  num_communities <- length(unique(membership(louvain)))
  
  metricas <- c(patient,num_communities,diameter,edge_density,num_nodes,modularity_score, mean_degree)
  metrics_summary <- rbind(metricas, metrics_summary)
  
}

cols <- c("patient","num_communities", "diameter" , "edge_density" ,"num_nodes", "modularity_score", "mean_degree")
colnames(metrics_summary) <- cols

metrics_summary <- metrics_summary %>% mutate(patient = as.numeric(patient)) %>% arrange(patient)

write.csv(metrics_summary, file.path(MetricstoVariables,"metrics_summary.csv"))



# 4. Cutoff determination and Clustering by metric 

## Youden Function 

reference <- Total_data$BinaryFirstEvent_NoCensored #arturo 

Youden_Function <- function(cutoff, variable, reference) {
  pred <- ifelse(variable >= cutoff, 1, 0)  
  TP <- sum(pred == 1 & reference == 1)
  FP <- sum(pred == 1 & reference == 0)
  TN <- sum(pred == 0 & reference == 0)
  FN <- sum(pred == 0 & reference == 1)
  
  sensibility <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  

  youden <- sensibility + specificity - 1
  return(youden)
}

## Patients clustering by Youden's cutoff point per metric 

BinaryClusters_byMetric <- data.frame(Patient = metrics_summary$patient) #pregunta a arturo
CutoffPoints <- data.frame()


for (x in 2:ncol(metrics_summary)) { #metrics_summary[,1] not metric; patients' ID
  
  # Arguments for Youden_Function
  variable <- as.numeric(metrics_summary[[x]])
  #reference <- Total_data$BinaryFirstEvent_NoCensored #arturo 
  
  col_name <- colnames((metrics_summary[x]))
  
  # Sequence of cutoff points
  cutoff_points <- seq(min(variable, na.rm = TRUE), max(variable, na.rm = TRUE), length.out = 100)
  
  # Youden indices for each cutoff point
  youden_index <- sapply(cutoff_points, Youden_Function, variable = variable, reference = reference)
  
  # Optimal cutoff point
  optimal_cutoff <- cutoff_points[which.max(youden_index)]
 
  # Cluster assignment 
  Clusters <- ifelse(variable >= optimal_cutoff, "1", "2")
  
  BinaryClusters_byMetric[[x]] <- Clusters
  
  # Column name assignation 
  colnames(BinaryClusters_byMetric)[ncol(BinaryClusters_byMetric)] <- paste0(col_name, "_youden")
  
  MetricCutoff <- data.frame(Metric = col_name, optimal_cutoff= optimal_cutoff )
  CutoffPoints <- rbind(CutoffPoints, MetricCutoff)
  
}

write.csv(BinaryClusters_byMetric, file.path(MetricstoVariables,"BinaryClusters_byMetric.csv"))
write.csv(CutoffPoints, file.path(MetricstoVariables,"CutoffPoints.csv"))



# 5 Evaluation of first event significance 


##Patient order corroboration 
all.equal(BinaryClusters_byMetric$Patient,Total_data$SSN_ORDER)

## BinaryClusters_byMetric and Survival Variables bind, !is.na filter

FirstEvent_Clustred <- BinaryClusters_byMetric %>% 
  mutate(BinaryFirstEvent_NoCensored = Total_data$BinaryFirstEvent_NoCensored) %>%
  mutate(EventFreeSurvDays = Total_data$`Event Free Survival Time in Days`) %>% 
  filter(!is.na(EventFreeSurvDays))              

surv_object <- Surv(time = FirstEvent_Clustred$EventFreeSurvDays, event = FirstEvent_Clustred$BinaryFirstEvent_NoCensored)


## Evaluación de significancia por metrica 

#num_communities_youden just 1 group, survdiff log_rank_communities not made 
# Mantel-Cox test
log_rank_communities <- survdiff(surv_object ~ num_communities_youden, data = FirstEvent_Clustred) # Error 1 cluster 
log_rank_diameter <- survdiff(surv_object ~ diameter_youden, data = FirstEvent_Clustred) #p =.8
log_rank_edge_density <- survdiff(surv_object ~ edge_density_youden, data = FirstEvent_Clustred) #p= .01
log_rank_nodes <- survdiff(surv_object ~ num_nodes_youden, data = FirstEvent_Clustred) #p=.5
log_rank_modularity <- survdiff(surv_object ~ modularity_score_youden, data = FirstEvent_Clustred) #p= .4
log_rank_mean_degree <- survdiff(surv_object ~ mean_degree_youden, data = FirstEvent_Clustred) # p=.01 



LogRankList <- list(log_rank_diameter, log_rank_edge_density, log_rank_nodes, log_rank_modularity, log_rank_mean_degree)

LogRankMetrics<- data.frame(
  Group = character(),
  N = integer(),
  Observed = numeric(),
  Expected = numeric(),
  O_E2_E = numeric(),
  O_E2_V = numeric(),
  Chisq = numeric(),
  DF = integer(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop para iterar sobre cada objeto de survdiff en LogRankList
for (i in 1:length(LogRankList)) {
  surv_object <- LogRankList[[i]]
  
  # Extraer los valores por grupo
  group_names <- names(surv_object$n)  # Nombre de cada grupo
  N <- surv_object$n                   # Tamaño de cada grupo
  Observed <- surv_object$obs          # Eventos observados
  Expected <- surv_object$exp          # Eventos esperados
  O_E2_E <- (Observed - Expected)^2 / Expected
  O_E2_V <- (Observed - Expected)^2 / surv_object$var
  
  # Extraer estadísticas generales
  Chisq <- surv_object$chisq           # Chi-cuadrado
  DF <- length(group_names) - 1        # Grados de libertad
  p_value <- pchisq(Chisq, DF, lower.tail = FALSE)  # Valor p
  
  # Almacenar en el data frame LogRankMetrics
  for (j in seq_along(group_names)) {
    LogRankMetrics <- rbind(LogRankMetrics, data.frame(
      Group = group_names[j],
      N = N[j],
      Observed = Observed[j],
      Expected = Expected[j],
      O_E2_E = O_E2_E[j],
      O_E2_V = O_E2_V[j],
      Chisq = ifelse(j == 1, Chisq, NA),  # Chi-cuadrado en la primera fila del grupo
      DF = ifelse(j == 1, DF, NA),        # DF en la primera fila del grupo
      p_value = ifelse(j == 1, p_value, NA)  # p-valor solo en la primera fila
    ))
  }
}

print(subset(LogRankMetrics, p_value < 0.05))
write.csv(LogRankMetrics, file.path(MetricstoVariables,"LogRankMetrics.csv"), row.names = FALSE)




############### WSS para kmeans ###############

#Metrics scaling
#Patient 137 has no Event Free Survival
metrics_summary <- metrics_summary[-137,]
NumericMtx <- apply(metrics_summary, 2, as.numeric)
ScaledMetrics <- scale(NumericMtx)
ScaledMetrics <- ScaledMetrics[,3:ncol(ScaledMetrics)]
dim(ScaledMetrics)

############################################
set.seed(123) 
kmeans_result <- kmeans(NormMetrics, centers = 2)
FirstEvent_Clustred$Kmeans <- kmeans_result$cluster


surv_object <- Surv(time = FirstEvent_Clustred$EventFreeSurvDays730, event = FirstEvent_Clustred$BinaryFirstEvent_NoCensored730)

#youden igual mean degree= edge density 


all.equal(FirstEvent_Clustred$mean_degree_youden, FirstEvent_Clustred$edge_density_youden)

EdgeD_MeanD730 <- survfit(surv_object ~ edge_density_youden, data = FirstEvent_Clustred)

# Crear el modelo Kaplan-Meier
KaplanEdge <- ggsurvplot(EdgeD_MeanD730, data = FirstEvent_Clustred, 
           pval = TRUE,             
           risk.table = TRUE,      
           xlab = "Event Free Survival in Days", 
           ylab = "First event",
           legend.title = "Edge Density / Mean Degree 730",
           censor = FALSE)

write.csv(KaplanEdge, file.path(MetricstoVariables,"KaplanEdge.csv"))





FirstEvent_Clustred <- FirstEvent_Clustred %>% mutate(EventFreeSurvDays730 = EventFreeSurvDays) %>% mutate(BinaryFirstEvent_NoCensored730 = BinaryFirstEvent_NoCensored)



FirstEvent_Clustred$BinaryFirstEvent_NoCensored730[FirstEvent_Clustred$EventFreeSurvDays730 > 730] <- 0
FirstEvent_Clustred$EventFreeSurvDays730[FirstEvent_Clustred$EventFreeSurvDays730 > 730] <- 730







Cutoffread.csv("~/TALL/MetricstoVariables//CutoffPoints.csv")











library(pROC)

# Supongamos que tienes el data.frame llamado 'data'
# Aquí un ejemplo de creación de un data.frame similar al tuyo
data <- data.frame(
  Avrge_degree = c("7.24375226367258", "9.15750915750916", "7.02000702000702", "7.56143667296786", "8.12345678901234"),
  BinaryNoCensored = c(0, 0, 0, 0, 1)
)

# Convertir Avrge_degree a tipo numérico
kmDF$Kmeans <- as.numeric(kmDF$Kmeans)

FirstEvent_Clustred$BinaryFirstEvent_NoCensored730 <- as.numeric(FirstEvent_Clustred$BinaryFirstEvent_NoCensored730)
metrics_summary$edge_density<- as.numeric(metrics_summary$edge_density)[-137]

dim(kmDF)
# Calcular la curva ROC
roc_result <- roc(FirstEvent_Clustred$BinaryFirstEvent_NoCensored730, metrics_summary$edge_density[-137])


dim(FirstEvent_Clustred)
# Calcular el AUC
auc_value <- auc(roc_result)
print(paste("AUC:", auc_value))

# Graficar la curva ROC
plot(roc_result, main = "Curva ROC", col = "blue")
abline(a = 0, b = 1, lty = 2, col = "red")

# Graficar la curva ROC
plot(roc_result, main = "AUC edgeden730", col = "purple")
text(x = 0.8, y = 0.8, labels = paste("AUC =", round(auc_value, 3)), col = "black", cex = 1)











