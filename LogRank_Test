#input: metrics, EventFreeSurvivaTime and BinaryFirstEvent 
#BinaryFirstEvent <- read.csv("~/TALL/Final/BinaryFirstEvent242.csv")
#metrics <- read.csv("~/TALL/Final/8metrics242.csv")

#metrics 

reference <- BinaryFirstEvent$BinaryFirstEvent3000
surv_object <- Surv(time = BinaryFirstEvent$EventFreeSurvDays, event = BinaryFirstEvent$BinaryFirstEvent3000)

Youden_Function <- function(variable, reference) {
    roc_obj <- roc(reference, variable, quiet = TRUE)
    auc <- roc_obj$auc
    youden_index <- roc_obj$sensitivities + roc_obj$specificities - 1
    optimal_cutoff <- roc_obj$thresholds[which.max(youden_index)]
    return(optimal_cutoff)
}
  

# DataFrame para los clusters binarios por métrica
BinaryClusters_byMetric <- data.frame(Patient = 1:nrow(metrics)) 
CutoffPoints <- data.frame()

# Calcular puntos de corte y asignar grupos binarios
for (x in 1:ncol(metrics)) { 
  variable <- as.numeric(metrics[[x]])
  col_name <- colnames(metrics[x])
  
  optimal_cutoff <- Youden_Function(variable, reference)
  Clusters <- ifelse(variable >= optimal_cutoff, "1", "2")
  BinaryClusters_byMetric[[x]] <- Clusters
  colnames(BinaryClusters_byMetric)[ncol(BinaryClusters_byMetric)] <- paste0(col_name, "_youden")
  
  MetricCutoff <- data.frame(Metric = col_name, optimal_cutoff = optimal_cutoff)
  CutoffPoints <- rbind(CutoffPoints, MetricCutoff)
}


# DataFrame para almacenar resultados
LogRankMetrics3000 <- data.frame(
  Metric = character(),
  N_Group = character(),
  Observed_Group1 = numeric(),
  Observed_Group2 = numeric(),
  Expected_Group1 = numeric(),
  Expected_Group2 = numeric(),
  O_E2_E = numeric(),
  O_E2_V = numeric(),
  Chisq = numeric(),
  p_value = numeric(),
  Range = character(),
  Mean_MAD_Group1 = character(),
  Mean_MAD_Group2 = character(),
  Youden_Cutoff = numeric(),
  stringsAsFactors = FALSE
)

# Evaluación de log-rank test por métrica

for (x in 1:ncol(metrics)) {
  variable_name <- colnames(metrics[x])
  group_col <- paste0(variable_name, "_youden")
  groups <- BinaryClusters_byMetric[[group_col]]
  
  # Realizar log-rank test
  log_rank <- survdiff(surv_object ~ groups, data = BinaryClusters_byMetric)
  
  # Calcular valores adicionales
  variable <- as.numeric(metrics[[variable_name]])
  cutoff <- CutoffPoints[CutoffPoints$Metric == variable_name, "optimal_cutoff"]
  
  range_var <- sprintf("%.4f-%.4f", min(variable, na.rm = TRUE), max(variable, na.rm = TRUE))
  group1 <- variable[groups == "1"]
  group2 <- variable[groups == "2"]
  
  mean_g1 <- mean(group1, na.rm = TRUE)
  mad_g1 <- mad(group1, na.rm = TRUE)
  mean_g2 <- mean(group2, na.rm = TRUE)
  mad_g2 <- mad(group2, na.rm = TRUE)
  
  percent_g1 <- round(sum(groups == "1") / length(groups) * 100, 1)
  percent_g2 <- round(sum(groups == "2") / length(groups) * 100, 1)
  
  # Extraer estadísticas del log-rank test
  N <- log_rank$n
  Observed <- log_rank$obs
  Expected <- log_rank$exp
  Var <- log_rank$var  # Varianza de cada grupo
  Chisq <- log_rank$chisq
  p_value <- pchisq(Chisq, df = 1, lower.tail = FALSE)  # Grado de libertad = 1 (2 grupos)
  
  # Calcular O_E2_E y O_E2_V para ambos grupos
  O_E2_E_Group1 <- (Observed[1] - Expected[1])^2 / Expected[1]
  O_E2_V_Group1 <- (Observed[1] - Expected[1])^2 / Var[1]
  
  O_E2_E_Group2 <- (Observed[2] - Expected[2])^2 / Expected[2]
  O_E2_V_Group2 <- (Observed[2] - Expected[2])^2 / Var[2]
  
  # Calcular los porcentajes
  percent_g1 <- (N[1] / sum(N)) * 100
  percent_g2 <- (N[2] / sum(N)) * 100
  
  # Agregar resultados al dataframe
  LogRankMetrics3000 <- rbind(LogRankMetrics3000, data.frame(
    Metric = variable_name,
    Range = range_var,
    Youden_Cutoff = cutoff,
    Group1 = sprintf("%d (%.1f%%)", N[1], percent_g1),
    Group2 = sprintf("%d (%.1f%%)", N[2], percent_g2),
    Mean_MAD_Group1 = sprintf("%.5f (%.5f)", mean_g1, mad_g1),  # Aquí se especifican 5 decimales
    Mean_MAD_Group2 = sprintf("%.5f (%.5f)", mean_g2, mad_g2),  # Aquí se especifican 5 decimales
    Observed_Group1 = Observed[1],
    Observed_Group2 = Observed[2],
    Expected_Group1 = Expected[1],
    Expected_Group2 = Expected[2],
    O_E2_E_Group1 = O_E2_E_Group1,
    O_E2_V_Group1 = O_E2_V_Group1,
    O_E2_E_Group2 = O_E2_E_Group2,
    O_E2_V_Group2 = O_E2_V_Group2,
    Chisq = Chisq,
    p_value = p_value
  ))
}


write.csv(LogRankMetrics3000, file.path("~/TALL/Final/final/LogRankMetrics3000.csv"), row.names = FALSE)
write.csv(BinaryClusters_byMetric, file.path("~/TALL/Final/final/BinaryClusters3000.csv"), row.names = FALSE)

##################################################################################

reference <- BinaryFirstEvent$BinaryFirstEvent1825
surv_object <- Surv(time = BinaryFirstEvent$EventFreeSurvDays1825, event = BinaryFirstEvent$BinaryFirstEvent1825)





Youden_Function <- function(variable, reference) {
  roc_obj <- roc(reference, variable, quiet = TRUE)
  auc <- roc_obj$auc
  youden_index <- roc_obj$sensitivities + roc_obj$specificities - 1
  optimal_cutoff <- roc_obj$thresholds[which.max(youden_index)]
  return(optimal_cutoff)
}


# DataFrame para los clusters binarios por métrica
BinaryClusters_byMetric <- data.frame(Patient = 1:nrow(metrics)) 
CutoffPoints <- data.frame()

# Calcular puntos de corte y asignar grupos binarios
for (x in 1:ncol(metrics)) { 
  variable <- as.numeric(metrics[[x]])
  col_name <- colnames(metrics[x])
  
  optimal_cutoff <- Youden_Function(variable, reference)
  Clusters <- ifelse(variable >= optimal_cutoff, "1", "2")
  BinaryClusters_byMetric[[x]] <- Clusters
  colnames(BinaryClusters_byMetric)[ncol(BinaryClusters_byMetric)] <- paste0(col_name, "_youden")
  
  MetricCutoff <- data.frame(Metric = col_name, optimal_cutoff = optimal_cutoff)
  CutoffPoints <- rbind(CutoffPoints, MetricCutoff)
}


# DataFrame para almacenar resultados
LogRankMetrics1825 <- data.frame(
  Metric = character(),
  N_Group = character(),
  Observed_Group1 = numeric(),
  Observed_Group2 = numeric(),
  Expected_Group1 = numeric(),
  Expected_Group2 = numeric(),
  O_E2_E = numeric(),
  O_E2_V = numeric(),
  Chisq = numeric(),
  p_value = numeric(),
  Range = character(),
  Mean_MAD_Group1 = character(),
  Mean_MAD_Group2 = character(),
  Youden_Cutoff = numeric(),
  stringsAsFactors = FALSE
)

# Evaluación de log-rank test por métrica

for (x in 1:ncol(metrics)) {
  variable_name <- colnames(metrics[x])
  group_col <- paste0(variable_name, "_youden")
  groups <- BinaryClusters_byMetric[[group_col]]
  
  # Realizar log-rank test
  log_rank <- survdiff(surv_object ~ groups, data = BinaryClusters_byMetric)
  
  # Calcular valores adicionales
  variable <- as.numeric(metrics[[variable_name]])
  cutoff <- CutoffPoints[CutoffPoints$Metric == variable_name, "optimal_cutoff"]
  
  range_var <- sprintf("%.4f-%.4f", min(variable, na.rm = TRUE), max(variable, na.rm = TRUE))
  group1 <- variable[groups == "1"]
  group2 <- variable[groups == "2"]
  
  mean_g1 <- mean(group1, na.rm = TRUE)
  mad_g1 <- mad(group1, na.rm = TRUE)
  mean_g2 <- mean(group2, na.rm = TRUE)
  mad_g2 <- mad(group2, na.rm = TRUE)
  
  percent_g1 <- round(sum(groups == "1") / length(groups) * 100, 1)
  percent_g2 <- round(sum(groups == "2") / length(groups) * 100, 1)
  
  # Extraer estadísticas del log-rank test
  N <- log_rank$n
  Observed <- log_rank$obs
  Expected <- log_rank$exp
  Var <- log_rank$var  # Varianza de cada grupo
  Chisq <- log_rank$chisq
  p_value <- pchisq(Chisq, df = 1, lower.tail = FALSE)  # Grado de libertad = 1 (2 grupos)
  
  # Calcular O_E2_E y O_E2_V para ambos grupos
  O_E2_E_Group1 <- (Observed[1] - Expected[1])^2 / Expected[1]
  O_E2_V_Group1 <- (Observed[1] - Expected[1])^2 / Var[1]
  
  O_E2_E_Group2 <- (Observed[2] - Expected[2])^2 / Expected[2]
  O_E2_V_Group2 <- (Observed[2] - Expected[2])^2 / Var[2]
  
  # Calcular los porcentajes
  percent_g1 <- (N[1] / sum(N)) * 100
  percent_g2 <- (N[2] / sum(N)) * 100
  
  # Agregar resultados al dataframe
  LogRankMetrics1825 <- rbind(LogRankMetrics1825, data.frame(
    Metric = variable_name,
    Range = range_var,
    Youden_Cutoff = cutoff,
    Group1 = sprintf("%d (%.1f%%)", N[1], percent_g1),
    Group2 = sprintf("%d (%.1f%%)", N[2], percent_g2),
    Mean_MAD_Group1 = sprintf("%.5f (%.5f)", mean_g1, mad_g1),  # Aquí se especifican 5 decimales
    Mean_MAD_Group2 = sprintf("%.5f (%.5f)", mean_g2, mad_g2),  # Aquí se especifican 5 decimales
    Observed_Group1 = Observed[1],
    Observed_Group2 = Observed[2],
    Expected_Group1 = Expected[1],
    Expected_Group2 = Expected[2],
    O_E2_E_Group1 = O_E2_E_Group1,
    O_E2_V_Group1 = O_E2_V_Group1,
    O_E2_E_Group2 = O_E2_E_Group2,
    O_E2_V_Group2 = O_E2_V_Group2,
    Chisq = Chisq,
    p_value = p_value
  ))
}


write.csv(LogRankMetrics1825, file.path("~/TALL/Final/final/LogRankMetrics1825.csv"), row.names = FALSE)
write.csv(BinaryClusters_byMetric, file.path("~/TALL/Final/final/BinaryClusters1825.csv"), row.names = FALSE)

########################################################



reference <- BinaryFirstEvent$BinaryFirstEvent730
surv_object <- Surv(time = BinaryFirstEvent$EventFreeSurvDays730, event = BinaryFirstEvent$BinaryFirstEvent730)


Youden_Function <- function(variable, reference) {
  roc_obj <- roc(reference, variable, quiet = TRUE)
  auc <- roc_obj$auc
  youden_index <- roc_obj$sensitivities + roc_obj$specificities - 1
  optimal_cutoff <- roc_obj$thresholds[which.max(youden_index)]
  return(optimal_cutoff)
}


# DataFrame para los clusters binarios por métrica
BinaryClusters_byMetric <- data.frame(Patient = 1:nrow(metrics)) 
CutoffPoints <- data.frame()

# Calcular puntos de corte y asignar grupos binarios
for (x in 1:ncol(metrics)) { 
  variable <- as.numeric(metrics[[x]])
  col_name <- colnames(metrics[x])
  
  optimal_cutoff <- Youden_Function(variable, reference)
  Clusters <- ifelse(variable >= optimal_cutoff, "1", "2")
  BinaryClusters_byMetric[[x]] <- Clusters
  colnames(BinaryClusters_byMetric)[ncol(BinaryClusters_byMetric)] <- paste0(col_name, "_youden")
  
  MetricCutoff <- data.frame(Metric = col_name, optimal_cutoff = optimal_cutoff)
  CutoffPoints <- rbind(CutoffPoints, MetricCutoff)
}


# DataFrame para almacenar resultados
LogRankMetrics730 <- data.frame(
  Metric = character(),
  N_Group = character(),
  Observed_Group1 = numeric(),
  Observed_Group2 = numeric(),
  Expected_Group1 = numeric(),
  Expected_Group2 = numeric(),
  O_E2_E = numeric(),
  O_E2_V = numeric(),
  Chisq = numeric(),
  p_value = numeric(),
  Range = character(),
  Mean_MAD_Group1 = character(),
  Mean_MAD_Group2 = character(),
  Youden_Cutoff = numeric(),
  stringsAsFactors = FALSE
)

# Evaluación de log-rank test por métrica

for (x in 1:ncol(metrics)) {
  variable_name <- colnames(metrics[x])
  group_col <- paste0(variable_name, "_youden")
  groups <- BinaryClusters_byMetric[[group_col]]
  
  # Realizar log-rank test
  log_rank <- survdiff(surv_object ~ groups, data = BinaryClusters_byMetric)
  
  # Calcular valores adicionales
  variable <- as.numeric(metrics[[variable_name]])
  cutoff <- CutoffPoints[CutoffPoints$Metric == variable_name, "optimal_cutoff"]
  
  range_var <- sprintf("%.4f-%.4f", min(variable, na.rm = TRUE), max(variable, na.rm = TRUE))
  group1 <- variable[groups == "1"]
  group2 <- variable[groups == "2"]
  
  mean_g1 <- mean(group1, na.rm = TRUE)
  mad_g1 <- mad(group1, na.rm = TRUE)
  mean_g2 <- mean(group2, na.rm = TRUE)
  mad_g2 <- mad(group2, na.rm = TRUE)
  
  percent_g1 <- round(sum(groups == "1") / length(groups) * 100, 1)
  percent_g2 <- round(sum(groups == "2") / length(groups) * 100, 1)
  
  # Extraer estadísticas del log-rank test
  N <- log_rank$n
  Observed <- log_rank$obs
  Expected <- log_rank$exp
  Var <- log_rank$var  # Varianza de cada grupo
  Chisq <- log_rank$chisq
  p_value <- pchisq(Chisq, df = 1, lower.tail = FALSE)  # Grado de libertad = 1 (2 grupos)
  
  # Calcular O_E2_E y O_E2_V para ambos grupos
  O_E2_E_Group1 <- (Observed[1] - Expected[1])^2 / Expected[1]
  O_E2_V_Group1 <- (Observed[1] - Expected[1])^2 / Var[1]
  
  O_E2_E_Group2 <- (Observed[2] - Expected[2])^2 / Expected[2]
  O_E2_V_Group2 <- (Observed[2] - Expected[2])^2 / Var[2]
  
  # Calcular los porcentajes
  percent_g1 <- (N[1] / sum(N)) * 100
  percent_g2 <- (N[2] / sum(N)) * 100
  
  # Agregar resultados al dataframe
  LogRankMetrics730 <- rbind(LogRankMetrics730, data.frame(
    Metric = variable_name,
    Range = range_var,
    Youden_Cutoff = cutoff,
    Group1 = sprintf("%d (%.1f%%)", N[1], percent_g1),
    Group2 = sprintf("%d (%.1f%%)", N[2], percent_g2),
    Mean_MAD_Group1 = sprintf("%.5f (%.5f)", mean_g1, mad_g1),  # Aquí se especifican 5 decimales
    Mean_MAD_Group2 = sprintf("%.5f (%.5f)", mean_g2, mad_g2),  # Aquí se especifican 5 decimales
    Observed_Group1 = Observed[1],
    Observed_Group2 = Observed[2],
    Expected_Group1 = Expected[1],
    Expected_Group2 = Expected[2],
    O_E2_E_Group1 = O_E2_E_Group1,
    O_E2_V_Group1 = O_E2_V_Group1,
    O_E2_E_Group2 = O_E2_E_Group2,
    O_E2_V_Group2 = O_E2_V_Group2,
    Chisq = Chisq,
    p_value = p_value
  ))
}


write.csv(LogRankMetrics730, file.path("~/TALL/Final/final/LogRankMetrics730.csv"), row.names = FALSE)
write.csv(BinaryClusters_byMetric, file.path("~/TALL/Final/final/BinaryClusters730.csv"), row.names = FALSE)

#################################################################





reference <- BinaryFirstEvent$BinaryFirstEvent365
surv_object <- Surv(time = BinaryFirstEvent$EventFreeSurvDays3655, event = BinaryFirstEvent$BinaryFirstEvent365)


Youden_Function <- function(variable, reference) {
  roc_obj <- roc(reference, variable, quiet = TRUE)
  auc <- roc_obj$auc
  youden_index <- roc_obj$sensitivities + roc_obj$specificities - 1
  optimal_cutoff <- roc_obj$thresholds[which.max(youden_index)]
  return(optimal_cutoff)
}


# DataFrame para los clusters binarios por métrica
BinaryClusters_byMetric <- data.frame(Patient = 1:nrow(metrics)) 
CutoffPoints <- data.frame()

# Calcular puntos de corte y asignar grupos binarios
for (x in 1:ncol(metrics)) { 
  variable <- as.numeric(metrics[[x]])
  col_name <- colnames(metrics[x])
  
  optimal_cutoff <- Youden_Function(variable, reference)
  Clusters <- ifelse(variable >= optimal_cutoff, "1", "2")
  BinaryClusters_byMetric[[x]] <- Clusters
  colnames(BinaryClusters_byMetric)[ncol(BinaryClusters_byMetric)] <- paste0(col_name, "_youden")
  
  MetricCutoff <- data.frame(Metric = col_name, optimal_cutoff = optimal_cutoff)
  CutoffPoints <- rbind(CutoffPoints, MetricCutoff)
}


# DataFrame para almacenar resultados
LogRankMetrics365 <- data.frame(
  Metric = character(),
  N_Group = character(),
  Observed_Group1 = numeric(),
  Observed_Group2 = numeric(),
  Expected_Group1 = numeric(),
  Expected_Group2 = numeric(),
  O_E2_E = numeric(),
  O_E2_V = numeric(),
  Chisq = numeric(),
  p_value = numeric(),
  Range = character(),
  Mean_MAD_Group1 = character(),
  Mean_MAD_Group2 = character(),
  Youden_Cutoff = numeric(),
  stringsAsFactors = FALSE
)

# Evaluación de log-rank test por métrica

for (x in 1:ncol(metrics)) {
  variable_name <- colnames(metrics[x])
  group_col <- paste0(variable_name, "_youden")
  groups <- BinaryClusters_byMetric[[group_col]]
  
  # Realizar log-rank test
  log_rank <- survdiff(surv_object ~ groups, data = BinaryClusters_byMetric)
  
  # Calcular valores adicionales
  variable <- as.numeric(metrics[[variable_name]])
  cutoff <- CutoffPoints[CutoffPoints$Metric == variable_name, "optimal_cutoff"]
  
  range_var <- sprintf("%.4f-%.4f", min(variable, na.rm = TRUE), max(variable, na.rm = TRUE))
  group1 <- variable[groups == "1"]
  group2 <- variable[groups == "2"]
  
  mean_g1 <- mean(group1, na.rm = TRUE)
  mad_g1 <- mad(group1, na.rm = TRUE)
  mean_g2 <- mean(group2, na.rm = TRUE)
  mad_g2 <- mad(group2, na.rm = TRUE)
  
  percent_g1 <- round(sum(groups == "1") / length(groups) * 100, 1)
  percent_g2 <- round(sum(groups == "2") / length(groups) * 100, 1)
  
  # Extraer estadísticas del log-rank test
  N <- log_rank$n
  Observed <- log_rank$obs
  Expected <- log_rank$exp
  Var <- log_rank$var  # Varianza de cada grupo
  Chisq <- log_rank$chisq
  p_value <- pchisq(Chisq, df = 1, lower.tail = FALSE)  # Grado de libertad = 1 (2 grupos)
  
  # Calcular O_E2_E y O_E2_V para ambos grupos
  O_E2_E_Group1 <- (Observed[1] - Expected[1])^2 / Expected[1]
  O_E2_V_Group1 <- (Observed[1] - Expected[1])^2 / Var[1]
  
  O_E2_E_Group2 <- (Observed[2] - Expected[2])^2 / Expected[2]
  O_E2_V_Group2 <- (Observed[2] - Expected[2])^2 / Var[2]
  
  # Calcular los porcentajes
  percent_g1 <- (N[1] / sum(N)) * 100
  percent_g2 <- (N[2] / sum(N)) * 100
  
  # Agregar resultados al dataframe
  LogRankMetrics365 <- rbind(LogRankMetrics365, data.frame(
    Metric = variable_name,
    Range = range_var,
    Youden_Cutoff = cutoff,
    Group1 = sprintf("%d (%.1f%%)", N[1], percent_g1),
    Group2 = sprintf("%d (%.1f%%)", N[2], percent_g2),
    Mean_MAD_Group1 = sprintf("%.5f (%.5f)", mean_g1, mad_g1),  # Aquí se especifican 5 decimales
    Mean_MAD_Group2 = sprintf("%.5f (%.5f)", mean_g2, mad_g2),  # Aquí se especifican 5 decimales
    Observed_Group1 = Observed[1],
    Observed_Group2 = Observed[2],
    Expected_Group1 = Expected[1],
    Expected_Group2 = Expected[2],
    O_E2_E_Group1 = O_E2_E_Group1,
    O_E2_V_Group1 = O_E2_V_Group1,
    O_E2_E_Group2 = O_E2_E_Group2,
    O_E2_V_Group2 = O_E2_V_Group2,
    Chisq = Chisq,
    p_value = p_value
  ))
}



write.csv(LogRankMetrics365, file.path("~/TALL/Final/final/LogRankMetrics365.csv"), row.names = FALSE)
write.csv(BinaryClusters_byMetric, file.path("~/TALL/Final/final/BinaryClusters365.csv"), row.names = FALSE)


View(LogRankMetrics365)
View(LogRankMetrics730)
View(LogRankMetrics1825)
View(LogRankMetrics3000)

pvalues <- data.frame( metrics = LogRankMetrics3000 , 
                       p365 = LogRankMetrics365$p_value, p730= LogRankMetrics730$p_value,
                       p1825 = LogRankMetrics1825$p_value, p3000 = LogRankMetrics3000$p_value )

pvalues <- cbind(pvalues, LogRankMetrics3000[,4:11])

pval <- write.csv(pvalues, stdout(), row.names = FALSE)

kable(LogRankMetrics365, format = "latex", caption = "Mantel-Cox Year 1")
kable(LogRankMetrics730, format = "latex", caption = "Mantel-Cox Year 2")
kable(LogRankMetrics1825, format = "latex", caption = "Mantel-Cox Year 5")
kable(LogRankMetrics3000, format = "latex", caption = "Mantel-Cox Year 9+")


o_E <- cbind(LogRankMetrics365$Observed_Group1 - LogRankMetrics365$Expected_Group1, LogRankMetrics730$Observed_Group1 - LogRankMetrics730$Expected_Group1, 
              LogRankMetrics1825$Observed_Group1 - LogRankMetrics1825$Expected_Group1, LogRankMetrics3000$Observed_Group1 - LogRankMetrics3000$Expected_Group1)

write.csv(o_E,stdout(), row.names = FALSE)
