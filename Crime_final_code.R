# =============================================================================
# Title:  Determinants of Crime Rate in US Cities (Updated)
# Author: Tommaso Zipoli, Luca Marchesi, Luca, Orlando
# Date:   14/06/2025
# =============================================================================

# ---- Section 1: Data Loading and Scaling ----
library(cluster)
library(ggcorrplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(FNN)
library(readxl)

# Load and prepare data
data_raw <- read_excel("/Applications/Rproject LMEC/data/dati crime.xlsx")
row.names(data_raw) <- data_raw$cityCode
selected_vars <- c("crmrte", "popden", "lawexpc", "polpc", "offarea", "lpcinc")
data <- data_raw %>% select(all_of(selected_vars)) %>% scale() %>% as.data.frame()

# ---- Section 2: PCA analysis
cov.m <- cov(data)
round(cov.m, digit=2)
cor.m <- cor(data)
round(cor.m, digit=2)
plot_correlation_matrix <- function(data, title) {
  cor_matrix <- cor(data, use = "complete.obs")
  cor_upper <- cor_matrix
  cor_upper[lower.tri(cor_upper)] <- NA
  ggcorrplot(
    cor_upper,
    method = "square",
    type = "full",
    lab = TRUE,
    lab_size = 3,
    tl.cex = 8.3,
    colors = c("#6D9EC1", "white", "#E46726"),
    title = title
  )
}
plot_correlation_matrix(data, "Matrice di Correlazione")
pca <- princomp(data, cor = TRUE)
str(pca)
summary(pca)
eigenval <- (pca$sdev)^2

df <- data.frame(PC = 1:length(eigenval), Eigenvalue = eigenval)

ggplot(df, aes(x = PC, y = Eigenvalue)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(color = "black", size = 2) +
  geom_text(aes(label = round(Eigenvalue, 4)), vjust = -0.8, size = 3.5) +
  labs(title = "Eigenvalues by Principal Component",
       x = "Principal Component", y = "Eigenvalue") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

D <- matrix(0, 6, 6)
diag(D) <- eigenval
C <- pca$loadings
View(C)
loadings <- C %*% sqrt(D)
stand.coeff <- loadings %*% diag(1 / eigenval)
data.s <- scale(data)
score <- data.s %*% stand.coeff

# Compute values
pr.var <- pca$sdev^2
pve <- pr.var / sum(pr.var)
cum_pve <- cumsum(pve)
n <- length(pve)

# Creazione del data frame
df <- data.frame(
  PC = factor(1:n),
  Explained = pve,
  Cumulative = cum_pve
)

# Costruzione del grafico senza istogramma e con curva per la varianza spiegata
ggplot(df, aes(x = PC)) +
  geom_line(aes(y = Explained, group = 1), color = "#50C878", linewidth = 1) +  # Linea rossa per varianza spiegata
  geom_point(aes(y = Explained), color = "#50C878", size = 2) +  # Punti rossi per varianza spiegata
  geom_line(aes(y = Cumulative, group = 1), color = "black", linewidth = 0.8) +  # Linea nera per la cumulativa
  geom_point(aes(y = Cumulative), color = "black", size = 2) +
  scale_y_continuous(
    name = "Proportion of Variance",
    limits = c(0, 1.001),
    breaks = seq(0, 1, 0.1),
    sec.axis = dup_axis(name = NULL)
  ) +
  labs(
    title = "Explained and Cumulative Variance by Principal Component",
    x = "Principal Component"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )

# Extract first two loadings
loadings <- pca$loadings[, 1:2]
load_df <- as.data.frame(loadings)
load_df$Variable <- rownames(load_df)

# Plot Loadings: PC1 vs PC2
plot(load_df$Comp.1, load_df$Comp.2,
     xlab = "Loadings PC1", ylab = "Loadings PC2",
     main = "Variable Loadings: PC1 vs PC2",
     type = "n", xlim = c(-1, 1), ylim = c(-1, 1))
text(load_df$Comp.1, load_df$Comp.2, labels = load_df$Variable, col = "black")

scores_manual <- data.s %*% stand.coeff[, 1:2]

# Step 3: Plot degli score PC1–PC2 per ogni osservazione
score_df <- as.data.frame(scores_manual)
colnames(score_df) <- c("PC1", "PC2")

plot(score_df$PC1, score_df$PC2,
     xlab = "PC1", ylab = "PC2",
     main = "Observations Scores: PC1 vs PC2",
     col = "black", pch = 20)
abline(h = 0, v = 0, col = "gray70")


# ---- Section 3: Factor analysis ----

f1 <- factanal(~ crmrte + popden + lawexpc + polpc + offarea + lpcinc,
               data = data, factors = 1, rotation = "none")
# Fattoriale con 2 fattori
f2 <- factanal(~ ., data = data, factors = 2, rotation = "none")

# Estrazione corretta dei risultati
Chisq2 <- round(f2$STATISTIC, 3)
df2 <- f2$dof
pval2 <- round(f2$PVAL, 4)

# Stampa leggibile
cat("Model with 2 Factors:\n")
cat("Chi-Squared =", Chisq2, "\n")
cat("Degrees of Freedom =", df2, "\n")
cat("p-value =", pval2, "\n")

# Tabella per 1–2 fattori
results <- data.frame(
  Factors = 1:2,
  ChiSq = round(c(f1$STATISTIC, f2$STATISTIC), 3),
  df     = c(f1$dof, f2$dof),
  p_value = round(c(f1$PVAL, f2$PVAL), 4)
)

print(results)
loadings(f2)

# Communalities
comm <- 1 - f2$uniquenesses
percVar <- (sum(comm) / 6) * 100

repcorr <- loadings(f2) %*% t(loadings(f2))
residuals <- round(cor.m - repcorr, 2)

print(comm)
print(percVar)
print(residuals)

# Rotazioni
library(GPArotation)
load_mat <- loadings(f2)

rot_varimax <- Varimax(load_mat)
print(rot_varimax$loadings)

rot_quartimax <- quartimax(load_mat)
print(rot_quartimax$loadings)

rot_oblimin <- oblimin(load_mat)
print(rot_oblimin$loadings)
print(rot_oblimin$Phi)

# Scores
fB <- factanal(~ ., data = data, factors = 2, rotation = "none", scores = "Bartlett")$scores
fT <- factanal(~ ., data = data, factors = 2, rotation = "none", scores = "regression")$scores

head(fB)
head(fT)

# Plot dei factor scores (Bartlett)
plot(fB[, 1], fB[, 2], pch = '.', main = "Factor Scores - Bartlett",
     xlab = "Factor 1", ylab = "Factor 2")
text(fB[, 1], fB[, 2], labels = rownames(data), cex = 0.9, font = 2)

# Tabella riassuntiva dei loadings (2 fattori)
L_raw <- as.matrix(loadings(f2))
L_varimax <- as.matrix(rot_varimax$loadings)
L_quartimax <- as.matrix(rot_quartimax$loadings)
L_oblimin <- as.matrix(rot_oblimin$loadings)

tab_all <- cbind(
  round(L_raw, 3),
  round(L_varimax, 3),
  round(L_quartimax, 3),
  round(L_oblimin, 3)
)

colnames(tab_all) <- c(
  "Raw_F1", "Raw_F2",
  "Varimax_F1", "Varimax_F2",
  "Quartimax_F1", "Quartimax_F2",
  "Oblimin_F1", "Oblimin_F2"
)

print(tab_all)

# ---- Section 4: Hierarchical Cluster Analysis ----
k <- 2
distance <- dist(as.matrix(data))

plot_and_cut <- function(hc, method_name, k) {
  windows()
  plot(hc, main = paste(method_name, "- k =", k))
  rect.hclust(hc, k = k, border = c("gray40", "#50C878"))
  cutree(hc, k = k)
}

hc.c <- hclust(distance)
member.c <- plot_and_cut(hc.c, "Complete Linkage", k)

hc.s <- hclust(distance, method = "single")
member.s <- plot_and_cut(hc.s, "Single Linkage", k)

hc.w <- hclust(distance, method = "ward.D2")
member.w <- plot_and_cut(hc.w, "Ward Method", k)

hc.ce <- hclust(distance, method = "centroid")
member.ce <- plot_and_cut(hc.ce, "Centroid Method", k)

# ---- Section 5: Cross-Tabulations ----
cat("\n--- Cross-tabulations (k =", k, ") ---\n")

cat("\nComplete vs Single:\n")
print(table(Complete = member.c, Single = member.s))

cat("\nComplete vs Ward:\n")
print(table(Complete = member.c, Ward = member.w))

cat("\nComplete vs Centroid:\n")
print(table(Complete = member.c, Centroid = member.ce))

# ---- Section 6: Cluster Centroids (on Scaled Data) ----
cat("\n--- Cluster Centroids (Complete Linkage, k =", k, ") ---\n")
agg_c <- aggregate(data, by = list(cluster = member.c), FUN = mean)
print(agg_c)

# ---- Section 7: K-Means Clustering on Scaled Data ----
set.seed(123)
kmeans_result <- kmeans(data, centers = k, nstart = 25)
data$cluster <- as.factor(kmeans_result$cluster)

# Scatterplots (on scaled data)
vars_to_plot <- c("popden", "lawexpc", "polpc", "offarea", "lpcinc")
plot_list <- lapply(vars_to_plot, function(var) {
  ggplot(data, aes_string(x = var, y = "crmrte", color = "cluster")) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray40") +
    stat_ellipse(aes_string(fill = "cluster"), geom = "polygon", alpha = 0.2, color = NA) +
    scale_color_manual(values = c("1" = "gray40", "2" = "#50C878")) +
    scale_fill_manual(values = c("1" = "gray80", "2" = "#7AC78A")) +
    labs(title = paste("Scaled Crime Rate vs", var),
         x = var,
         y = "Scaled Crime Rate") +
    theme_minimal() +
    theme(legend.position = "bottom")
})

windows()
empty_plot <- ggplot() + theme_void()
grid.arrange(grobs = c(plot_list, list(empty_plot)), nrow = 3, ncol = 2,
             top = "Crime Rate vs Other Variables by K-Means Cluster (Scaled)")

# ---- Section 7: PCA Representation (K-Means Clusters) ----
pca <- prcomp(data[, selected_vars], center = TRUE, scale. = FALSE)
pca_coords <- data.frame(pca$x[, 1:2], cluster = data$cluster)

# Create grid background for PCA cluster zones
x_range <- seq(min(pca_coords$PC1), max(pca_coords$PC1), length.out = 200)
y_range <- seq(min(pca_coords$PC2), max(pca_coords$PC2), length.out = 200)
grid <- expand.grid(PC1 = x_range, PC2 = y_range)

centroids <- aggregate(. ~ cluster, data = pca_coords, FUN = mean)
knn_model <- get.knnx(data = centroids[, c("PC1", "PC2")], query = grid, k = 1)
grid$cluster <- as.factor(centroids$cluster[knn_model$nn.index])

# Plot PCA
windows()
ggplot() +
  geom_tile(data = grid, aes(x = PC1, y = PC2, fill = cluster), alpha = 0.07) +
  geom_point(data = pca_coords, aes(x = PC1, y = PC2, color = cluster), size = 2) +
  stat_ellipse(data = pca_coords, aes(x = PC1, y = PC2, fill = cluster),
               geom = "polygon", alpha = 0.2, color = NA) +
  scale_fill_manual(values = c("1" = "gray80", "2" = "#7AC78A")) +
  scale_color_manual(values = c("1" = "gray40", "2" = "#50C878")) +
  labs(title = "K-Means Clustering in the PCA Space (Scaled Data)",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  theme(legend.position = "none")
