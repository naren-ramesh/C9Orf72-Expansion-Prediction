# ============================================================
# Least Absolute Shrinkage and Selection Operator Logistic 
# Regression to predict pathogenic C9orf72 repeat expansion status.
# Files needed: C9orf72 CpG Annotations.csv, C9orf72 CpGs.csv
# These files were provided on GitHub and can be used to replicate these results.
# ============================================================

# ============================================================
# Libraries: Load all packages needed for analysis and plotting
# ============================================================
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)   # EPICv2 array annotation
library(IlluminaHumanMethylationEPICv2manifest)        # EPICv2 manifest file
library(RColorBrewer)                                  # Color palettes
library(missMethyl)                                    # Methylation-specific utilities
library(minfiData)                                     # Example methylation data
library(ggplot2)                                       # Visualization
library(DMRcate)                                       # Differential methylation regions
library(Gviz)                                          # Genomic data visualization
library(rpart)                                         # Decision tree modeling
library(ggplot2)                                       # Repeated load (safe but redundant)
library(randomForest)                                  # Random forest modeling
library(rpart.plot)                                    # Plotting decision trees
library(rpart)                                         # Repeated load
library(caret)                                         # Machine learning utilities
library(tidyr)                                         # Data wrangling
library(glmnet)                                        # Regularized regression (Lasso/Elastic Net)
library("viridis")                                     # Color scales
library(ggsignif)                                      # Significance bars in ggplot

# ============================================================
# Data input: C9orf72 CpG annotation and methylation matrix
# ============================================================
c9orfanno <- read.csv("C9orf72 CpG Annotations.csv", row.names = 1)
dat       <- read.csv("C9orf72 CpGs.csv")

# ============================================================
# Initialize empty data frames to store accuracy and error rates
# for each array platform: EPICv2, EPICv1, 450K, and 27K
# ============================================================
all.acc.dat   <- as.data.frame(matrix(nrow = 100, ncol = 4))
all.typei.dat <- as.data.frame(matrix(nrow = 100, ncol = 4))
all.typeii.dat<- as.data.frame(matrix(nrow = 100, ncol = 4))
colnames(all.acc.dat) <- colnames(all.typei.dat) <- colnames(all.typeii.dat) <-
  c("EPICv2", "EPICv1", "Methyl450K", "Methyl27K")

# ============================================================
# EPICv2 model: First pass training/testing with 70/30 split
# ============================================================
set.seed(123)  # For reproducibility
train_index <- sample(seq_len(nrow(dat)), size = 0.7 * nrow(dat))
train_data  <- dat[train_index, ]
test_data   <- dat[-train_index, ]

# Prepare predictors (X) and response (y)
x_train <- model.matrix(exp ~ ., data = train_data)[, -1]
y_train <- train_data$exp
x_test  <- model.matrix(exp ~ ., data = test_data)[, -1]
y_test  <- test_data$exp

# Fit L1-regularized logistic regression (Lasso) with cross-validation
set.seed(123)
cv_model <- cv.glmnet(x_train, y_train, family = "binomial",
                      alpha = 1, nfolds = 5)

# Inspect model performance and coefficients
cat("Best lambda:", cv_model$lambda.min, "\n")
print(cv_model)
best_lambda  <- cv_model$lambda.min
coefficients <- coef(cv_model, s = best_lambda)
print(coefficients)

# Identify features with non-zero coefficients (important predictors)
important_vars <- rownames(coefficients)[which(coefficients != 0)]
cat("Important variables:", important_vars, "\n")

# Predictions and accuracy assessment
predictions.prob <- predict(cv_model, s = best_lambda, newx = x_test, type = "response")
predictions      <- predict(cv_model, s = best_lambda, newx = x_test, type = "class")
accuracy         <- mean(predictions == y_test)
cat("Test Accuracy:", accuracy, "\n")

# Confusion matrix
conf_matrix <- table(Actual = y_test, Predicted = predictions)
print(conf_matrix)

# ============================================================
# EPICv2 model: 100 resamples for accuracy, Type I, and Type II error rates
# ============================================================
all.acc   <- c()
all.typei <- c()
all.typeii<- c()
for (i in 1:100) {
  set.seed(i)
  train_index <- sample(seq_len(nrow(dat)), size = 0.7 * nrow(dat))
  train_data  <- dat[train_index, ]
  test_data   <- dat[-train_index, ]
  
  x_train <- model.matrix(exp ~ ., data = train_data)[, -1]
  y_train <- train_data$exp
  x_test  <- model.matrix(exp ~ ., data = test_data)[, -1]
  y_test  <- test_data$exp
  
  cv_model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
  best_lambda <- cv_model$lambda.min
  predictions <- predict(cv_model, s = best_lambda, newx = x_test, type = "class")
  
  all.acc <- c(all.acc, mean(predictions == y_test))
  conf_matrix <- table(Actual = y_test, Predicted = predictions)
  
  FP <- conf_matrix[2,1]  # False Positives
  TN <- conf_matrix[2,2]  # True Negatives
  FN <- conf_matrix[1,2]  # False Negatives
  TP <- conf_matrix[1,1]  # True Positives
  
  all.typei  <- c(all.typei,  FP / (FP + TN))   # Type I error rate
  all.typeii <- c(all.typeii, FN / (TP + FN))   # Type II error rate
}
summary(all.acc)
summary(all.typei)
summary(all.typeii)

# Store results in summary data frames
all.acc.dat[, 1]   <- all.acc
all.typei.dat[, 1] <- all.typei
all.typeii.dat[, 1]<- all.typeii

# Plot EPICv2 distributions for accuracy and error rates
plot(density(all.acc),   main = "EPICv2 Accuracy")
plot(density(all.typei), main = "EPICv2 Type I Error Rate")
plot(density(all.typeii),main = "EPICv2 Type II Error Rate")

# Identify important coefficients linked to array annotation
important.coefs <- cbind(c9orfanno[, c("Methyl450_Loci", "Methyl27_Loci", "EPICv1_Loci")],
                         model.coef = matrix(coefficients)[-1, ])
important.coefs[important.coefs$model.coef != 0, ]

# ============================================================
# EPICv1 model: First pass training/testing with 70/30 split
# ============================================================
set.seed(123)  # For reproducibility
epicv1.dat   <- dat[, c(1, which(c9orfanno$EPICv1_Loci != "") + 1)]
train_index  <- sample(seq_len(nrow(epicv1.dat)), size = 0.7 * nrow(epicv1.dat))
train_data   <- epicv1.dat[train_index, ]
test_data    <- epicv1.dat[-train_index, ]

# Prepare predictors and response
x_train <- model.matrix(exp ~ ., data = train_data)[, -1]
y_train <- train_data$exp
x_test  <- model.matrix(exp ~ ., data = test_data)[, -1]
y_test  <- test_data$exp

# Fit Lasso logistic regression
set.seed(123)
cv_model    <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
best_lambda <- cv_model$lambda.min
coefficients <- coef(cv_model, s = best_lambda)

# Inspect model and features
cat("Best lambda:", best_lambda, "\n")
print(cv_model)
print(coefficients)
important_vars <- rownames(coefficients)[which(coefficients != 0)]
cat("Important variables:", important_vars, "\n")

# Predictions and accuracy
predictions.prob <- predict(cv_model, s = best_lambda, newx = x_test, type = "response")
predictions      <- predict(cv_model, s = best_lambda, newx = x_test, type = "class")
accuracy         <- mean(predictions == y_test)
cat("Test Accuracy:", accuracy, "\n")
conf_matrix <- table(Actual = y_test, Predicted = predictions)
print(conf_matrix)

# ============================================================
# EPICv1: 100 resamples
# ============================================================
all.acc   <- c()
all.typei <- c()
all.typeii<- c()
for (i in 1:100) {
  set.seed(i)
  epicv1.dat <- dat[, c(1, which(c9orfanno$EPICv1_Loci != "") + 1)]
  train_index <- sample(seq_len(nrow(epicv1.dat)), size = 0.7 * nrow(epicv1.dat))
  train_data  <- epicv1.dat[train_index, ]
  test_data   <- epicv1.dat[-train_index, ]
  
  x_train <- model.matrix(exp ~ ., data = train_data)[, -1]
  y_train <- train_data$exp
  x_test  <- model.matrix(exp ~ ., data = test_data)[, -1]
  y_test  <- test_data$exp
  
  cv_model    <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
  best_lambda <- cv_model$lambda.min
  predictions <- predict(cv_model, s = best_lambda, newx = x_test, type = "class")
  
  all.acc <- c(all.acc, mean(predictions == y_test))
  conf_matrix <- table(Actual = y_test, Predicted = predictions)
  
  FP <- conf_matrix[2,1]
  TN <- conf_matrix[2,2]
  FN <- conf_matrix[1,2]
  TP <- conf_matrix[1,1]
  
  all.typei  <- c(all.typei,  FP / (FP + TN))
  all.typeii <- c(all.typeii, FN / (TP + FN))
}
summary(all.acc)
summary(all.typei)
summary(all.typeii)

# Store results and plot
all.acc.dat[, 2]   <- all.acc
all.typei.dat[, 2] <- all.typei
all.typeii.dat[, 2]<- all.typeii
plot(density(all.acc),   main = "EPICv1 Accuracy")
plot(density(all.typei), main = "EPICv1 Type I Error Rate")
plot(density(all.typeii),main = "EPICv1 Type II Error Rate")


# ============================================================
# 450K model: First pass training/testing
# ============================================================
set.seed(123)
methyl450k.dat <- dat[, c(1, which(c9orfanno$Methyl450_Loci != "") + 1)]
train_index    <- sample(seq_len(nrow(methyl450k.dat)), size = 0.7 * nrow(methyl450k.dat))
train_data     <- methyl450k.dat[train_index, ]
test_data      <- methyl450k.dat[-train_index, ]

x_train <- model.matrix(exp ~ ., data = train_data)[, -1]
y_train <- train_data$exp
x_test  <- model.matrix(exp ~ ., data = test_data)[, -1]
y_test  <- test_data$exp

set.seed(123)
cv_model    <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
best_lambda <- cv_model$lambda.min
coefficients <- coef(cv_model, s = best_lambda)
cat("Best lambda:", best_lambda, "\n")
print(cv_model)
print(coefficients)
important_vars <- rownames(coefficients)[which(coefficients != 0)]
cat("Important variables:", important_vars, "\n")

predictions.prob <- predict(cv_model, s = best_lambda, newx = x_test, type = "response")
predictions      <- predict(cv_model, s = best_lambda, newx = x_test, type = "class")
accuracy         <- mean(predictions == y_test)
cat("Test Accuracy:", accuracy, "\n")
conf_matrix <- table(Actual = y_test, Predicted = predictions)
print(conf_matrix)

# ============================================================
# 450K: 100 resamples
# ============================================================
all.acc   <- c()
all.typei <- c()
all.typeii<- c()
for (i in 101:200) {
  set.seed(i)
  methyl450k.dat <- dat[, c(1, which(c9orfanno$Methyl450_Loci != "") + 1)]
  train_index <- sample(seq_len(nrow(methyl450k.dat)), size = 0.7 * nrow(methyl450k.dat))
  train_data  <- methyl450k.dat[train_index, ]
  test_data   <- methyl450k.dat[-train_index, ]
  
  x_train <- model.matrix(exp ~ ., data = train_data)[, -1]
  y_train <- train_data$exp
  x_test  <- model.matrix(exp ~ ., data = test_data)[, -1]
  y_test  <- test_data$exp
  
  cv_model    <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
  best_lambda <- cv_model$lambda.min
  predictions <- predict(cv_model, s = best_lambda, newx = x_test, type = "class")
  
  all.acc <- c(all.acc, mean(predictions == y_test))
  conf_matrix <- table(Actual = y_test, Predicted = predictions)
  
  FP <- conf_matrix[2,1]
  TN <- conf_matrix[2,2]
  FN <- conf_matrix[1,2]
  TP <- conf_matrix[1,1]
  
  all.typei  <- c(all.typei,  FP / (FP + TN))
  all.typeii <- c(all.typeii, FN / (TP + FN))
}
summary(all.acc)
summary(all.typei)
summary(all.typeii)

# Store results and plot
all.acc.dat[, 3]   <- all.acc
all.typei.dat[, 3] <- all.typei
all.typeii.dat[, 3]<- all.typeii
plot(density(all.acc),   main = "Methyl450K Accuracy")
plot(density(all.typei), main = "Methyl450K Type I Error Rate")
plot(density(all.typeii),main = "Methyl450K Type II Error Rate")


# ============================================================
# 27K model: First pass training/testing
# ============================================================
set.seed(123)
methyl27k.dat <- dat[, c(1, which(c9orfanno$Methyl27_Loci != "") + 1)]
train_index   <- sample(seq_len(nrow(methyl27k.dat)), size = 0.7 * nrow(methyl27k.dat))
train_data    <- methyl27k.dat[train_index, ]
test_data     <- methyl27k.dat[-train_index, ]

x_train <- model.matrix(exp ~ ., data = train_data)[, -1]
y_train <- train_data$exp
x_test  <- model.matrix(exp ~ ., data = test_data)[, -1]
y_test  <- test_data$exp

set.seed(123)
cv_model    <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
best_lambda <- cv_model$lambda.min
coefficients <- coef(cv_model, s = best_lambda)
cat("Best lambda:", best_lambda, "\n")
print(cv_model)
print(coefficients)
important_vars <- rownames(coefficients)[which(coefficients != 0)]
cat("Important variables:", important_vars, "\n")

predictions.prob <- predict(cv_model, s = best_lambda, newx = x_test, type = "response")
predictions      <- predict(cv_model, s = best_lambda, newx = x_test, type = "class")
accuracy         <- mean(predictions == y_test)
cat("Test Accuracy:", accuracy, "\n")
conf_matrix <- table(Actual = y_test, Predicted = predictions)
print(conf_matrix)

# ============================================================
# 27K: 100 resamples
# ============================================================
all.acc   <- c()
all.typei <- c()
all.typeii<- c()
for (i in 1:100) {
  set.seed(i)
  methyl27k.dat <- dat[, c(1, which(c9orfanno$Methyl27_Loci != "") + 1)]
  train_index <- sample(seq_len(nrow(methyl27k.dat)), size = 0.7 * nrow(methyl27k.dat))
  train_data  <- methyl27k.dat[train_index, ]
  test_data   <- methyl27k.dat[-train_index, ]
  
  x_train <- model.matrix(exp ~ ., data = train_data)[, -1]
  y_train <- train_data$exp
  x_test  <- model.matrix(exp ~ ., data = test_data)[, -1]
  y_test  <- test_data$exp
  
  cv_model    <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
  best_lambda <- cv_model$lambda.min
  predictions <- predict(cv_model, s = best_lambda, newx = x_test, type = "class")
  
  all.acc <- c(all.acc, mean(predictions == y_test))
  conf_matrix <- table(Actual = y_test, Predicted = predictions)
  
  # Special handling if a class is absent in predictions
  if (sum(predictions == "expanded") == 0) {
    conf_matrix <- cbind(expanded = c(0, 0), conf_matrix)
  }
  
  FP <- conf_matrix[2,1]
  TN <- conf_matrix[2,2]
  FN <- conf_matrix[1,2]
  TP <- conf_matrix[1,1]
  
  all.typei  <- c(all.typei,  FP / (FP + TN))
  all.typeii <- c(all.typeii, FN / (TP + FN))
}
summary(all.acc)
summary(all.typei)
summary(all.typeii)

# Store results and plot
all.acc.dat[, 4]   <- all.acc
all.typei.dat[, 4] <- all.typei
all.typeii.dat[, 4]<- all.typeii

# ============================================================
# Combine accuracy results for all platforms for violin plots
# ============================================================
all.acc.dat.long <- data.frame(
  value = c(all.acc.dat$EPICv2,
            all.acc.dat$EPICv1,
            all.acc.dat$Methyl450K,
            all.acc.dat$Methyl27K),
  group = factor(rep(c("EPICv2", "EPICv1", "Methyl450K", "Methyl27K"), each = 100),
                 levels = c("EPICv2", "EPICv1", "Methyl450K", "Methyl27K"))
)

ggplot(all.acc.dat.long, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, color = c("white", "gray", "gray", "black"), outlier.colour = NA) +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  theme_bw() +
  labs(title = "Accuracy of Predictor Models Based on Array Type",
       x = "Array Type", y = "Accuracy") +
  geom_signif(comparisons = list(c("EPICv2", "EPICv1"),
                                 c("EPICv2", "Methyl450K"),
                                 c("EPICv1", "Methyl450K"),
                                 c("Methyl450K", "Methyl27K"),
                                 c("EPICv1", "Methyl27K"),
                                 c("EPICv2", "Methyl27K")),
              y_position = c(1.015, 1.027, 1.041, 1.051, 1.063, 1.075),
              tip_length = 0.01)

# ============================================================
# Combine Type I results for all platforms for violin plots
# ============================================================
all.typei.dat.long <- data.frame(
  value = c(all.typei.dat$EPICv2,
            all.typei.dat$EPICv1,
            all.typei.dat$Methyl450K,
            all.typei.dat$Methyl27K),
  group = factor(rep(c("EPICv2", "EPICv1", "Methyl450K", "Methyl27K"), each = 100)))
all.typei.dat.long$group <- factor(all.typei.dat.long$group, levels = c("EPICv2", "EPICv1", "Methyl450K", "Methyl27K"))
head(all.typei.dat.long)


ggplot(all.typei.dat.long, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, color = c("white", "gray", "gray", "black"), outlier.colour = NA) +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  theme_bw() +
  labs(title = "Type I Error Rate of Predictor Models Based on Array Type",
       x = "Array Type",
       y = "Type I Error Rate") +
  geom_signif(comparisons = list(c("EPICv2", "EPICv1"),
                                 c("EPICv2", "Methyl450K"),
                                 c("EPICv1", "Methyl450K"),
                                 c("Methyl450K", "Methyl27K"),
                                 c("EPICv1", "Methyl27K"),
                                 c("EPICv2", "Methyl27K")),
              y_position = c(0.018, 0.022, 0.026, 0.051, 0.054, 0.057), 
              tip_length = 0.01) 

# ============================================================
# Combine Type II results for all platforms for violin plots
# ============================================================
all.typeii.dat.long <- data.frame(
  value = c(all.typeii.dat$EPICv2,
            all.typeii.dat$EPICv1,
            all.typeii.dat$Methyl450K,
            all.typeii.dat$Methyl27K),
  group = factor(rep(c("EPICv2", "EPICv1", "Methyl450K", "Methyl27K"), each = 100)))
all.typeii.dat.long$group <- factor(all.typeii.dat.long$group, levels = c("EPICv2", "EPICv1", "Methyl450K", "Methyl27K"))

ggplot(all.typeii.dat.long, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, color = c("white", "gray", "gray", "black"), outlier.colour = NA) +
  scale_fill_viridis(option = "mako", discrete = TRUE) +
  theme_bw() +
  labs(title = "Type II Error Rate of Predictor Models Based on Array Type",
       x = "Array Type",
       y = "Type II Error Rate") +
  geom_signif(comparisons = list(c("EPICv2", "EPICv1"),
                                 c("EPICv2", "Methyl450K"),
                                 c("EPICv1", "Methyl450K"),
                                 c("Methyl450K", "Methyl27K"),
                                 c("EPICv1", "Methyl27K"),
                                 c("EPICv2", "Methyl27K")),
              y_position = c(0.45, 1.15, 1.23, 1.29, 1.36, 1.43),
              tip_length = 0.01)