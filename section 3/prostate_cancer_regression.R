rm(list = ls())
library(stargazer) ; library(stats) ; library(tidyverse) ; library(readr) 
library(corrplot) ; library(psych) ; library(glmnet) 
set.seed(123) # set seed for reproducibility
# Prostate cancer data
# LEAST SQUARES -----------
# Data come from Stamey et al. (1989). 
# Look a correlation between the level of prostate-specific antigen and
# a number of clinical measures in men who were about to receive a radical prostatectomy. 
# The variables are :
# log cancer volume (lcavol) 
# log prostate weight (lweight)
# age
# log of the amount of benign prostatic hyperplasia (lbph)
# seminal vesicle invasion (svi)
# log of capsular penetration (lcp),
# Gleason score (gleason)
# percent of Gleason scores 4 or 5 (pgg45).
# We first load the data and examine it. 
data <- read_delim("Linear Methods for Regression/prostate.data.txt", delim = "\t") # import the data
names(data)[1] <- "ID" # rename the ID column
head(data) ; summary(data) # visualize
vars <- c("lcavol", "lweight", "age", "lbph", "svi", "lcp", "gleason", "pgg45") # features
data <- data %>% mutate_at(append(vars, "lpsa"), as.numeric) # convert the the variables to numeric.
data <- data %>% mutate_at(vars, ~ scale(.) %>% as.vector) # standardize the features
corr_matrix <- cor(data[,vars]) # correlation matrix
corrplot(corr_matrix, method = "square", addCoef.col = 'black')  # associated plot
# scatter plot matrices (SPLOM), 
# with bivariate scatter plots below the diagonal, 
# histograms on the diagonal, 
#and the Pearson correlation above the diagonal. 
pairs.panels(data[, vars], hist.col = "lightblue", density = F, smooth = F, 
             ellipses = F) 
# We fit a linear model to the log of prostate-specific antigen, lpsa, after
# Randomly split the data into training set of size 67 and test set of size 30. 
train_data <- subset(data, train == "TRUE") ; nrow(train_data)
#train_data <- train_data %>% mutate_at(vars, ~ scale(.) %>% as.vector) # standardize the features
formula <- reformulate(vars, response = "lpsa") ; formula
model <- lm(formula, data = train_data)
summary(model)
# we compute the associated MSPE on the test data
test_data <- subset(data, train != "TRUE") ; nrow(test_data)
#test_data <- test_data %>% mutate_at(vars, ~ scale(.) %>% as.vector) # standardize the features
predictions <- predict(model, newdata = test_data)
MSPE_1 <- sum((predictions - test_data$lpsa)^2)/nrow(test_data) ; MSPE_1 # approx 0.52
std_1 <- sqrt(var((predictions - test_data$lpsa)^2)/nrow(test_data)) ; std_1 # std deviation
# In contrast, if we use the mean lpsa in the training set, 
# we would get a MSPE of
MSPE_2 <- sum((mean(train_data$lpsa) - test_data$lpsa)^2) / nrow(test_data) ; MSPE_2 # 1.06
MSPE_1/MSPE_2 - 1 # linear model reduces "base error rate" by about 50\%

# RIDGE ----------
X <- model.matrix(formula, data = train_data)[,-1] # design matrix from train data
y <- train_data$lpsa # outcome vector from train data
ridge <- cv.glmnet(X, y, nfolds = 10, trace.it = 1, alpha = 0) # Ridge with CV10
coef(ridge, s = "lambda.min") # look at the model associated to lambda minimizing CV10
X_test <- model.matrix(formula, data = test_data)[,-1] # design matrix from test data
y_hat_ridge <- predict(ridge, newx = X_test, s = "lambda.min") # fitted values on test data
y_test <- test_data$lpsa # outcome vector from test data
MSPE_ridge <- sum((y_hat_ridge - y_test)^2) / nrow(test_data) ; MSPE_ridge # MSPE
std_ridge  <- as.numeric(sqrt(var((y_hat_ridge - y_test)^2) / nrow(test_data))) ; std_ridge # sd
MSPE_ridge / MSPE_1 - 1
# ridge reduces test error of LS estimates by a small amount. (around 5%)

df_lambda <- function(lambda, X){ # function to retrieve the effective df associated to each lambda
  hat_lambda <- X %*% solve(t(X) %*% X + lambda * diag(ncol(X))) %*% t(X)
  return(tr(hat_lambda))
}

lambdas <- ridge$glmnet.fit$lambda # the sequence of lambdas
coefs <- as.matrix(t(coef(ridge$glmnet.fit))) # coefs associated with each lambda value
coefs <- coefs[, -1]  # we don't care abobut the intercept
dfs <- sapply(lambdas, function(l) df_lambda(l, X)) # effective df at each lambda value
df_min <- df_lambda(ridge$lambda.min, X) # effective df at lambda min
df_1se <- df_lambda(ridge$lambda.1se, X) # effective df at lambda.1se
plot_data <- as.data.frame(coefs) # Reshape to long format for ggplot
plot_data$df <- dfs
plot_data_long <- plot_data %>%
  pivot_longer(cols = -df, names_to = "variable", values_to = "coefficient")
ggplot(plot_data_long, aes(x = df, y = coefficient, color = variable)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = df_min, linetype = "dashed", color = "black") +
  geom_vline(xintercept = df_1se, linetype = "dashed", color = "black") + 
  annotate("text", x = df_min, y = max(plot_data_long$coefficient), 
           label = paste0("df = ", round(df_min, 2)), hjust = -0.1, size = 3.5) +
  labs(x = "df(lambda)", y = "Coefficients",
       title = "Ridge Coefficient Paths vs. Effective Degrees of Freedom",
       color = "Variable") +
  theme_minimal()
# the first dashed line is the effective df for lambda 1.se
# the second dashed line is the effective df for lambda.min. 


# CV error data, aligned to the same lambda sequence as before
cv_data <- data.frame(
  lambda = ridge$lambda,
  df     = dfs,          # reuse the df vector you already computed
  cvm    = ridge$cvm,    # mean CV error at each lambda
  cvlo   = ridge$cvlo,   # lower bound of the CV error (mean - 1 SE)
  cvup   = ridge$cvup    # upper bound of the CV error (mean + 1 SE)
)

ggplot(cv_data, aes(x = df, y = cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), color = "grey70", width = 0.1) +
  geom_line(color = "steelblue", linewidth = 0.6) +
  geom_point(color = "steelblue", size = 1.5) +
  geom_vline(xintercept = df_min, linetype = "dashed", color = "black") +
  geom_vline(xintercept = df_1se, linetype = "dashed", color = "black") +
  annotate("text", x = df_min, y = max(cv_data$cvup), 
           label = paste0("df.min = ", round(df_min, 2)), hjust = -0.1, size = 3.5) +
  annotate("text", x = df_1se, y = max(cv_data$cvup) * 0.95, 
           label = paste0("df.1se = ", round(df_1se, 2)), hjust = -0.1, size = 3.5) +
  labs(x = "df(lambda)", y = "Cross-Validation Error (MSE)",
       title = "Ridge CV Error vs. Effective Degrees of Freedom") +
  theme_minimal()