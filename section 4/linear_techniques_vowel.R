rm(list = ls())
library(readr) ; library(tidyverse) ; library(corrplot) ; library(MASS)

train <- read.csv("Linear Methods for Classification/vowel_train.csv", sep = ",")
train$x.10. <- gsub("\\\\", "", train$x.10.)
train$x.10. <- as.numeric(gsub("}", "", train$x.10.))
train$x.10 <- train$x.10. ; train$x.10. <- NULL

test <- read.csv("Linear Methods for Classification/vowel_test.csv", sep = ",")
test$x.10. <- gsub("\\\\", "", test$x.10.)
test$x.10. <- as.numeric(gsub("}", "", test$x.10.))
test$x.10 <- test$x.10. ; test$x.10. <- NULL
head(train) ; summary(train)

corr_matrix <- cor(train[,-1]) # correlation matrix
corrplot(corr_matrix, method = "square", addCoef.col = 'black')  # associated plot

train[paste0("y_", 1:11)] <- lapply(1:11, function(k) { as.integer(train$y == k)})
test[paste0("y_", 1:11)] <- lapply(1:11, function(k) { as.integer(test$y == k)})

x_vars <- paste0("x.", 1:10) # string vector of the name of the features

# Regression Indicator Matrix ------

models <- vector("list", 11) # create an empty list to store the 11 regressions

for (k in 1:11) { # loop through the 11 groups
  y_var <- paste0("y_", k) # string vector of the name of group k 
  
  formula_k <- as.formula(
    paste(y_var, "~", paste(x_vars, collapse = " + ")) # formula regress y_k on the features
  )
  
  models[[k]] <- lm(formula_k, data = train) # fit the reg
}

fitted_train <- as.data.frame(
  lapply(models, fitted) # retrieve the fitted values for each model
) # we store them in a data.frame

names(fitted_train) <- paste0("y_", 1:11) # rename the columns for clarity

# the predicted outcome is the one with the highest fitted values 
prediction_train <- max.col(fitted_train, ties.method = "first")

# Training error rate is ~0.48
train_error_rate <- mean(train$y != prediction_train) ; train_error_rate

# we obtain the fitted values of the model on the test data
fitted_test <- as.data.frame(
  lapply(models, function(model) predict(model, newdata = test))
)

names(fitted_test) <- paste0("y_", 1:11)

prediction_test <- max.col(fitted_test, ties.method = "first")

# the test error rate is ~0.67
test_error_rate <- mean(test$y != prediction_test) ; test_error_rate

# LDA ------ 

lda_model <- lda(y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10,
                 data = train) ; lda_model

lda_prediction_train <- predict(lda_model, newdata = train)
lda_prediction_test <- predict(lda_model, newdata = test)

# training error rate is ~ 0.32
lda_train_error_rate <- mean(train$y != lda_prediction_train$class) ; lda_train_error_rate
# test error rate is ~ 0.56
lda_test_error_rate <- mean(test$y != lda_prediction_test$class) ; lda_test_error_rate

# QDA ---------

qda_model <- qda(y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10,
                 data = train) ; qda_model

qda_prediction_train <- predict(qda_model, newdata = train)
qda_prediction_test <- predict(qda_model, newdata = test)

# training error rate is ~ 0.01
qda_train_error_rate <- mean(train$y != qda_prediction_train$class) ; qda_train_error_rate
# test error rate is ~ 0.53
qda_test_error_rate <- mean(test$y != qda_prediction_test$class) ; qda_test_error_rate

