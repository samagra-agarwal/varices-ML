library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(caret)
library(reshape2)
library(xgboost)
library(SHAPforxgboost)

d <- read_xlsx(path = "varices.xlsx")

t.train <- d %>% filter(group == 1)
t.test <- d %>% filter(group == 2)
validation <- d %>% filter(group==3)

baveno <- d %>% filter (group>1)  %>% filter(PLT<150 | LSMi>20)

t.train <-  t.train %>% mutate(ETIOLOGY = factor(ETIOLOGY), SEX = factor(SEX), HIGH.RISK = factor(HIGH.RISK), 
                                           BLEED = factor(BLEED))
t.test<-  t.test %>% mutate(ETIOLOGY = factor(ETIOLOGY), SEX = factor(SEX), HIGH.RISK = factor(HIGH.RISK), 
                                        BLEED = factor(BLEED))
validation<-  validation %>% mutate(ETIOLOGY = factor(ETIOLOGY), SEX = factor(SEX), HIGH.RISK = factor(HIGH.RISK), 
                                    BLEED = factor(BLEED))

baveno <- baveno %>% mutate(ETIOLOGY = factor(ETIOLOGY), SEX = factor(SEX), HIGH.RISK = factor(HIGH.RISK), 
                            BLEED = factor(BLEED))


library(ipred)
library(doSNOW)
library(e1071)

# 10-fold crossvalidartion to be perforemd 3 times and perform grid search
train.control = trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "grid")

#grid search for xgboost hyperparameters
tune.grid = expand.grid(eta = c(0.05, 0.075, 0.1), nrounds = c(50, 75, 100), max_depth = 6:8, 
                        min_child_weight = c(2.0, 2.25, 2.5), colsample_bytree = c(0.3, 0.4, 0.5),
                        gamma = 0, subsample = 1)

#multithreading
cl = makeCluster(spec = 4, type = "SOCK")

#register multithread cluster
registerDoSNOW(cl)

#actual model construction
t.train2 <- t.train[,-12]
caret.cv = train(BLEED ~ ., data = t.train2, method = "xgbTree", 
                 tuneGrid = tune.grid, trControl = train.control)

stopCluster(cl)

caret.cv

pred2 = predict(caret.cv, t.train)
pred3 = predict(caret.cv, t.test)
pred4 = predict(caret.cv, validation)
pred5 = predict(caret.cv, baveno)

confusionMatrix(pred2, t.train$BLEED)
confusionMatrix(pred3, t.test$BLEED)
confusionMatrix(pred4, validation$BLEED)
confusionMatrix(pred5, baveno$BLEED)

d <- d %>% mutate(hr2 = 1 - as.numeric(HIGH.RISK))
d <- d %>% mutate(hrfactor  = factor(hr2), blfactor = factor(BLEED))
confusionMatrix(d$hrfactor, d$blfactor)

t.test$pred = as.numeric(pred3)
t.test <- t.test %>%mutate(stratifier = ifelse(pred==1 & HIGH.RISK==1, 4, 
                                                           ifelse(HIGH.RISK==1, 3, ifelse(pred==1, 2, 1))))


validation$pred = as.numeric(pred4)
validation <- validation %>%mutate(stratifier = ifelse(pred==1 & HIGH.RISK==1, 4, 
                                                       ifelse(HIGH.RISK==1, 3, ifelse(pred==1, 2, 1))))


t.train$pred = as.numeric(pred2)
t.train <- t.train %>%mutate(stratifier = ifelse(pred==1 & HIGH.RISK==1, 4, 
                                                             ifelse(HIGH.RISK==1, 3, ifelse(pred==1, 2, 1))))

baveno$pred = as.numeric(pred5)
baveno <- baveno %>%mutate(stratifier = ifelse(pred==1 & HIGH.RISK==1, 4, 
                                               ifelse(HIGH.RISK==1, 3, ifelse(pred==1, 2, 1))))

s.train <- Surv(event = t.train$BLEED==0, time = t.train$TIME.TO.BLEED)
s.test <- Surv(event = t.test$BLEED==0, time = t.test$TIME.TO.BLEED)
s.valid <- Surv(event = validation$BLEED==0, time = validation$TIME.TO.BLEED)
s.all <- Surv(event = d$BLEED==0, time = d$TIME.TO.BLEED)
s.baveno <- Surv(event = baveno$BLEED ==0, time = baveno$TIME.TO.BLEED)

summary(sf.train, times = c(360, 1080))
summary(sf.test, times = c(360, 1080))
summary(sf.valid, times = c(360, 1080))
summary(survfit(s.all ~ d$HIGH.RISK), times = c(360, 1080))
summary(survfit(s.baveno ~ baveno$stratifier), times = c(360, 1080))

sf.train = survfit(s.train ~ t.train$stratifier, data = t.train)
sf.test = survfit(s.test ~ t.test$stratifier, data = t.test)
sf.valid = survfit(s.valid ~ validation$stratifier, data = validation)
sf.baveno <- survfit(s.baveno ~ baveno$stratifier, data = baveno)


p1 <- ggsurvplot(sf.test, break.time.by = 360, xlim= c(0, 1830), risk.table = TRUE, 
                 surv.median.line = "hv", pval = TRUE, risk.table.height = 0.35,ylab = "Bleed-free survival",
                 legend.labs = c("True low-risk", "EGD low-risk/ ML high-risk", 
                                 "EGD high-risk/ ML low-risk", "True high-risk"), 
                 title = "Internal Validation cohort")
p2 <- ggsurvplot(sf.train, break.time.by = 360, xlim= c(0, 1830), risk.table = TRUE, 
                 surv.median.line = "hv", pval = TRUE, risk.table.height = 0.35,ylab = "Bleed-free survival",
                 legend.labs = c("True low-risk", "EGD low-risk/ ML high-risk", 
                                 "EGD high-risk/ ML low-risk", "True high-risk"), title = "Derivation cohort")
p3 <- ggsurvplot(sf.valid, break.time.by = 360, xlim= c(0, 1830), risk.table = TRUE, 
                 surv.median.line = "hv", pval = TRUE, risk.table.height = 0.35, ylab = "Bleed-free survival",
                 legend.labs = c("True low-risk", "EGD low-risk/ ML high-risk", 
                                 "EGD high-risk/ ML low-risk", "True high-risk"), 
                 title = "External validation cohort")

p4 <- ggsurvplot(sf.baveno, break.time.by = 360, xlim= c(0, 1830), risk.table = TRUE, 
                 surv.median.line = "hv", pval = TRUE, risk.table.height = 0.35, ylab = "Bleed-free survival",
                 legend.labs = c("True low-risk", "EGD low-risk/ ML high-risk", 
                                 "EGD high-risk/ ML low-risk", "True high-risk"), 
                 title = "Baveno high-risk validation cohort")

p2
p1
p3



ggsave("validation.png", print(p1), width = 9.5, height = 4.5, units = "in", dpi = 600)
ggsave("derivation.png", print(p2), width = 9.5, height = 4.5, units = "in", dpi = 600)
ggsave("external.png", print(p3), width = 9.5, height = 4.5, units = "in", dpi = 600)
ggsave("baveno.png", print(p4), width = 9.5, height = 4.5, units = "in", dpi = 600)


train$cohort = x




y_var <-  "BLEED"
dataX <- t.train2 %>% select(-y_var, -"group")
# hyperparameter tuning results
param_dart <- list(objective = "binary:logistic",  # For regression
                   # booster = "dart",
                   nrounds = 100,
                   eta = 0.075,
                   max_depth = 7,
                   gamma = 0,
                   subsample = 1,
                   colsample_bytree = 0.4
)

mod <- xgboost::xgboost(data = data.matrix(dataX), 
                        label = as.matrix(t.train2[[y_var]]), 
                        xgb_param = param_dart, nrounds = param_dart$nrounds,
                        verbose = FALSE, nthread = parallel::detectCores() - 2,
                        early_stopping_rounds = 8)

# To return the SHAP values and ranked features by mean|SHAP|
shap_values <- shap.values(xgb_model = mod, X_train = data.matrix(dataX))
# The ranked features by mean |SHAP|
shap_long <- shap.prep(xgb_model = mod, X_train = data.matrix(dataX))
new_labels <- list(AGE = "Age in years", SEX = "Gender", ETIOLOGY = "Etiology", LSMi = "Baseline LSM",
                   HIGH.RISK = "Endoscopic classification", HB = "Hemoglobin", BIL = "Bilirubin",
                   PLT = "Platelet count", ALB = "Albumin", INR = "INR", Creatinine = "Creatinine")
fplot <- shap.plot.summary(shap_long)

ggsave("shap.png", print(fplot), width = 8.5, height = 5.6, units = "in", dpi = 600)

lapply(names(dataX), shap.plot.dependence, data_long = shap_long)

library(DiagrammeR)
gr <- xgb.plot.multi.trees(mod, render = TRUE)
ggsave("tree.png", print(gr), width = 8.5, height = 5.6, units = "in", dpi = 600)
