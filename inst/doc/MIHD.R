### R code from vignette source 'MIHD.Rnw'

###################################################
### code chunk number 1: MIHD.Rnw:33-34
###################################################
library('MIHD')


###################################################
### code chunk number 2: MIHD.Rnw:37-40
###################################################
data(datagaussian)
data(databinary)
data(datapoisson)


###################################################
### code chunk number 3: MIHD.Rnw:43-45
###################################################
sum(is.na(datagaussian))
sum(is.na(datagaussian))/dim(datagaussian)[1]


###################################################
### code chunk number 4: MIHD.Rnw:50-51
###################################################
imp=MIdurr(data=datagaussian,method="lasso",family="gaussian",m=5)


###################################################
### code chunk number 5: MIHD.Rnw:54-57
###################################################
imp$n.missing
imp$col.missing
head(imp$impute)


###################################################
### code chunk number 6: MIHD.Rnw:62-67 (eval = FALSE)
###################################################
## MIdurr(data=datagaussian,method="adaptive lasso",family="gaussian")
## MIdurr(data=datagaussian,method="blasso",family="gaussian")
## MIdurr(data=databinary,method="lasso",family="binary")
## MIdurr(data=databinary,method="adaptive lasso",family="binary")
## MIdurr(data=databinary,method="blasso",family="binary")


###################################################
### code chunk number 7: MIHD.Rnw:70-72 (eval = FALSE)
###################################################
## MIiurr(data=datagaussian,method="lasso",family="gaussian")
## MIiurr(data=databinary,method="adaptive lasso",family="binary")


###################################################
### code chunk number 8: MIHD.Rnw:79-80
###################################################
data(GMdata)


###################################################
### code chunk number 9: MIHD.Rnw:84-85
###################################################
GM.imp=GMdurr(data=GMdata)


###################################################
### code chunk number 10: MIHD.Rnw:90-92
###################################################
imp=as.mice(data=GMdata,imp=GM.imp)
head(imp$nmis)


###################################################
### code chunk number 11: MIHD.Rnw:95-96
###################################################
imp$imp$z1


###################################################
### code chunk number 12: MIHD.Rnw:99-100 (eval = FALSE)
###################################################
## complete(imp)


###################################################
### code chunk number 13: MIHD.Rnw:105-106
###################################################
fit=with(imp,lm(z1~z2+z4))


###################################################
### code chunk number 14: MIHD.Rnw:109-110
###################################################
print(pool(fit))


###################################################
### code chunk number 15: MIHD.Rnw:113-114
###################################################
round(summary(pool(fit)),3)


