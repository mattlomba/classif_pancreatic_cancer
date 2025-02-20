#Mattia Lombardi#

#project of Computational Statistics:
#Model based Classification of pancreatic cancer

#libraries:
library(tidyverse)
library(ggplot2)
library(mclust)
library(visdat) #missed data visualization
library(ggcorrplot)
library(GGally) 
library(ggthemes)
library(factoextra) #pca library
library(caret)

#-----------------------------------------
#PCA and Exploration Data Analysis (EDA)
#-----------------------------------------

setwd("C:/Users/39349/OneDrive/Desktop/UNIMIB/ALTRO/dataset")

data<-read.csv("pancreatic_cancer.csv")
#visualize data
str(data)
View(data)

#cambio in numerica la variabile sex M-F
data$sex<-ifelse(data$sex=="M",1,2)
table(data$diagnosis)/590
table(data$sex)/590
#proportion are ok 

#missing data
anyNA(data) #TRUE
vis_miss(data)

#purtroppo le variabili REG1A e plasma_CA19_9
#hanno molti valori mancanti
#per fare un analisi preliminare delle PCA 
#consideriamo solamente le 209 osservazioni complete

data_comp<-data[complete.cases(data),]
#considero variabili solo numeriche
data_num<-data_comp[,-c(1,2,3,6,7,8)]
data_lab<-as.factor(data_comp[,6])
table(data_lab)/209

#possiamo notare che i missing values provenivano
#prevalentemente dal gruppo di controllo o con diagnosi benigna
#dati mancanti NMAR
str(data_num)

#correlazione
ggcorrplot(cor(data_num), hc.order = TRUE, type = "lower",lab = TRUE,
           colors = c("green3", "white", "darkorange2"))

#PCA
pca<-princomp(data_num,cor=T)
pca

#rappr. della % di varianza spiegata dalle components
fviz_eig(pca, barfill = "sandybrown",barcolor="black",main = "Principal Component Analysis")
sum(pca$sdev[1:6]^2)/8
#le prime 6 components spiegano 91% varianza within

#grafici non presenti nel report: 
#1)rappr. del grafico rappresentante le 2 componenti principali suddivise per sesso
groups <- as.factor(ifelse(data_num$sex==1,"M","F"))

fviz_pca_ind(pca,
             col.ind = groups, 
             palette = c("#FC4E07","#00AFBB" ),
             addEllipses = TRUE,
             legend.title = "Sex",
             repel = TRUE,
             title=c("2Dim PCA-Sex")
)

#2)rappresentazione grafica 2PCA dimensions - creatinine
fviz_pca_ind(pca, col.ind=data_num$creatinine,
             legend.title = "Creatinine",
             repel = TRUE,
             title=c("2Dim PCA-Creatinine")
)+scale_color_gradient2(low="white", mid="slateblue",
                        high="orange", midpoint=0.6)

#costruisco il grafico dei loadings
scores<-pca$loadings[,1:6]^2
scores
data_scores<-as.vector(scores)

comp<- c ( rep("Comp.1",8),rep("Comp.2",8),rep("Comp.3",8),rep("Comp.4",8),rep("Comp.5",8),rep("Comp.6",8))
variables<-c("age","sex","plasma_CA19_9","creatinine","LYVE1","REG1B","TFF1","REG1A")
data_loadings<-data.frame(comp,variables,data_scores)
data_loadings

#grafico loadings
ggplot(data_loadings, aes(fill=variables, y=data_scores, x=comp)) + 
  geom_bar(position="stack", stat="identity")+
  labs(title = "Loadings of different components",y="loadings")+
  scale_fill_brewer(palette = "Paired")+
  theme_few()

#data_num
#mostro le funz di densit? delle variabili
#suddivise per ogni componente
#install.packages("caret")
library(caret)
levels(data_lab)<-c("control","benign","malignant")

featurePlot(data_num[,-c(3,8)],data_lab,plot="density",
            scales=list(x = list(relation="free"),y = list(relation="free")),
            adjust = 1.5,pch="|",auto.key = list(columns = 3))

#-----------------
#EDDA CLASSIFIER
#-----------------

#rileggo i dati per evitare eventuali errori
data<-read.csv("pancreatic_cancer.csv")


#data manipulation:
data$sex<-ifelse(data$sex=="M",1,2)

#modifico le etichette in diagnosi senza metastasi e con cancro
data$diagnosis<-ifelse(data$diagnosis %in% c(1,2),1,2)

#considero solo le variabili complete e numeriche
dada<-data[,-c(1,2,3,6,7,8,9,14)]
#etichette reali
data_lab<-as.factor(data$diagnosis)

#per ottenere sempre stessi risultati
set.seed(123)

#training e test set
test_lab = sample(1: 590 ,trunc(0.25* 590),replace = F)

train<-as.data.frame(dada[-test_lab,])
data_train_lab<-data_lab[-test_lab]

#dataframe per inserire risultati
results = as.data.frame ( matrix ( nrow = 500 , ncol = 2))
colnames(results)=c("model","cv")

library(Rmixmod)
clas.train=mixmodLearn(train,data_train_lab,
                       models = mixmodGaussianModel(family="all",equal.proportions=FALSE),
                       criterion="BIC")
clas.train

for ( i in 1:500){
  clas.train=mixmodLearn(train,data_train_lab,
                         models = mixmodGaussianModel(family="all",equal.proportions=FALSE),
                         criterion=c("CV","BIC"))
  
  model=clas.train@bestResult@model
  cv=clas.train@bestResult@criterionValue[1]
  results[i,]=c(model,cv)
}

#ordino i risultati
results$cv = as.numeric ( results$cv )
results<-results[order(results$cv,method="radix"),]

#seleziono miglior modello
best.mod = mixmodLearn (train,data_train_lab,
                        models = mixmodGaussianModel(listModels=results$model[1]))
best.mod["bestResult"]

#prediction
pred.values=mixmodPredict(data=dada[test_lab,],classificationRule=best.mod["bestResult"]) 

Uncertainty=apply(pred.values@proba,1,min)
types=factor(pred.values@partition,levels =c(1,2),labels = c("B","M"))

act = factor ( data_lab [ test_lab ] , levels = c ("1","2") ,
               labels = c (" B ","M"))

prd = factor ( pred.values@partition , levels = c ("1","2") ,
               labels = c (" B ","M"))


library(caret)
conf.matrix = confusionMatrix ( data = act , reference = prd )
conf.matrix

tab_conf<-as_tibble(conf.matrix$table)

#library per costruire una confuzion matrix carina
library(cvms)

matrix<-table(tibble("target"=act,
                     "prediction"=prd))

plot_confusion_matrix(as_tibble(matrix),
                      target_col = "target", 
                      prediction_col = "prediction",
                      counts_col = "n",
                      add_normalized = FALSE,
                      add_col_percentages = FALSE,
                      add_row_percentages = FALSE,
                      palette = "Greens")

#(ADI)
adjustedRandIndex (act,prd)

#visualizzo grafico con unit? misclassificate e livello incertezza
library(mclust)
mis<-classError(act,prd)

#creo un dataset contenente variabili misclassified
error_data<-data[test_lab[mis$misclassified],]
error_data<-cbind(error_data,prd[mis$misclassified])
#ppca graphic of predicted variables
data_pred<-dada[test_lab,]

pca<-princomp(data_pred,cor=T)
prd = factor ( pred.values@partition , levels = c ("1","2") ,
               labels = c (" B ","M"))

miss_class<-mis$misclassified
data_pred$diagnosis<-as.character(prd)
data_pred$diagnosis[miss_class]<-"ErrorClass"
data_pred$diagnosis<-as.factor(data_pred$diagnosis)

#grafico tramite factoextra che visualizza le prime due PCA
fviz_pca_ind(pca, col.ind=data_pred$diagnosis,
             geom = c("point"),alpha=0.8,
             pointsize=Uncertainty,
             legend.title = "Diagnosis",
             repel = TRUE,
             title=c("2Dim PCA-EDDA Classifier"))+
  scale_color_manual(values = c("palegreen2","gray","sandybrown"))

#----------------------------------------------
#FINITE MIXTURE of REGRESSION MODEL (GAUSSIAN)
#----------------------------------------------

#rileggo i dati
data<-read.csv("pancreatic_cancer.csv")

#data manipulation
data$sex<-ifelse(data$sex=="M",1,2)
#in questo caso consideriamo k =3

#PPCA imputation method
library(missMDA)
#ppca funz solo con num variabili
#ipotizza normalit? delle variabili

data_num<-data[-c(1,2,3,6,7,8)]
nb<-estim_ncpPCA(data_num,scale=T)
comp<-imputePCA(data_num,ncp = 5,scale=T)

data.pca<-abs(as.data.frame(comp$completeObs))
#View(data.pca)

data_num2<-data[,-c(1,2,3,7,8,9,14)]

#data visualization (non presente nel report)
#la variabilit? aumenta man mano che la diagnosi peggiora
ggplot(data=data_num2, mapping=aes(x=creatinine,y=LYVE1))+
  geom_point(aes(colour=diagnosis))+
  geom_smooth(method="lm", se=F, size=1.5)+
  theme_light()

library(mixtools)
library(flexmix)
# flexmix(formula, data, k = NULL, cluster = NULL, 
#         model = NULL, concomitant = NULL,.....)

#aggiungo la variabile diagnosis
#nostra variabile risposta della regressione
datapi<-cbind(data[c(6)],data.pca)

#rimuovo outliers
data_ben=datapi[datapi$diagnosis==c(2),]

(outliers_tf<-boxplot.stats(data_ben$TFF1)$out)
(outliers_pla<-boxplot.stats(data_ben$plasma_CA19_9)$out)

data_fin<-datapi[-which(data_ben$TFF1 %in% outliers_tf),]
data_fin2<-data_fin[-which(data_ben$plasma_CA19_9 %in% outliers_pla),]
data_fin3<-data_fin2[-which(data_ben$plasma_CA19_9 %in% outliers_pla),]

datapi<-data_fin3
#str(datapi)
#algoritmo controllo E-M
itermax <- 1000

#entropia
eicval=1
eics<-matrix(nrow=itermax,ncol=2)
#bic
bicval <- Inf
bics<-matrix(nrow=itermax,ncol=2)
#n.b: da considerare abs(), potrebbe assumere valori negativi 

#meccanismo di controllo per EM
#(tempistiche: circa 5 min)

#per ridurre il tempo ho confrontato solamente entropia 
#tramite un metodo di MC, non molto ufficiale, 
#se EIC si discosta molto dalla soglia 0.6 la classificazione 
#inizia ad accorpare le due componenti M-B 

for (iter in 1: itermax)
{
  fit <-flexmix(diagnosis~.,data=datapi,k=3)
  eics[iter,]<-c(iter,EIC(fit))
  if (EIC(fit)<eicval)
  {
    eicval<-EIC(fit)
    bestfit <- fit
  }
}

#results:
summary(bestfit)
plot(bestfit)
BIC(bestfit)
EIC(bestfit)
ICL(bestfit)  #ICL=entropia+BIC
KLdiv(bestfit) #distance Kullback Leiber

fit=bestfit

lab<-as.factor(datapi[,1])
pred<-as.factor(fit@cluster)

#riaggiusto i livelli etichette
levels(pred)<-c(2,3,1)
#pred
#lab
levels(pred)<-c("B","M","C")
levels(lab)<-c("C","B","M")

conf.matrix = confusionMatrix ( data = lab, reference = pred)
conf.matrix
#viasualizza graficamente confuzion matrix
tab_conf<-as_tibble(conf.matrix$table)

matrix<-table(tibble("target"=lab,
                     "prediction"=pred))

plot_confusion_matrix(as_tibble(matrix),
                      target_col = "target", 
                      prediction_col = "prediction",
                      counts_col = "n",
                      add_normalized = FALSE,
                      add_col_percentages = FALSE,
                      add_row_percentages = FALSE,
                      palette = "Oranges")

adjustedRandIndex(lab,pred)
fit@components
#varianza delle components
#come si notava dal grafico iniziale 
#sono pi? elevate man mano che la diagnosi ? pi? grave
(sigma<-parameters(fit)[10,])

#analisi della regressione
summary(refit(fit))

fit@components
fit@df #gradi libert?
#parametri
par<-parameters(besttfit)

#creo dataset con dati predicted
data_pred<-cbind(datapi[,-1],pred)
levels(pred)<-c("B","M","C")

#visualizzo funzione di densit? per etichette previste
#considero LIVE1 in quantouna delle variabili pi? significative 
ggplot(data_pred,aes(x=LYVE1,color=pred))+
  geom_density(lwd=1)+
  xlim(-10,30)+
  scale_color_manual(values=c("green3", "sandybrown", "steelblue"))+
  labs ( title = " FMM density - Pred")+
  theme_bw()


datapi$diagnosis<-as.factor(datapi$diagnosis)
levels(datapi$diagnosis)<-c("C","B","M")
#visualizzo distribuzione dei gruppi reali
ggplot(datapi,aes(x=LYVE1,color=diagnosis))+
  geom_density(lwd=1)+
  xlim(-10,30)+
  theme_bw()+
  scale_color_manual(values = c("steelblue","palegreen2","sandybrown"))+
  labs ( title = " FMM density - Real")

