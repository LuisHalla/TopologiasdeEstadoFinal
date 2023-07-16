
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TÉCNICAS DE APRENDIZAJE SUPERVISADO APLICADAS A   #
#     CLASIFICACIÓN DE TOPOLOGÍAS DE ESTADO FINAL     #
#                 EN EL EXPERIMENTO SBND              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### PACKAGES ####

# Principales
library(caret) #KNN
library(MASS) #LDA QDA
library(ggplot2)

# Gráficos
library(ggExtra)
library(GGally)
library(ggpubr)
library(cvms) #Matriz de confusión

# Otros
library(dplyr)
library(tidyr)
library(MLmetrics)
library(modeest)


#### DATA ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# TOPOLOGÍAS:                              # 
#                                          #
# Sucesos sin piones                       #
# Sucesos con un (y solo un) pion neutro   #
# sucesos con un (y solo un) pion cargado  #
# Sucesos con más de un pion               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Lectura del fichero de datos
data <- read.csv('data.csv', header = T)

# Se definen las 4 categorías a partir de las columnas nPi0 y nPiCh
levels = c('NoPi', 'Pi0', 'PiCh', 'Other')

topos <- list(data$nPi0==0 & data$nPiCh==0,
              data$nPi0==1 & data$nPiCh==0,
              data$nPi0==0 & data$nPiCh==1,
              data$nPi0 + data$nPiCh > 1)

# Se añade columna con las categorías
Cat <- vector()
for(i in 1:4) Cat[topos[[i]]] <- levels[i]

data <- data.frame(data, 'Cat' = factor(Cat, levels))


# Divisón en conjuntos de entrenamiento (80%) y de test (20%)
index <- createDataPartition(data$Cat, p = 0.8, list = FALSE)
train <- data[index, ]
test <- data[-index, ]


# Creación subsets para análisis exploratorio
data2t <- subset(data, subset = Cat == 'NoPi' | Cat == 'Other')
data2t$Cat <- factor(data2t$Cat, labels = c('NoPi', 'Other'))

# Se crea una nueva columna para almacenar las predicciones
data$Cat_pred <- data$Cat
data2t$Cat_pred <- data2t$Cat


# Para mejor visualización de maxOpenAngle y maxShowerLength, se eliminan los
# sucesos sin cascadas al representar (todos los sucesos sin cascadas tienen
# ambas variables igual a cero)
mOA <- subset(data, subset = nShowers!=0)


# Creación subsets para logística
data12 <- subset(data, subset = Cat == 'NoPi' | Cat == 'Pi0')
train12 <- data12[index,]
test12 <- data12[-index,]

data13 <- subset(data, subset = Cat == 'NoPi' | Cat == 'PiCh')
train13 <- data13[index,]
test13 <- data13[-index,]

data14 <- subset(data, subset = Cat == 'NoPi' | Cat == 'Other')
train14 <- data14[index,]
test14 <- data14[-index,]


# Creación subsets para KNN
data_knn <- data
data_knn[5:10] <- scale(data_knn[5:10], center = TRUE, scale = TRUE)

train_knn <- data_knn[index, ]
test_knn <- data_knn[-index, ]

#### FUNCIONES ####

# Crear matriz de confusión
cm <- function(targets, predictions, positive=2){
  
  confusion_matrix(targets=targets,
                   predictions=predictions,
                   positive=positive,
                   metrics=list("Accuracy" = TRUE))
}


# Mostrar métricas
cm_metrics <- function(cm, binary=T){
  
  Acc <- ifelse(binary, cm$Accuracy, cm$`Overall Accuracy`)
  
  list('Accuracy'=Acc,
       'Sensitivity'=cm$Sensitivity,
       'Specificity'=cm$Specificity,
       'PPV'=cm$`Pos Pred Value`,
       'NPV'=cm$`Neg Pred Value`,
       'F1'=cm$F1,
       'Kappa'=cm$Kappa)
}


# Mostrar matriz
cm_plot <- function(cm, binary=T, class='Other'){
  
  if (binary) class_order <- c(class, 'NoPi')
  else class_order <- c('Other', 'PiCh', 'Pi0', 'NoPi')
  
  plot_confusion_matrix(cm,
                        class_order=class_order,
                        diag_percentages_only=T,
                        add_counts=T,
                        add_normalized=F,
                        palette="Greens")
}


#### EXPLORATORIO: EVENTOS Y ESTADÍSTICOS ####

# Eventos
total <- nrow(data)
table <- table(data$Cat)

porcentaje <- 100 * table / total

# OUTPUT
total
table
porcentaje


# Estadísticos
vars <- c('nShowers', 'maxShowerLength', 'maxOpenAngle',
          'nTracks', 'maxTrackLength', 'nStubs')

# Bucle que calcula los parámetros para cada una de las variables reconstruidas
# y los almacena en un dataframe. Impone NA a la modas de las variables continuas
estadisticos <- list()

for (i in 5:10) {
  estadisticos[[i-4]] <- data.frame('Media' = tapply(data[[i]], data$Cat, mean),
                                    'Mediana' = tapply(data[[i]], data$Cat, median),
                                    'SD' = tapply(data[[i]], data$Cat, sd),
                                    'Moda' = tapply(data[[i]], data$Cat, mfv))
  if (i %in% c(6,7,9)) estadisticos[[i-4]]$Moda <- NA
}
names(estadisticos) <- vars

# OUTPUT
estadisticos


#### EXPLORATORIO: GRÁFICOS DE CAJAS ####

cajas <- list()

cajas[[1]] <- ggplot(data, aes(x=Cat, y=nShowers, color=Cat)) +
  geom_boxplot() + ylim(0,15)

cajas[[2]] <- ggplot(data, aes(x=Cat, y=maxShowerLength, color=Cat)) +
  geom_boxplot() + ylim(0,200)

cajas[[3]] <- ggplot(data, aes(x=Cat, y=nTracks, color=Cat)) +
  geom_boxplot()

cajas[[4]] <- ggplot(data, aes(x=Cat, y=maxOpenAngle, color=Cat)) +
  geom_boxplot()

figure <- ggarrange(cajas[[1]], cajas[[2]],
                    cajas[[3]], cajas[[4]],
                    ncol=2, nrow=2,
                    common.legend=T, legend="top")

# OUTPUT
figure


#### EXPLORATORIO: HISTOGRAMAS Y DENSIDADES ####

overlap <- list()

overlap[[1]] <- ggplot(data, aes(nShowers, fill=Cat)) +
                geom_histogram(binwidth=1) + xlim(0,10)
overlap[[2]] <- ggplot(data, aes(nTracks, fill=Cat)) +
                geom_histogram(binwidth=1) + xlim(0,10)
overlap[[3]] <- ggplot(data, aes(nStubs, fill=Cat)) +
                geom_histogram(binwidth=1) + xlim(0,14)

overlap[[4]] <- ggplot(mOA, aes(maxShowerLength, fill=Cat, colour=Cat)) +
                geom_density(position="stack") + xlim(0,175)
overlap[[5]] <- ggplot(data, aes(maxTrackLength, fill=Cat, colour=Cat)) +
                geom_density(position="stack") + xlim(0,500)
overlap[[6]] <- ggplot(mOA, aes(maxOpenAngle, fill=Cat, colour=Cat)) +
                geom_density(position="stack")

figure <- ggarrange(overlap[[1]], overlap[[4]],
                    overlap[[2]], overlap[[5]],
                    overlap[[3]], overlap[[6]],
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "top")

# Energía (GeV)
total_E <- ggplot(data, aes(trueE)) + geom_density(fill="black")
cat_E <- ggplot(data, aes(trueE, fill=Cat)) + geom_density()

# OUTPUT
figure
total_E
cat_E


#### CLASIFICACIÓN PRIMITIVA: ELECCIÓN DE VARIABLES ####

dispersion <- ggpairs(data2t,
              columns = c(5,8,10),
              aes(color=Cat),
              legend = 1,
              upper = "blank",
              diag = list(continuous = wrap("barDiag",
                                            binwidth = 1,
                                            color="black")))

Showers_Tracks <- ggMarginal(ggplot(data2t,
                                    aes(x=nShowers, y=nTracks, color=Cat)) +
                               geom_point(size=2) +
                               xlim (-2,15) + ylim (-2,15) +
                               stat_ellipse(type = "norm", level = 0.95),
                             type="histogram",
                             binwidth = 1, groupFill=T, alpha=1)

# OUTPUT
dispersion
Showers_Tracks


#### CLASIFICACIÓN PRIMITIVA: 2 VARIABLES Y 2 TOPOLOGÍAS ####

max_x <- max(data2t$nShowers)
max_y <- max(data2t$nTracks)
table2t <- table(data2t$Cat)

for(x in 0:max_x) for(y in 0:max_y){
  
  data2t_xy <- data2t$nShowers == x & data2t$nTracks == y
  
  table_xy <- table(data2t[data2t_xy,]$Cat)
  max <- which.max(table_xy / table2t)
  
  if (sum(table_xy)>0) data2t[data2t_xy,]$Cat_pred <- factor(names(max))
}

conf_mat <- cm(targets=data2t$Cat,
               predictions=data2t$Cat_pred,
               positive='NoPi')

conf_mat_metrics <- cm_metrics(conf_mat)
conf_mat_plot <- cm_plot(conf_mat)

# OUTPUT
conf_mat_metrics
conf_mat_plot


#### CLASIFICACIÓN PRIMITIVA: 2 VARIABLES Y 4 TOPOLOGÍAS ####

max_x <- max(data$nShowers)
max_y <- max(data$nTracks)
table <- table(data$Cat)

for(x in 0:max_x) for(y in 0:max_y){
  
  data_xy <- data$nShowers == x & data$nTracks == y
  
  table_xy <- table(data[data_xy,]$Cat)
  max <- which.max(table_xy / table)
  
  if (sum(table_xy)>0) data[data_xy,]$Cat_pred <- factor(names(max))
}

conf_mat <- cm(targets=data$Cat,
               predictions=data$Cat_pred)

conf_mat_metrics <- cm_metrics(conf_mat, binary=F)
conf_mat_plot <- cm_plot(conf_mat, binary=F)

# OUTPUT
conf_mat_metrics
conf_mat_plot


#### APRENDIZAJE SUPERVISADO: LDA ####

lda_all <- lda(Cat ~ nShowers + nTracks + nStubs +
                     maxShowerLength + maxTrackLength + maxOpenAngle,
               data = train)

lda_pred <- predict(lda_all, test)

conf_mat <- cm(targets=test$Cat,
               predictions=lda_pred$class)

conf_mat_metrics <- cm_metrics(conf_mat, binary=F)
conf_mat_plot <- cm_plot(conf_mat, binary=F)

# OUTPUT
conf_mat_metrics
conf_mat_plot


#### APRENDIZAJE SUPERVISADO: QDA ####

qda_all <- qda(Cat ~ nShowers + nTracks + nStubs +
                     maxShowerLength + maxTrackLength + maxOpenAngle,
               data = train)

qda_pred <- predict(qda_all, test)

conf_mat <- cm(targets=test$Cat,
               predictions=qda_pred$class)

conf_mat_metrics <- cm_metrics(conf_mat, binary=F)
conf_mat_plot <- cm_plot(conf_mat, binary=F)

# OUTPUT
conf_mat_metrics
conf_mat_plot


#### APRENDIZAJE SUPERVISADO: LOGÍSTICA NoPi vs Pi0 ####

# Model selection
logit12 <- glm(Cat ~ nShowers + nTracks + nStubs +
                 maxShowerLength + maxTrackLength + maxOpenAngle,
               data = train12, family = binomial)

fstep12 <- step(logit12, trace=F)

# Training
logit12 <- glm(Cat ~ nShowers + nTracks + nStubs +
                 maxShowerLength,
               data = train12, family = binomial)

# Test
pred <- predict(logit12, test12, type='response') < .5

test12$Cat_pred[pred] <- 'NoPi'
test12$Cat_pred[!pred] <- 'Pi0'

conf_mat <- cm(targets=factor(test12$Cat),
               predictions=factor(test12$Cat_pred))

conf_mat_metrics <- cm_metrics(conf_mat)
conf_mat_plot <- cm_plot(conf_mat, class='Pi0')

# OUTPUT
summary(fstep12)
summary(logit12)

conf_mat_metrics
conf_mat_plot


#### APRENDIZAJE SUPERVISADO: LOGÍSTICA NoPi vs PiCh ####

# Model selection
logit13 <- glm(Cat ~ nShowers + nTracks + nStubs +
                 maxShowerLength + maxTrackLength + maxOpenAngle,
               data = train13, family = binomial)

fstep13 <- step(logit13, trace=F)

# Training (todas las variables permanecen en el modelo)

# Test
pred <- predict(logit13, test13, type='response') < .5

test13$Cat_pred[pred] <- 'NoPi'
test13$Cat_pred[!pred] <- 'PiCh'

conf_mat <- cm(targets=factor(test13$Cat),
               predictions=factor(test13$Cat_pred))

conf_mat_metrics <- cm_metrics(conf_mat)
conf_mat_plot <- cm_plot(conf_mat, class='PiCh')

# OUTPUT
summary(fstep13)
summary(logit13)

conf_mat_metrics
conf_mat_plot


#### APRENDIZAJE SUPERVISADO: LOGÍSTICA NoPi vs Other ####

# Model selection
logit14 <- glm(Cat ~ nShowers + nTracks + nStubs +
                 maxShowerLength + maxTrackLength + maxOpenAngle,
               data = train14, family = binomial)

fstep14 <- step(logit14, trace=F)

# Training
logit14 <- glm(Cat ~ nShowers + nTracks + nStubs +
                 maxShowerLength  + maxOpenAngle,
               data = train14, family = binomial)

# Test
pred <- predict(logit14, test14, type='response') < .5

test14$Cat_pred[pred] <- 'NoPi'
test14$Cat_pred[!pred] <- 'Other'

conf_mat <- cm(targets=factor(test14$Cat),
               predictions=factor(test14$Cat_pred))

conf_mat_metrics <- cm_metrics(conf_mat)
conf_mat_plot <- cm_plot(conf_mat, class='Other')

# OUTPUT
summary(fstep14)
summary(logit14)

conf_mat_metrics
conf_mat_plot

#### APRENDIZAJE SUPERVISADO: KNN ####

# Training: hyper-parameter CV tuning
training_control <- trainControl(method="cv", number=10)
grid <- data.frame(k = seq(3, 101, by=2))

knn_cv <- train(Cat ~ nShowers + nTracks + nStubs +
                  maxTrackLength + maxShowerLength + maxOpenAngle, 
                data=train_knn,
                method="knn",
                trControl=training_control,
                metric="Kappa",
                tuneGrid=grid)

# Test
test_knn$Cat_pred <- predict(knn_cv, newdata=test_knn, type="raw")

conf_mat <- cm(targets=test_knn$Cat,
               predictions=test_knn$Cat_pred)

conf_mat_metrics <- cm_metrics(conf_mat, binary=F)
conf_mat_plot <- cm_plot(conf_mat, binary=F)

# OUTPUT
knn_cv$bestTune

conf_mat_metrics
conf_mat_plot
