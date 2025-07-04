# Análisis de datos 

## Normalización

En los estudios metagenómicos, la normalización de datos es un paso fundamental para garantizar comparaciones biológicas precisas entre muestras. 
Esto se debe a que el número total de secuencias obtenidas por muestra puede variar considerablemente debido a factores técnicos, como diferencias en la eficiencia de preparación de librerías, profundidad de secuenciación o incluso sesgos introducidos durante el proceso de PCR (Leek et al., 2010). 
Estas variaciones no están relacionadas con diferencias biológicas reales, pero pueden distorsionar análisis posteriores, especialmente al estimar la abundancia relativa de microorganismos dando como resultado un efecto de lote.
Este efecto suele ser responsable de la mayor fuente de expresión diferencial al combinar datos. Puede enmascarar cualquier diferencia biológica real y llevar a conclusiones erróneas (Leek et al., 2010; Ritchie et al., 2015).
Por ello, se aplican métodos de normalización que permiten ajustar las cuentas de secuencias a una escala comparable. Este proceso asegura que las diferencias observadas reflejen cambios biológicamente relevantes y no artefactos técnicos.

### Escenario inicial: Datos sin normalizar
Problema: La profundidad total de secuenciación es diferente entre muestras (10k, 18k y 27k reads). Esto puede hacer parecer que Microbio 1 es más abundante en la muestra C solo por artefacto técnico.
```
Muestra A       Muestra B       Muestra C
[Microbio 1]     [Microbio 1]     [Microbio 1]
   5000             8000            12000
[Microbio 2]     [Microbio 2]     [Microbio 2]
   3000             6000             9000
[Microbio 3]     [Microbio 3]     [Microbio 3]
   2000             4000             6000

Abundancia bruta (reads)
|
| Microbio 1      ████        ▓▓▓▓▓▓▓▓       ░░░░░░░░░░░░
|                 ████        ▓▓▓▓▓▓▓▓       ░░░░░░░░░░░░
|                 ████        ▓▓▓▓▓▓▓▓       ░░░░░░░░░░░░
| Microbio 2      ██            ▓▓▓▓           ░░░░░░░░░
|                 ██            ▓▓▓▓           ░░░░░░░░░
| Microbio 3      ██            ▓▓▓▓           ░░░░░░░░░
|
+--------------------------------------------------------
                   Muestra A     Muestra B      Muestra C
```
Interpretación errónea posible: "Microbio 1 parece más abundante en la muestra C simplemente porque se secuenció más".



### Aplicar normalización (por ejemplo, proporción relativa)
Ahora las muestras son comparables: se eliminó el sesgo de profundidad.
```
Muestra A       Muestra B       Muestra C
[Microbio 1]     [Microbio 1]     [Microbio 1]
   50%              44.4%           44.4%
[Microbio 2]     [Microbio 2]     [Microbio 2]
   30%              33.3%           33.3%
[Microbio 3]     [Microbio 3]     [Microbio 3]
   20%              22.2%           22.2%


Abundancia relativa (%)
|
| Microbio 1      ████████     ████████      ████████
| Microbio 2      ██████       ███████       ███████
| Microbio 3      ████         ████          ████
|
+--------------------------------------------------------
                   Muestra A     Muestra B      Muestra C
```
La proporción de cada microbio es muy similar entre muestras. Las diferencias previas eran artefactos técnicos, no biológicos.


## Métodos de normalización

Hay tres etapas principales de normalización de ARN-seq que debes considerar:

**1. Normalización dentro de una muestra ("within-sample normalization")**
Este tipo de normalización ajusta los datos dentro de una sola muestra, permitiendo comparar la expresión relativa de genes entre sí en esa muestra. Su objetivo es corregir:
La longitud del gen: Genes más largos tienden a captar más lecturas simplemente por su tamaño, no necesariamente por estar más expresados (Mortazavi et al., 2008).
La profundidad de secuenciación: Dos muestras con diferente número de lecturas no pueden compararse sin ajustar este sesgo técnico.
Métricas como RPKM, FPKM y TPM realizan esta corrección. Sin embargo, estas medidas no son adecuadas para comparar entre muestras, porque no controlan completamente el efecto de la composición total del transcriptoma (Zhao et al., 2021).

**2. Normalización entre muestras dentro de un dataset ("between-sample normalization")**
Este segundo nivel busca que las muestras dentro de un mismo experimento puedan compararse entre sí, eliminando efectos como:

* Diferente número de lecturas por muestra
* Diferencias en la composición global del transcriptoma o genoma
* Recordemos que RNA-seq mide abundancia relativa, no absoluta. Esto significa que si en una muestra hay un grupo de genes muy sobreexpresados, pueden "ocultar" o distorsionar el nivel relativo de otros genes, aunque no hayan cambiado (Robinson & Oshlack, 2010).

Métodos como TMM (Trimmed Mean of M-values), DESeq2's median-of-ratios, quantile normalization, o Total Sum Scaling (TSS), se usan para resolver estos problemas (Bolstad et al., 2003; Love et al., 2014).

**3. Normalización entre datasets ("cross-study normalization" o batch correction)**
Cuando se combinan datos de distintos estudios (por ejemplo, proyectos secuenciados en distintos laboratorios o años), aparecen efectos de lote (batch effects). Estos efectos pueden deberse a:

* El centro o plataforma de secuenciación
* El protocolo experimental
* El momento en que se realizó el experimento

Estos efectos suelen ser la mayor fuente de variabilidad en estudios integrativos, incluso más que las diferencias biológicas reales (Leek et al., 2010). Métodos como:
Combat (Johnson et al., 2007), SVA (Surrogate Variable Analysis; Leek et al., 2012), RUV (Remove Unwanted Variation), se aplican para corregir efectos conocidos y desconocidos en estudios multi-dataset.

-- Métodos más comúnes:

## *Dentro de muestras*

**Counts Per Million (CPM):** es un método de normalización de datos de conteo que ajusta los valores para tener en cuenta la profundidad total de secuenciación.
Normaliza los datos de ARN-seq según la profundidad de secuenciación, pero no según la longitud del gen. Por lo tanto, aunque se trata de un método de normalización intramuestral, la normalización por CPM no es adecuada para comparaciones intramuestrales de la expresión génica.
Se pueden realizar comparaciones entre muestras cuando se utiliza CPM junto con métodos de normalización intraconjunto de datos.
Para cada muestra, se calcula el total de reads. Cada valor en la tabla se divide entre ese total y luego se multiplica por 1 millón, para escalar los valores a un mismo nivel.

Limitaciones
* No corrige por longitud génica: Esto lo diferencia de métodos como RPKM o TPM, que sí lo hacen.
* Puede estar sesgado por genes muy expresos: Si uno o pocos genes dominan la muestra, pueden distorsionar el resto.
* No es adecuado para análisis diferencial riguroso

**Fragments per kilobase of transcript per million fragments mapped (FPKM) para datos de extremos emparejados y Reads per kilobase of transcript per million reads mapped (RPKM) para datos de un solo extremo:**
Mide la abundancia de un transcrito corrigiendo por la longitud del gen o transcrito (en kilobases) y la profundidad total de reads mapeados (en millones). FPKM es similar a RPKM pero En lugar de contar reads, se cuentan fragmentos únicos (Mortazavi et al., 2008).
Las unidades FPKM/RPKM comparan mejor la expresión génica dentro de una misma muestra (Zhao, Ye y Stanton, 2020) y no entre muestras.

Limitaciones
* No controla bien la variabilidad biológica. No es adecuado para análisis estadístico riguroso de diferencias de expresión.
* Sesgo por genes muy expresos puede afectar la normalización global.
* No es recomendado para análisis diferencial avanzado.

**Transcripts Per Million (TPM):**  Al igual que RPKM/FPKM, corrige por longitud del gen/transcrito y profundidad total de secuenciación. Pero lo hace en un orden diferente, lo cual mejora la estabilidad del método cuando se comparan múltiples genes o condiciones.
Divide cada conteo entre la longitud del gen (en kb) y normaliza por el total de reads por kb (escala a millones). 

Limitaciones
* No controla por variabilidad biológica o sesgos técnicos globales: Si hay genes muy expresos en una muestra, pueden afectar la percepción del resto.
* No es adecuado para análisis diferencial estadísticamente robusto.
* No se recomienda si el objetivo es detectar cambios absolutos entre condiciones.

## *Entre muestras*

**Rarefacción:** es una técnica que se usa para igualar el número de lecturas (reads) entre diferentes muestras. La rarefacción consiste en submuestrear aleatoriamente cada muestra para igualar todas a una profundidad de secuenciación común y baja , llamada profundidad de rarefacción.
Si existen tres muestras con valores de 30 kb, 40 kb y 10 kb la rarefacción tomaría solo 10,000 reads de cada una.
La rarefacción se aplica principalmente a datos de conteo, donde cada fila representa una categoría (taxón, gen, OTU, etc.) y cada columna representa una muestra. Algunos ejemplos comunes son:

Limitaciones
* Pérdida de información: descartamos muchos reads, lo cual puede ser un problema si algunas muestras tienen pocos datos inicialmente.
* No ideal para datos muy dispersos (con muchas categorías poco frecuentes).
* No reemplaza a métodos estadísticos más avanzados como normalización por tamaño de biblioteca o modelos basados en distribuciones (ej: Poisson, Negative Binomial).

**Total Sum Scaling (TSS):**  es un método sencillo y ampliamente utilizado para transformar datos de conteo en proporciones relativas entre muestras.
Para cada muestra, se suma todos los conteos (número total de lecturas). Luego, cada conteo individual (por ejemplo, de una especie o ASV/OTU) se divide entre el total de lectas de esa muestra.
El resultado es una proporción relativa, también llamada abundancia relativa.

Limitaciones
* Sesgo por taxa dominantes: Si hay un taxón muy abundante, puede distorsionar la percepción del resto.
* No controla bien el sesgo de bajo recuento: En muestras con pocos reads totales, pequeños cambios pueden afectar mucho las proporciones.
* No es adecuado para análisis estadísticos inferenciales avanzados: No considera la variabilidad biológica ni técnica.

**Trimmed Mean of M-values (TMM):** es un método de normalización global diseñado para corregir diferencias técnicas entre muestras en estudios de RNA-seq o metagenómica funcional.
A diferencia de métodos como RPKM/FPKM o TPM , que son útiles para visualizar tendencias, TMM está pensado para análisis estadísticamente válidos y es el método predilecto en paquetes como edgeR.
Calcula M-values (log-ratio) y A-values (promedio log-intensidad) entre dos muestras, eliminar genes extremadamente altos o bajos (trimming), calcular el promedio ponderado de los M-values restantes y Aplicar el factor de corrección a los conteos originales.

Limitaciones
* No corrige por longitud del gen/transcrito
* Puede ser inestable con pocas muestras o sin replicados.
* No funciona bien con datos muy dispersos o con alto número de genes no expresos.

**Quantile:** La normalización por cuantiles es un método estadístico que hace que todas las distribuciones de expresión (o abundancia) sean comparables entre muestras, forzando a que todas tengan la misma distribución estadística (media, varianza, percentiles, etc.).
Se toma una tabla de datos con genes/filas y muestras/columnas, para cada muestra, se ordenan los valores de menor a mayor. Se calcula el promedio de todos los valores correspondientes a cada posición (es decir, el promedio del primer valor de todas las muestras, luego el segundo, etc.).
Estos promedios se reasignan a cada muestra según su rango original (esto mantiene el orden relativo dentro de cada muestra, pero iguala la distribución global).

Limitaciones
* No corrige por longitud génica ni por profundidad de secuenciación.
* Puede distorsionar diferencias biológicas reales si hay muchos cambios globales entre condiciones.
* No es adecuado para análisis diferencial basado en modelos estadísticos como DESeq2 o edgeR, ya que altera la estructura de varianza-media.
* Solo debe aplicarse a datos transformados (ej: log2), no a conteos brutos.

**DESeq2:** calcula el factor de normalización por profundidad de secuenciación. Utiliza la mediana del cociente por gen (median-of-ratios method). Luego estima la dispersión para cada taxón.
Calcula la media geométrica para cada taxón (fila) entre todas las muestras. Para cada muestra, divide cada conteo por la media geométrica de ese taxón. El factor de normalización de cada muestra es la mediana de esos cocientes (ignorando genes con media geométrica cero).

Esto corrige por diferencias de profundidad de secuenciación y composición global.

## Transformación
En muchos análisis estadísticos tradicionales (como regresión, ANOVA, PCA, etc.), se asumen ciertos supuestos sobre la distribución de los datos. 
Sin embargo, los datos biológicos (especialmente los de abundancias ómicas) rara vez cumplen estos supuestos en su forma cruda. Por ejemplo, los datos de abundancia:

* Son altamente asimétricos (skewed), con muchos ceros y pocos valores muy altos.
* Pueden estar en escalas logarítmicas naturales, por lo que transformarlos ayuda a hacerlos más comparables.
* Tienen varianza que depende de la media (problema típico de datos de conteo como RNA-seq, metagenómica o microbiota).

Esto puede llevar a resultados engañosos, especialmente en análisis multivariante, regresión o tests de hipótesis. La transformación de datos busca corregir estos problemas, permitiendo que los métodos estadísticos comunes (como ANOVA, t-test, PCA o modelos lineales) se apliquen de forma más confiable (Quinn, et al., 2018).
Los principales supuestos que intentamos satisfacer al aplicar una transformación son:

* Normalidad : Muchas pruebas estadísticas asumen que los datos siguen una distribución normal.
* Homocedasticidad : La varianza debe ser similar entre grupos para evitar sesgos en los resultados.
* Linealidad : En modelos predictivos o correlacionales, se espera una relación lineal entre variables.
* Estabilidad de varianza : Especialmente importante en datos donde la varianza aumenta con la media (como conteos de reads).


#### Tipos comunes de transformaciones

* Transformación logarítmica: log(x + 1) para evitar valores negativos o cero. Reduce la asimetría (skewness), comprimir valores extremos, estabilizar la varianza. hace que los datos sean más "lineales".
* Transformación raíz cuadrada: Útil cuando los datos son conteos discretos. Especialmente útil para datos de conteo con muchos valores bajos y pocos altos.
* Hellinger (raíz cuadrada de proporciones): Frecuentemente usada en ecología y microbioma, sobre todo cuando los datos son proporcionales (como luego de TSS). Mejora la aplicabilidad de métodos lineales como PCA o clustering.
* Variance Stabilizing Transformation (VST): Propia de paquetes como DESeq2, esta transformación estabiliza la varianza sin necesidad de eliminar genes con baja expresión. DESeq2 transforma cada valor de conteo con una fórmula no lineal derivada del modelo binomial negativo.

###  Importante: transformar no es lo mismo que normalizar
Normalización ajusta por factores técnicos (profundidad de secuenciación, longitud del gen, etc.). Transformación adapta los datos a los requisitos estadísticos o matemáticos del método que vas a aplicar.
A menudo se hacen ambos pasos, y en ese orden: primero normalizar, luego transformar.


## Ejemplo de normalización, transformación, análisis de diferenciación y PCA
Obtención de los datos
```R
library(tidyverse)

# Leer datos
filtered_data <- read.csv("data.csv", header = T, sep = ",")
filtered_data <- filtered_data[,-1]

# Crear la matriz de abundancia absoluta
head(filtered_data)

abundance_matrix <- filtered_data %>%
  group_by(sample, name) %>%
  summarise(total_percentage = sum(reads, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = total_percentage)
abundance_matrix$name <- trimws(abundance_matrix$name)

# Obtener fenodata
feno <- filtered_data %>%
  select(sample, date) %>%
  distinct() %>%
  mutate(
    date = as.factor(str_remove(date, " ")),
    date2 = as.factor(c(rep("Hibernation", 13), rep("Breeding season", 35)))
  ) %>%
  column_to_rownames("sample")
```
<details> 
   <summary> Resumen general del análisis con DESeq2 </summary>
Tomar matriz de conteos crudos (abundance_matrix)
→ Filas = taxones, columnas = muestras.
Se construye un objeto DESeqDataSet (dds) con diseño experimental (~ date2)

Si se ejecuta DESeq(dds):
 * Normaliza por tamaño de biblioteca (usa size factors). 
 * Estima dispersión de cada taxón (cuánto varía su abundancia entre muestras dentro del mismo grupo).
 * Ajusta modelos negativos binomiales por taxón (los datos de conteo no siguen una distribución normal, sino una binomial negativa).

</details>

Veamos un ejemplo con Deseq:
```R
# Se requiere Biocmanager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Cargar paquetes
library(BiocManager)
library(DESeq2)
library(vsn)
library(hexbin)

### Exploracion
meanSdPlot(as.matrix(abundance_matrix[,2:49]), xlab = "Taxa", ylab = "Desviación estandar")
boxplot(abundance_matrix[,2:49], las = 2, cex.axis = 0.7)

```
Los datos descargados no han sido previamente normalizados ni han sufrido una transformación. En datos de secuenciación crudos, los genes con medias altas suelen tener varianza muy alta, y los genes poco expresados, varianza muy baja.

Se procede a realizar una normalización por Deseq y a ajustar el modelo binomial negativo de la distribución de la abundancia de conteos.
```R
######### Analisis de abundancia diferencial ##########

# se crea un data frame con la matriz de abundancia
abundance_matrix <- as.data.frame(abundance_matrix)
row.names(abundance_matrix) <- abundance_matrix$name 
abundance_matrix <- abundance_matrix[,-1]

# crea objeto deseq
dds <- DESeqDataSetFromMatrix(countData = abundance_matrix,
                              colData = feno,
                              design = ~ date2)  # ~1 no factor modelado
```
Se transforma los datos con VST (varianceStabilizingTransformation): La varianza deja de depender fuertemente de la media, algo que sí ocurre en datos de conteo sin transformar. Se obtiene una matriz con varianza aproximadamente constante lista para PCA, clustering, ANOVA exploratoria, etc.
* blind = TRUE: ignora el diseño experimental al transformar si solo se quiere explorar la variabilidad global. Si se quisiera comparar grupos, blind = FALSE.
* 

```R
# Transformacion vst (Variance Stabilizing Transformation)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean")
rlog_counts <- assay(vsd)

# verifica si existen NA y los nombres de las columnas
sum(is.na(rlog_counts))
colnames(rlog_counts)

# Evalua la estabilización de la varianza
boxplot(rlog_counts, las = 2, cex.axis = 0.7)
meanSdPlot(rlog_counts, xlab = "Taxa", ylab = "Desviación estandar")

```

### Pruebas estadísticas para diferencias entre grupos.

Cuando se recopilan datos para varias variables en las mismas unidades muestrales, siempre es posible examinar las variables una por una en lo que respecta a las pruebas de significancia. Desafortunadamente, este enfoque simple presenta una desventaja debido al uso repetido de pruebas de significancia, cada una de las cuales tiene cierta probabilidad de llevar a una conclusión errónea.

Para evaluar si las varianzas entre grupos son iguales (suposición de homocedasticidad), utilizamos la prueba de Levene.
- **Hipótesis nula ($H_0$):** Las varianzas son iguales entre los grupos.
- **Hipótesis alternativa ($H_1$):** Al menos una varianza es diferente.
  
```R
rlog_counts # matriz
feno # factor

# Prueba t univariante de dos muestras. Para un solo taxón (por ejemplo, la fila 1)
library(car)
leveneTest(rlog_counts[1, ], feno$date2) # prueba de homocedasticidad
t.test(rlog_counts[1, ] ~ feno$date2, var.equal = FALSE)
```
Es posible determinar cuál de estas variables, si alguna, parece haber tenido valores medios diferentes para los grupos de comparación. Sin embargo, además de estas pruebas, también podría ser interesante saber si todas las variables consideradas en conjunto sugieren una diferencia entre los grupos comparados.
Para responder a esta pregunta se necesita una prueba multivariante. Una posibilidad es la prueba T2 de Hotelling. El estadístico utilizado es entonces una generalización del estadístico t. 

Para evaluar si las matrices de covarianza entre grupos son iguales (suposición de homocedasticidad), utilizamos la prueba de Box.
- **Hipótesis nula ($H_0$):** Las matrices de covarianza son iguales entre los grupos.
- **Hipótesis alternativa ($H_1$):** Las matrices de covarianza son distintas entre los grupos.
  
```R
#Box’s M-test (homogeneidad de matrices de covarianza)
library(biotools)
# Requiere que cada nivel del factor tenga n <= p. Asegurarse de que p < n_grupo.
boxM(t(rlog_counts[1:12, ]), feno$date2)

# Hotelling’s T² (multivariado, dos grupos)
# Subset de taxa para cumplir n > p, por ejemplo los 12 primeros

library(Hotelling)
library(corpcor)
hotel_resp <- hotelling.test(t(rlog_counts[1:12, ]) ~ feno$date2)
hotel_resp
```
Casos que involucran la clasificación de observaciones en función de *uno o más factores*, y que incluyen *múltiples variables respuesta*, pueden analizarse mediante una generalización del análisis de varianza conocida como **Análisis Multivariado de Varianza (MANOVA)**. 

```R
# MANOVA (multivariado, dos grupos).
response <- t(rlog_counts[1:13, ]) # Subset para reducir dimensionalidad
manova_model <- manova(response ~ feno$date2)
summary(manova_model, test = "Wilks")
```
Utilizando DESeq2: DESeq2 compara los grupos definidos por date2 (por ejemplo, "Breeding season" vs "Hibernation") usando el modelo ajustado. Para esto se usa una prueba de Wald, que evalúa si el log2 fold change (log2FC) entre los grupos es significativamente distinto de cero. El estadístico de Wald se calcula dividiendo el LFC por su error estándar, y este valor se compara con una distribución normal estándar para obtener un valor p.

```R
# nombres de los coeficientes disponibles en el objeto con los resultados
resultsNames(dds)

# Se extraen los valores de interés que son distintos significativamente para el contraste especificado
res <- results(dds, name = "date2_Breeding_season_vs_Hibernation_season", alpha = 0.05) # Breeding como categoría de referencia y ajuste por FDR en 0.05

# filtro solo los significativos luego del ajuste por FDR
dif_spec <- res[!is.na(res$pvalue) & res$pvalue < 0.05,]
dif_spec
```

### Principal Component Analysis y análisis de cluster
La escala de las variables numéricas es homogénea, y los datos son del mismo tipo (conteos de abundancia) por lo que nos basaremos en la matriz de covarianzas de los datos centrados.

#### Covarianza

La covarianza entre dos variables aleatorias $X$ y $Y$ se define como:

$$
\text{Cov}(X, Y) = \frac{1}{n - 1} \sum_{i=1}^{n} (X_i - \bar{X})(Y_i - \bar{Y})
$$

- $X_i$, $Y_i$: valores individuales de cada variable.

- $Xˉ$, $Yˉ$: medias de las variables X y Y, respectivamente.

- $n$: número total de observaciones.

- Se usa $n−1$ para calcular la covarianza muestral; si es poblacional, usaría $n$.


```R

dat_scaled <- scale(t(rlog_counts), center = T, scale = F) 

```
#### Por qué scale = FALSE? Relación entre covarianza y correlación
A partir de una matriz de covarianzas, podemos obtener la matriz de correlaciones. El coeficiente de correlación entre dos variables $j$ y $k$ se define como:

$$
r_{jk} = \frac{c_{jk}}{s_j s_k}
$$

donde:
- $r_{jk}$: coeficiente de correlación
- $c_{jk}$: covarianza entre las variables $j$ y $k$
- $s_j$, $s_k$: desviaciones estándar de las variables $j$ y $k$

Calculamos la matriz de varianzas ajustando a dividir por n en vez de por (n-1) para que los resultados sean comparables con la matriz de correlación.

* En estadística inferencial, al estimar la covarianza poblacional a partir de una muestra, se usa n−1 para corregir el sesgo (“grados de libertad”). En PCA, sin embargo, no estamos estimando, sino describiendo la estructura interna de esos mismos datos. Queremos la “matriz de momentos de segundo orden” sin corrección de sesgo.

```R
n <- dim(abundance_matrix)[1]
S <-cov(dat_scaled)
S1 <- cov(dat_scaled)*(n-1)/n
show(S)
```

Calculando las componentes principales con funciones básicas

#### mediante la diagonalización de la matriz de varianza

```R
# calculando eigenvalues y eigenvectores
EIG <- eigen(S)
show(EIG)

# multiplicar la matriz original por la matriz de vectores propios
eigenVecs1 <- EIG$vectors
PCAS1 <- dat_scaled %*% eigenVecs1
head(PCAS1)

# Calcular porcentaje de variabilidad explicado.
vars1<- EIG$values/sum(EIG$values)
round(vars1,3)

# Graficar
xlabel <- paste("PCA1 ", round(vars1[1]*100, 2),"%" )
ylabel <- paste("PCA2 ", round(vars1[2]*100,2),"%" )
plot(PCAS1[,1], PCAS1[,2], main = "2 primeras PCs",
     xlab=xlabel, ylab=ylabel)

```

Calculando las componentes principales con funciones de R
```R
library(ggplot2)
library(mclust)

# PCA en datos transformados y normalizados
pca_data <- prcomp(t(rlog_counts))
pca_data$x


### Analisis de clustering (PCA basado en abundancias funcionales)
varianza_explicada <- pca_data$sdev^2 / sum(pca_data$sdev^2)

# Mostrar la varianza explicada
print(varianza_explicada)
head(varianza_explicada, 2)
sum(varianza_explicada)

# Clusters GMM (Gaussian Mixture Models): Modela los datos como una mezcla de distribuciones gaussianas.
# A diferencia del k-means, GMM permite clusters con diferentes tamaños, orientaciones y densidades.

set.seed(123)
# convertir a data.frame
pca_data <- as.data.frame(pca_data$x)

# ajuste del modelo
gmm_result <- Mclust(pca_data[, c("PC1", "PC2")])

# volver factor el grupo de clasificación
pca_data$group <- as.factor(gmm_result$classification)

gmm_result$BIC

# Grafico de PCA 
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) + 
  geom_point(size = 3) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2, linewidth = 1) +
  geom_text(aes(label = rownames(pca_data)), color = "black", hjust = 1.25) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  labs(title = paste("Log-likelihood: ", round(gmm_result$loglik,2),
                     "\nBayesian Information Criterion: ", round(gmm_result$bic),
                     "\n (Next top BIC: VEV 2.24     VVV 1.15     VVE -0.47"), 
       x = "Principal Component 1 (33%)", y = "Principal Component 2 (20.5%)", color = "Group")


# Exportar a 16X 10
```

## Ejemplo de rarefaccion, análisis de distancias y Principal Coordinates Analysis con VEGAN

```R
# cargar paquetes

library(permute)
library(lattice)
library(vegan)
library(boot)

# Formatear matriz

abundance_matrix_2 <- filtered_data %>%
  group_by(sample, name) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = total_reads, values_fill = 0) %>%
  column_to_rownames("name")

dat <- t(abundance_matrix_2)

# Rarefaccion
# Funcion que genera puntos de rarefaccion para cada muestra
rarefaction_df <- function(mat, step = 100, max_depth = NULL) {
  results <- list()
  for (i in 1:nrow(mat)) {
    sample_name <- rownames(mat)[i]
    counts <- mat[i, ]
    counts <- counts[which(counts > 0)]
    total_reads <- sum(counts)
    
    max_reads <- if (is.null(max_depth)) total_reads else min(max_depth, total_reads)
    depths <- seq(1, max_reads, by = step)
    
    richness <- sapply(depths, function(d) {
      rarefy(counts, sample = d)
    })
    
    df <- data.frame(
      Sample = sample_name,
      Depth = depths,
      Richness = richness
    )
    
    results[[i]] <- df
  }
  
  do.call(rbind, results)
}

# Profundidad de corte (probar distintos puntos de corte)
min_depth <- 6000

# Generar el dataframe para ggplot
rarefaction_data <- rarefaction_df(dat, step = 100, max_depth = min_depth)

ggplot(rarefaction_data, aes(x = Depth, y = Richness, color = Sample)) +
  geom_line(size = 1) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(color = "black"),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title = element_text(size = 20, face = "bold")) +
  labs(title = "",
       x = "Number of reads",
       y = "Estimated richness (number of orders)") 

# Guardar a 16X10
```
## Distancias
### Distancia de Mahalanobis
Cuando existen más variables (p) que observaciones (n), el cálculo estándar de la distancia de Mahalanobis falla porque la matriz de covarianza $Σ$ es singular (no invertible), un caso extremadamente común en metagenómica.

$$D²(x) = (x − μ)ᵀ Σ⁻¹ (x − μ)$$

Donde:
- `x`: vector de observaciones
- `μ`: vector de medias
- `Σ`: matriz de covarianza
- `Σ⁻¹`: inversa de la matriz de covarianza

Si p >= n (número de variables >= número de muestras),
entonces Sigma es singular → Sigma⁻¹ no existe.

Se puede solucionar a través de la reducción del número de dimensiones utilizando un subconjunto de taxa (por ejemplo, los más variables o relevantes), o utilizando una estimación regularizada de la matriz de covarianza, que es invertible incluso si p > n.

```R
library(corpcor)

# Transponer: filas = observaciones
rlog_counts_t <- t(rlog_counts)

# Calcular matriz de covarianza regularizada (shrinkage)
cov_shrink <- cov.shrink(rlog_counts_t)

# Media
mu <- colMeans(rlog_counts_t)

# Calcular distancia de Mahalanobis regularizada
mahal_dist <- mahalanobis(rlog_counts_t, center = mu, cov = cov_shrink)

```
### Análisis de distancia para presencias ausencias (diversidad beta)
En metagenómica, una situación común es aquella en la que la similitud o distancia entre dos elementos debe basarse en una lista de sus presencias y ausencias. Existen medidas de similitud de uso común como:


- Índice de coincidencia simple:  
  $$a / (a + d + n)$$

- Índice de Ochiai:  
  $$a / √[(a + b)(a + c)]$$

- Índice de Dice-Sørensen:  
  $$2a / (2a + b + c)$$

- Índice de Jaccard:  
  $$a / (a + b + c)$$

Todos estos varían de cero (sin similitud) a uno (similitud completa), por lo que las medidas de distancia complementarias se pueden calcular restando los índices de similitud de uno.


```R
# Normalizacion por TSS -Total Sum Scaling- y transformacion
tss_matrix <- sweep(abundance_matrix, 1, rowSums(abundance_matrix), "/") # Normalización
hellinger_data <- decostand(tss_matrix, method = "hellinger") # Suaviza las diferencias grandes entre abundancias

# Transformación a presencia/ausencia para índices binarios
pa_matrix <- (abundance_matrix > 0) * 1

# Analisis de diversidad Beta. vegdist ya devuelve 1 - Índice de similitud
bray_dist <- vegdist(t(hellinger_data), method = "bray") # abundancias
jaccard_dist <- vegdist(t(pa_matrix), method = "jaccard") # presencia/ausencia
euclidean_dist <- dist(t(hellinger_data), method = "euclidean") # abundancias normalizadas, muy sensible a diferencias extremas

# a matrices de distancias
bray_matrix <- as.matrix(bray_dist)
jaccard_matrix <- as.matrix(jaccard_dist)
euclidean_matrix <- as.matrix(euclidean_dist)
```


### Principal Coordinates Analysis
La ordenación describe el proceso de producir un pequeño número de variables que pueden utilizarse para describir la relación entre un grupo de objetos, a partir de una matriz de distancias o similitudes entre los objetos, o de los valores de algunas variables medidas en cada objeto. Un ejemplo es la utilización de los primeros dos componentes principales construidos con anterioridad. En biología, se utiliza con mayor frecuencia para resumir las relaciones entre diferentes especies, determinadas a partir de su abundancia en diferentes ubicaciones, o, alternativamente, para resumir las relaciones entre diferentes ubicaciones a partir de la abundancia de diferentes especies en dichas ubicaciones. 

El análisis de coordenadas principales (Principal Coordinates Analysis) es otra manera de resumir las relaciones de interés. El análisis de coordenadas principales utiliza un enfoque de valores propios que puede considerarse una generalización del análisis de componentes principales. Parte de una matriz de distancias entre varios objetos y buscan determinar ejes de ordenación de manera que las posiciones de los objetos en una configuración t-dimensional no coinciden con las distancias o similitudes originales.

```R
# Visualizacion PCoA (Principal Coordinates Analysis)
pcoa_res <- cmdscale(beta_dist, eig = TRUE, k = 2)

# Calcular varianza explicada
var_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)

# cambiar nombre de filas
row.names(pcoa_res$points) <- c(sprintf("%02d", seq(1,48)))

# Graficar
plot(pcoa_res$points, col = feno$date2, pch = 19,
     xlab = paste0("PCoA 1 (", var_explained[1], "%)"),
     ylab = paste0("PCoA 2 (", var_explained[2], "%)"),
     main = "")
text(pcoa_res$points[,1], pcoa_res$points[,2], 
     labels = rownames(pcoa_res$points), pos = 2, cex = 0.8)

# Elipses por grupo de fecha
ordiellipse(pcoa_res$points, feno$date2, kind = "sd", draw = "polygon",
            col = c("grey","red"), alpha = 50, label = FALSE,
            border = c("grey","red"))

# Graficar a 16 X 10

# Analisis estadistico PERMANOVA (Permutational Multivariate Analysis of Variance)
adonis2(bray_dist ~ date2, data = feno, permutations = 999)
```


#### Referencias

1. Bolstad, B. M., et al. (2003). A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics, 19(2), 185–193.
2. Johnson, W. E., et al. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics, 8(1), 118–127.
3. Leek, J. T., et al. (2010). Tackling the widespread and critical impact of batch effects in high-throughput data. Nature Reviews Genetics, 11(10), 733–739.
4. Leek, J. T., & Storey, J. D. (2012). A general framework for multiple testing dependence. PNAS, 109(33), 12409–12414.
5. Love, M. I., et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.
6. Mortazavi, A., et al. (2008). Mapping and quantifying mammalian transcriptomes by RNA-Seq. Nature Methods, 5(7), 621–628.
7. Ritchie, M. E., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47.
8. Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology, 11(3), R25.
9. Zhao, S., et al. (2021). Evaluation of two main RNA-seq approaches for gene quantification in clinical RNA sequencing. NPJ Genomic Medicine, 6(1), 1–10.
10. Zhao, S., Ye, Z. and Stanton, R. (2020) ‘Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols’. Rna, 26(8), pp.903-909. Available at: 10.1261/rna.074922.120.
11. Quinn, T. P., Erb, I., Richardson, M. F., & Crowley, T. M. (2018). Understanding sequencing data as compositions: an outlook and review. Bioinformatics, 34(16), 2870–2878. Available at: https://doi.org/10.1093/bioinformatics/bty175 
12. Manly BF, Alberto JA. Multivariate statistical methods: a primer. Chapman and Hall/CRC; 2016 Nov 3.
13. Jolliffe, I. T. (2002). Principal Component Analysis (2nd ed.). New York: Springer-Verlag.
14. Pinto, M. G. Algebra Lineal con aplicaciones 7ma Stanley I Grossman.
15. James R. Schott. Matrix analysis for statistics. Third edition. | Hoboken, New Jersey : John Wiley & Sons, 2016. 
