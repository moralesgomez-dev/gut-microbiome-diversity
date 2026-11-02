data(enterotype)
enterotype

# Contiene tres objetos: otu_table, sample_data y tax_table. Vamos a ir viendo cada una
OTU1 <- otu_table(enterotype)
sample1 <- sample_data(enterotype)
TAX1 <- tax_table(enterotype)

# Viendo OTU1 nos damos cuenta de que la tabla de abundancia viene en proporciones relativas

# nsamples() devuelve el numero de muestras
nsamples(OTU1)
nsamples(sample1)
nsamples(TAX1) # No tiene sentido, ya que tax_table no tiene muestras, sino taxones. Por lo tanto, no tiene nsamples.
nsamples(enterotype)

# ntaxa() devuelve el numero de taxones
ntaxa(OTU1)
ntaxa(sample1) # No tiene sentido, ya que sample_data no tiene taxones, sino muestras. Por lo tanto, no tiene ntaxa.
ntaxa(TAX1)
ntaxa(enterotype)

# sample_names() devuelve los nombres de las muestras
sample_names(OTU1)
sample_names(sample1)
sample_names(enterotype)

# sample_variables() devuelve los nombres de las variables de metadatos
sample_variables(OTU1) # No tiene sentido, ya que otu_table no tiene variables de metadatos, sino OTUs. Por lo tanto, no tiene sample_variables.
sample_variables(sample1)
sample_variables(TAX1) # No tiene sentido, ya que tax_table no tiene variables de metadatos, sino taxones. Por lo tanto, no tiene sample_variables.
sample_variables(enterotype)

# Comprobamos cuantos NA hay en los metadatos
sum(is.na(sample1))

meta_inicial <- data.frame(sample1)

# Da la casuliadad que las muestras que hemos mirado tiene NA en algunas variables. Vqamos a ver cuantos NA hay en esas variables
colSums(is.na(meta_inicial[, c("Enterotype", "SeqTech", "Nationality", "Gender", "Age", "ClinicalStatus")]))

# Solo SeqTech y Etnterotype estan completas. Las otras tienen mucho porcentaje de NA

# Vamos a ver si hay una coincidencia: 241 valores faltantes es un numero parecido al numero de muestras que no son sanger
table(meta_inicial$SeqTech, !is.na(meta_inicial$Nationality))

# Las 39 muestras con Nationality son todas Sanger

table(meta_inicial$SeqTech, !is.na(meta_inicial$Gender))
table(meta_inicial$SeqTech, !is.na(meta_inicial$ClinicalStatus))

# Ocurre igual con Gender y CLinical Status
# TODO README: Las variables demograficas y clinicas solo estan disponibles para 39 muestras Sangerm por lo que no se usan para el analisis principal

# Por ultimo vemos como se reparte las tecnologias y los enterotipos
table(meta_inicial$SeqTech, meta_inicial$Enterotype, useNA = "ifany")

# Los tres enterotipos apartecen en las tres tecnologias, pero el reparto no es igual en las tres.
# En sanger domina el 3 y el 1 es minoritario
# Por tanto, al comparar comunidades por enterotipo, habra que incluir seqtech en el modelo para n atribuir al enterotipo lo que e parte es de la tecnica

# Filtramos taxones raros y muestras pobres
otu <- as(otu_table(enterotype), "matrix")

# filter_taxa() filtra los taxones segun criterios de abundacia de OTU entre muestras. 
# Requiere un objeto phyloseq (enterotype en este caso) y flist (una funcion que devuelva un vector logico de longitud igual al numero de taxones).
# Basicamente esa funcion va a determionar si un taxon se queda o no en el objeto phyloseq.
# Por otro lado tenemos prune_samples() que filtra las muestras segun criterios de abundancia de OTU entre taxones. Podemos pasarle un vetor logico,
# procedente de una funcion que determine si una muestra se queda o no en el objeto phyloseq.
# Tambien podemos usar filter_taxa() con PRUNE=TRUE para que automaticamente elimine las muestras que no tengan ningun taxon tras el filtrado de taxones.

# Antes de eso, tenemos que elegir un criterio de filtrado.
# De antes ya sabemos que tenemos 553 taxones y 280 muestras

# Vamos a ver como se reparte la prevalencia, es decir, en cuantas muestras aparece cada taxon
# taxa_are_rows() es un slot del objeto que indica si los taxones estan en filas o columnas. 
taxa_are_rows(enterotype) # TRUE, por lo que los taxones estan en filas y las muestras en columnas.
# Vemos los nombres de los niveles taxonomicos, con rank_names() que devuelve un vector con los nombres de los niveles taxonomicos, en orden jerarquico.
rank_names(enterotype) 
# Vemos un resumen de la abundancia de las muestras
summary(colSums(otu)) # media y mediana de 0.999
summary(sample_sums(enterotype))

min(otu[otu > 0]) # Minimo de 0.00000001
quantile(otu[otu > 0], c(0.01, 0.05, 0.25, 0.5))
# 1%         5%         25%        50% 
# 6.4000e-07 2.0600e-06 7.9500e-06 2.0715e-05 

# Vamos a hacer una tabla de sensibilidad, donde cada celda es el numero de taxones que se quedarian.
# Comparamos diferentes valores de abundancia minima y fraccion de muestras.
detection <- c(0, 1e-5, 1e-4, 1e-3, 1e-2) # abundancia relativa mínima
prev_min  <- c(0.01, 0.05, 0.10, 0.20) # fracción de muestras
tab <- sapply(detection, function(d) {
  prev_d <- rowMeans(otu > d)               # fracción de muestras donde supera d
  sapply(prev_min, function(p) sum(prev_d >= p))
})
dimnames(tab) <- list(prev_min = prev_min, detection = detection)
tab

# Abundancia conservada por muestra tras el filtrado de taxones
combos <- expand.grid(d = c(1e-5, 1e-4, 1e-3), p = c(0.05, 0.10))
res <- t(mapply(function(d, p) {
  keep <- rowMeans(otu > d) >= p
  ab <- colSums(otu[keep, , drop = FALSE])   
  c(taxa = sum(keep), mediana = median(ab), minimo = min(ab))
}, combos$d, combos$p))
cbind(combos, res)

# # Filtro de taxones: abundancia relativa > 1e-4 en al menos el 5 % de las muestras (14 de 280).
# 553 -> 131 taxones. Por qué:
# - x > 0 no discrimina: 500 de 553 taxones pasan incluso exigiendo el 20 % de las muestras.
# - Retención mediana por muestra 0,9988 (mínimo 0,82); con d = 1e-3 o p = 10 % el mínimo cae a ~0,64-0,65.
# - p = 5 % (14 muestras) queda por debajo del grupo más pequeño (enterotipo 2, 35 muestras).
# Límite: la retención mide cuánta abundancia se pierde, no demuestra que lo eliminado sea ruido.

# Filtro de muestras: no se elimina ninguna (prune_samples no se aplica).
# Retención mínima 0,82 sin separación clara del resto; el nº de taxones detectados depende de
# la tecnología (Illumina 96, Sanger ~60, Pyro454 ~37), no de la calidad de la muestra.

# Aplicamos el filtro con filter_taxa() y prune_samples() para eliminar los taxones raros y las muestras pobres.
d <- 1e-4
p <- 0.05

flist <- function(x) mean(x > d) >= p
enterotype_filt <- filter_taxa(enterotype, flist, prune = TRUE)
ntaxa(enterotype_filt) # 131

# Para ver si eliminamos alguna muestra, vamos a ver cuanto conserva cada muestra, y cuantas de los 131 taxones tiene cada muestra.
retenida <- sample_sums(enterotype_filt)
detectados <- colSums(as(otu_table(enterotype_filt), "matrix") > 0)

summary(retenida)
sort(retenida)[1:5]
summary(detectados)
sort(detectados)[1:5]

hist(retenida, breaks = 30)
hist(detectados, breaks = 30)

# La mayoria esta por encima del 0.99 u ñla mediana es de 0.9988. Alguna muestra esta por debado pero conservan aun asi mas del 80%, razonable.
# El maximo de taxones detectados es 69 y el q3 tambien es 96. Esto puede significar que mas de ua cuarta parte de las muestras tiene exactamente 96 taxones.
# Parece ser que hay dos tipos de muestras, una con muchos taxones y otras con pocos. Lo comprobamos

meta <- data.frame(sample_data(enterotype_filt))
meta$retenida   <- retenida[rownames(meta)]
meta$detectados <- detectados[rownames(meta)]

sample_variables(enterotype_filt)               # mira si hay una variable de tecnología
table(meta$SeqTech)
tapply(meta$detectados, meta$SeqTech, summary)  # detectados según la tecnología

meta[names(sort(retenida)[1:3]), ]              # metadatos de las 3 de menor retención

# Segun la tecnologia utilizada cambia bastnate el conteo de taxones
# Esto es para tenerlo en cuenta para las siguientes fses de anlisis de composiciom
# TODO README: El numero de taxones varia en funcion de la tecnologia, por lo que las comparacioes de riqueza estan confundidas por la tecnica.

# Seguimos con el flujo noroma: Renormalizar

enterotype_filt <- transform_sample_counts(enterotype_filt, function(x) x / sum(x))
summary(sample_sums(enterotype_filt))   


saveRDS(enterotype_filt, "data/enterotype_filt.rds")


