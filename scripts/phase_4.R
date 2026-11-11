# Fase 4: ML ligero a traves de clustering
# Se aplicara clustering no supervisado con cluster::pam sobre la amtriz bray.curtis
# Colorear PCoA por cluster y comparar con la etiqueta de los metadatos mediante una tablad e contigencia

# Cargamos las librerias del setup
source("scripts/00_setup.R")
library(cluster)

# Cargamos el objeto ps de la fase 2 (sin el taxon -1)
ps_cl <- readRDS("data/ps_cl.rds")

# Calculamos de nuevo las distancias y el ord para despues
distance_bc <- distance(ps_cl, method = "bray")
ordination <- ordinate(ps_cl, method = "PCoA", distance = distance_bc)

# cluster::pam es una version mas robusta de un K-means 
# Tambien nos permite seleccionar el numero de clusteres con pr$silinfo$avg.widht
# Los parametros mas importantes son x, la matriz de distancias, y k, el numero de clusteres seleccionados

# Calculamos la silueta para varios k
ks <- 2:6
sil <- sapply(ks, function(k) pam(distance_bc, k = k)$silinfo$avg.width)
data.frame(k = ks, silueta_media = sil)

# Parace que el mejor K es el 3, si bien no difiere mucho de 4. Tomamos k=3
# Hacemos ya pam para k=3
fit <- pam(distance_bc, k = 3)
meta <- data.frame(sample_data(ps_cl))
# fit$clustering nos dice a que grupo pertenece cada muestra, segun lo clusterizado
# Alineamos por nombre de muestra
cl <- fit$clustering[rownames(meta)] 

# Vamos a ver si cada grupo hecho por cluster se concentra en algun enterotipo. Tambien se va a ahcer por tecnologia
# Construimos la tabla de contingencia
tab_ent <- table(cluster = cl, Enterotype = meta$Enterotype, useNA = "ifany")
tab_ent

# Calculamos el reparto de enterotipos dentro de cada cluster
round(prop.table(tab_ent, 1), 2)

# Hacemos lo mismo para SeqTech
tab_seqtech <- table(cluster = cl, SeqTech = meta$SeqTech, useNA = "ifany")
tab_seqtech
round(prop.table(tab_seqtech, 1), 2)

# Coloreamos el PCoA por cluster
sample_data(ps_cl)$cluster <- factor(cl[sample_names(ps_cl)])
grafico_pcoa_clustering <- plot_ordination(ps_cl, ordination = ordination, color = "cluster")
ggsave("results/figures/phase_4_pcoa_clustering.png",
       grafico_pcoa_clustering, width = 8, height = 6, dpi = 300)

# Vemos que hay cierta correspondencia entre cluster y enterotupos, sobre todo para cl3 y ent2, pero no hay una separacion uno-uno clara
# Segun Seqtech, wel cluster 2 si es claramente dominante la tec Pyro, pero los cl1 y 3 estan repartidos
# Esto puede sugerir que el cl basado en bray curtis capta la tecnologia de sencuenciacion.

# Por tanto, los clusters no reproducen limpiamnte los enterotipos y parecen estar influidos por la tecnologia. Igualmente esto no prueba por si solo que las asociaciones sean significativas.
