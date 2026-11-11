# Fase 2 — Composición y diversidad alfa

# Cargamos las librerias del setup
source("scripts/00_setup.R")

# Cargamos el objeto ps de la fase 1
ps <- readRDS("data/enterotype_filt.rds")

# Nos tenemos que quedar con los géneros mas abundantes.
# Como vimos en la fase 1, solo tenemos el rango género. Igualmente vemos de nuevo los rangos taxonomicos
meta <- data.frame(sample_data(ps))
ntaxa(ps) # 131 taxones
rank_names(ps) # Solo tenemos el rango taxonomico: Genus o Genero

ps_gen <- prune_taxa(taxa_names(ps) != "-1", ps)   # solo géneros clasificados
# Nos quuedamos solo con los generos clasificados, es decir, quitamos el taxon clasfiicado como -1, pero mantenemos el ps para que el test de diversidad alfa lo tenga en cuenta

# Vamos a construir una tabla taxonomica
taxonomia <- as(tax_table(ps_gen), "matrix")
head(taxonomia)
identical(unname(taxonomia[, "Genus"]), rownames(taxonomia))

# Como el nombre del taxon es identico al genero, no necesitamos ninguna tabla de corresponencia taxon <-> genero

# Seleccionamos los diez generos clasificados con mayor abundancia total --> taxa_sums + prune_taxa
# taxa_sums() devuelve el numero total de individuos observados de cada taxa. Es equivalente a rowSums o colSums
# prune_taxa() elimina los taxa no deseados, o seleccionados

# Obtenems la abundancia de los taxones
sum_taxones <- taxa_sums(ps_gen)
head(sum_taxones)

# Top 10 taxones/genero por abundancia. En este caso, como es identico, no hace falta hacer la correspondencia
top10 <- sort(sum_taxones, decreasing = TRUE)[1:10]
top10

# Recordemos que ese top 10 es medido por abundancia relativa acumulada en todas las muestras.

# Graficos de barras
# Usamos la funcion psmelt para graficar
ps_melt <- psmelt(ps)
head(ps_melt)
names(ps_melt)

# Creamos una columna en ps_metl para indicar si el taxon esta clasificado o no
ps_melt$Clasificacion <- ifelse(
  ps_melt$OTU %in% names(top10),
  "top10",
  "Otro"
)
head(ps_melt)

# Calculamos la media de Abundance por grupo y categoria
paso1 <- aggregate(Abundance ~ Sample + Enterotype + Clasificacion, data = ps_melt, FUN = sum)
res <- aggregate(Abundance ~ Enterotype + Clasificacion, data = paso1, FUN = mean)

# Tras hacer esto nos damos cuenta que no suman 1 los enterotipos, es decir, falta en torno al 20%, que seguramente proceda del taxon -1 que quitamos al principio

# Vamos a hacer unas comprobaciones del taxon -1
tax_table(ps)["-1", ]
sin_clasificar <- as(otu_table(ps), "matrix")["-1", ]
tapply(sin_clasificar, meta[names(sin_clasificar), "SeqTech"], summary)

# En pyro el taxon sin clasificar es 0, mientras que en Illumina y Sanger ocupa de media la mitad de la muestra
# Cada tecnologia parece que se proceso con un metodo de clasificacion diferente
# Con este resultado es mejor quitar el -1 y renormalizar, para que todas las muestras sean propociones sobre lo clasificado

# Vamos a mostrar tanto la parte clasificada como la no clasificada
# Volvemos a cargar el objeto ps, de cero
ps0 <- readRDS("data/enterotype_filt.rds")

ps_cl <- prune_taxa(taxa_names(ps0) != "-1", ps0) # Nos quitamos el taxon -1
min(sample_sums(ps_cl))  
ps_cl <- transform_sample_counts(ps_cl, function(x) x / sum(x))
summary(sample_sums(ps_cl)) 

saveRDS(ps_cl, "data/ps_cl.rds")

# Calculamos el top 10 con el ps clasificado. De esta manera el ranking no esta sesgado hacia las muestras de pyro
top10 <- names(sort(taxa_sums(ps_cl), decreasing = TRUE))[1:10]

# Repetimos la parte de psmelt y aggregate, pero con ps_cl
ps_cl_melt <- psmelt(ps_cl)
head(ps_cl_melt)
names(ps_cl_melt)

# Creamos una columna en ps_metl para indicar si el taxon esta clasificado o no
ps_cl_melt$Categoria <- ifelse(ps_cl_melt$OTU %in% top10,
                               as.character(ps_cl_melt$OTU),
                                "Otros")
head(ps_cl_melt)

# Calculamos la media de Abundance por grupo y categoria
paso1 <- aggregate(Abundance ~ Sample + Enterotype + SeqTech + Categoria,
                   data = ps_cl_melt, FUN = sum)
res   <- aggregate(Abundance ~ Enterotype + SeqTech + Categoria,
                   data = paso1, FUN = mean)

aggregate(Abundance ~ Enterotype + SeqTech, data = res, FUN = sum)   # deben dar 1

# Graficamos
# Fijamos el orden de las categorias
res$Categoria <- factor(res$Categoria, levels = c(top10, "Otros"))

grafico_composicion <- ggplot(
  res,
  aes(x = Enterotype, y = Abundance, fill = Categoria)
) +
  geom_col() +
  facet_wrap(~ SeqTech) +
  labs(x = "Enterotipo", y = "Abundancia relativa media (sobre lo clasificado)")

# ACLARACIONES IMPORTANTES PARA ESTA FASE
# 1. No usar el objeto clasificado en top 10 y otros. Hay que usar el objeto ps original
# 2. ps ya esta filtrado en 131 taxones. observed describira cuantos taxones quedaron tras ese filtrado y estan presentes en cada muestra, no la riqueza original antes del filtrado
# 3. Las abundacias son proporciones relativas. Esto permite caculas Shannon y Simpson, pero no se interpreta Chao1 ni ACE como estimaciones fiables de riqueza.
# 4. Hay que considerar el efecto de SeqTech al comparar grupos. Hay que tener en cuenta esta varibale en el diseño del analsis, o bien comparar dentro de cada tecnologia
# 5. Hay variables con las uqe no merece la pena trabjar el estudio, como las demograficas
# 6. Considerar el taxon -1

# Diversidad alfa
# Shannon y Simpson describen la diversidad de los taxones retenidos:
# Shannon combina riqueza y equitatividad; Simpson da mas peso a los taxones abundantes.
# estimate_richness puede advertir que no hay singletons; los resultados se interpretan
# unicamente para los generos retenidos tras el filtrado de la fase 1.

# Se va a llevar a cabo un análisis principal la diversidad solo de géneros clasificados, renormalizada.
# Tambien se va a incluir otro analsis incluyendo -1 como sensibilidad, dejando claro que es una categoría agregada de identidad incierta
# Por tanto e trabajara con ps-cl y ps0

# Shannon y Simpson se pueden calcular con las proporciones
div_cl <- estimate_richness(
  ps_cl,
  measures = c("Shannon", "Simpson")
)

div_sensibilidad <- estimate_richness(
  ps0,
  measures = c("Shannon", "Simpson")
)

# Observed: número de taxones presentes por muestra
observed_por_muestra <- function(ps_obj) {
  otu <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) otu <- t(otu)
  rowSums(otu > 0)
}

div_cl$Observed <- observed_por_muestra(ps_cl)[rownames(div_cl)]
div_sensibilidad$Observed <-
  observed_por_muestra(ps0)[rownames(div_sensibilidad)]

# Vamos a añadir variables a los df de riquezas
meta_cl <- data.frame(sample_data(ps_cl))

div_cl$Muestra <- rownames(div_cl)
div_cl$Enterotype <- meta_cl[
  div_cl$Muestra,
  "Enterotype"
]
div_cl$SeqTech <- meta_cl[
  div_cl$Muestra,
  "SeqTech"
]

meta_sens <- data.frame(sample_data(ps0))

div_sensibilidad$Muestra <- rownames(div_sensibilidad)
div_sensibilidad$Enterotype <- meta_sens[
  div_sensibilidad$Muestra,
  "Enterotype"
]
div_sensibilidad$SeqTech <- meta_sens[
  div_sensibilidad$Muestra,
  "SeqTech"
]

# Boxplots por enterotipo, separados por tecnología de secuenciación
box_shannon_cl <- ggplot(
  div_cl,
  aes(x = Enterotype, y = Shannon, fill = SeqTech)
) +
  geom_boxplot() +
  facet_wrap(~ SeqTech) +
  labs(
    title = "Shannon: géneros clasificados",
    x = "Enterotipo",
    y = "Índice de Shannon",
    fill = "Tecnología"
  ) +
  scale_fill_brewer(palette = "Set2")

box_simpson_cl <- ggplot(
  div_cl,
  aes(x = Enterotype, y = Simpson, fill = SeqTech)
) +
  geom_boxplot() +
  facet_wrap(~ SeqTech) +
  labs(
    title = "Simpson: géneros clasificados",
    x = "Enterotipo",
    y = "Índice de Simpson",
    fill = "Tecnología"
  ) +
  scale_fill_brewer(palette = "Set2")

box_shannon_sens <- ggplot(
  div_sensibilidad,
  aes(x = Enterotype, y = Shannon, fill = SeqTech)
) +
  geom_boxplot() +
  facet_wrap(~ SeqTech) +
  labs(
    title = "Shannon: incluyendo el taxón -1",
    x = "Enterotipo",
    y = "Índice de Shannon",
    fill = "Tecnología"
  ) +
  scale_fill_brewer(palette = "Set2")

box_simpson_sens <- ggplot(
  div_sensibilidad,
  aes(x = Enterotype, y = Simpson, fill = SeqTech)
) +
  geom_boxplot() +
  facet_wrap(~ SeqTech) +
  labs(
    title = "Simpson: incluyendo el taxón -1",
    x = "Enterotipo",
    y = "Índice de Simpson",
    fill = "Tecnología"
  ) +
  scale_fill_brewer(palette = "Set2")

box_shannon_cl
box_simpson_cl
box_shannon_sens
box_simpson_sens

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/phase_2_composicion_generos.png",
       grafico_composicion, width = 10, height = 6, dpi = 300)
ggsave("results/figures/phase_2_shannon_generos_clasificados.png",
       box_shannon_cl, width = 10, height = 6, dpi = 300)
ggsave("results/figures/phase_2_simpson_generos_clasificados.png",
       box_simpson_cl, width = 10, height = 6, dpi = 300)
ggsave("results/figures/phase_2_shannon_incluye_taxon_menos_1.png",
       box_shannon_sens, width = 10, height = 6, dpi = 300)
ggsave("results/figures/phase_2_simpson_incluye_taxon_menos_1.png",
       box_simpson_sens, width = 10, height = 6, dpi = 300)
