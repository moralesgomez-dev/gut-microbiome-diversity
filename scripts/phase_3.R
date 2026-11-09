# Fase 3: Diversidad beta

# Cargamos las librerias del setup
source("scripts/00_setup.R")

# Cargamos el objeto ps de la fase 2 (sin el taxon -1)
ps_cl <- readRDS("data/ps_cl.rds")

# Calculamos la distancia Bray-Curtis --> distance(), que toma un objeto ps y un metodo, en este caso bray
distance_bc <- distance(ps_cl, method = "bray")

# Tambien calculamos la ordenacion por PCoA, que representa en un mapa 2D respetando las distancias de bc
# Lo hacemos con ordinate(), que toma un objeto ps, un metodo, que será PCoA y las distancias calculadas
ordination <- ordinate(ps_cl, method = "PCoA", distance = distance_bc)

# Vamos a representarlo con plot_ordination(), que toma un ps, la ordenacion y mas configs de grafico
grafico_pcoa_seqtech <- plot_ordination(
  ps_cl, ordination = ordination, color = "SeqTech"
)
grafico_pcoa_enterotype <- plot_ordination(
  ps_cl, ordination = ordination, color = "Enterotype"
)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/phase_3_pcoa_seqtech.png",
       grafico_pcoa_seqtech, width = 8, height = 6, dpi = 300)
ggsave("results/figures/phase_3_pcoa_enterotype.png",
       grafico_pcoa_enterotype, width = 8, height = 6, dpi = 300)

# Vemos como segun el tipo que sea la muestra, por seqtech o enterotype, se agrupan de manera diferente y separada, generalmente

# Vamos a constrastar la separacion entre grupos con PERMANOVA
# Basicamente vemos si los centroides de los grupos difeiren mas de lo esperable al azar.
# Adonis no admite NA, por lo que tendremos que quitar los NA que hay presentes en Enterotype y recalcular distancias

ps_t <- subset_samples(ps_cl, !is.na(Enterotype))
dist_t <- distance(ps_t, method = "bray")
meta_t <- data.frame(sample_data(ps_t))
all(labels(dist_t) == rownames(meta_t))

# Adonis toma como formula una matriz de distancias (ya calculada) asi como variables que incluir en el analisis (~), data el df de origen de los datos,
set.seed(123)
adonis2(dist_t ~ Enterotype + SeqTech, data = meta_t,
        permutations = 999, by = "margin")

# Tras estos resultados podriamos concluir que  tanto el enterotipo como la tecnología se asocian con diferencias en la composición según Bray–Curtis; 
# la tecnología tiene un efecto marginal algo mayor en este modelo.

# Sin embargo, PERMANOVA puede salir significativo porque los grupos estan en distintos sitios (segun composicion) o porque unos grupos estan mas dispersos que otros.
# Por ello vamos a comprobar con betadisper, que mide la distancia de cada uestra al centro de su grupo, junto con permutest, que constrasta si esas dispersiones difieren

# betadisper toma las distancias, grupo que generalmente es el factor (enterotype en este caso). Tambien puede tomar otros parametros como tipo, ajuste de bias...
bd <- betadisper(dist_t, meta_t$Enterotype)

# Toma como objeto principal un de la clase betadisper (resultado de hacer betadisper) y el numero de permutaciones
permutest(bd, permutations = 999)
png("results/figures/phase_3_dispersion_enterotype.png",
    width = 2400, height = 1800, res = 300)
boxplot(bd)
dev.off()

# Comprobamos tambien si la dispersion difiere entre tecnologias
bd_seqtech <- betadisper(dist_t, meta_t$SeqTech)
permutest(bd_seqtech, permutations = 999)
png("results/figures/phase_3_dispersion_seqtech.png",
    width = 2400, height = 1800, res = 300)
boxplot(bd_seqtech)
dev.off()

# Una prueba significativa indica dispersion desigual entre grupos; no identifica
# por si sola que grupos difieren ni demuestra diferencias entre sus centroides.

# Interpretar el PERMANOVA con cautela si la dispersion es desigual.
