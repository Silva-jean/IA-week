library(tidyverse)
library(AGHmatrix)
library(asreml)

# 1. Ler os fenótipos e os marcadores
dados = "C:/Users/jeanP/Documents/MEGA/IA-week/data"
arquivo_mh = file.path(dados, "marcadores_hibridos.rds")
arquivo_matrizes = file.path(dados, "matrizes_genomicas.rds")

phen = read.csv(file.path(dados, "blue_data.csv")) %>%
  rename(env = Field_Location, hib = Pedigree)

mh = readRDS(arquivo_mh)
ml = readRDS(file.path(dados, "marcadores_linhagens.rds"))

# 2. Identificar os híbridos sem correspondência em mh
ids_fenotipos = unique(phen$hib)
ids_marcadores = rownames(mh)

stopifnot(!is.null(ids_marcadores), !anyDuplicated(ids_marcadores))

print(table(tem_marcadores = ids_fenotipos %in% ids_marcadores))
ids_sem_correspondencia = setdiff(ids_fenotipos, ids_marcadores)

# 3. Conferir os parentais dos híbridos sem correspondência
parentais = phen %>%
  filter(hib %in% ids_sem_correspondencia) %>%
  distinct(hib, DH1, DH2) %>%
  mutate(
    DH1_em_ml = DH1 %in% rownames(ml),
    DH2_em_ml = DH2 %in% rownames(ml)
  )

print(as.data.frame(parentais), row.names = FALSE)

# 4. Manter as observações dos híbridos com marcadores
phen = phen %>%
  filter(!is.na(hib), hib %in% ids_marcadores)

# 5. Alinhar os marcadores à ordem dos híbridos em phen
ids_hibridos = unique(phen$hib)
mh = mh[ids_hibridos, , drop = FALSE]

# 6. Conferir os IDs e a ordem
stopifnot(
  length(ids_hibridos) > 0L,
  identical(rownames(mh), ids_hibridos)
)
print(length(ids_hibridos))

# 7. Carregar as matrizes salvas ou recalcular quando os insumos mudarem
md5_mh = unname(tools::md5sum(arquivo_mh))
stopifnot(!is.na(md5_mh))

matrizes = NULL
if (file.exists(arquivo_matrizes)) {
  matrizes = tryCatch(
    readRDS(arquivo_matrizes),
    error = function(e) NULL
  )
}

matrizes_validas =
  is.list(matrizes) &&
  identical(matrizes$md5_mh, md5_mh) &&
  identical(rownames(matrizes$aditiva), ids_hibridos) &&
  identical(colnames(matrizes$aditiva), ids_hibridos) &&
  identical(rownames(matrizes$dominancia), ids_hibridos) &&
  identical(colnames(matrizes$dominancia), ids_hibridos)

if (matrizes_validas) {
  message("Matrizes genômicas carregadas de ", arquivo_matrizes)
} else {
  matrizes = list(
    md5_mh = md5_mh,
    aditiva = Gmatrix(mh, method = "VanRaden"),
    dominancia = Gmatrix(mh, method = "Vitezica")
  )
  saveRDS(matrizes, arquivo_matrizes)
  message("Matrizes genômicas calculadas e salvas em ", arquivo_matrizes)
}

matriz_aditiva = matrizes$aditiva
matriz_dominancia = matrizes$dominancia

print(dim(matriz_aditiva))
print(dim(matriz_dominancia))

# 8. Resumir os cruzamentos e a frequência dos parentais
cruzamentos = phen %>%
  filter(!is.na(DH1), !is.na(DH2)) %>%
  distinct(DH1, DH2)

freq_DH1 = cruzamentos %>%
  count(DH1, name = "n_cruzamentos", sort = TRUE)

freq_DH2 = cruzamentos %>%
  count(DH2, name = "n_cruzamentos", sort = TRUE)

freq_genitores = cruzamentos %>%
  mutate(id_cruzamento = row_number()) %>%
  pivot_longer(c(DH1, DH2), values_to = "genitor") %>%
  distinct(id_cruzamento, genitor) %>%
  count(genitor, name = "n_cruzamentos", sort = TRUE) %>%
  mutate(percentual_cruzamentos = 100 * n_cruzamentos / nrow(cruzamentos))

resumo_cruzamentos = tibble(
  n_cruzamentos = nrow(cruzamentos),
  n_genitores_DH1 = nrow(freq_DH1),
  n_genitores_DH2 = nrow(freq_DH2),
  n_genitores_total = nrow(freq_genitores)
)

print(resumo_cruzamentos)
print(head(freq_DH1))
print(head(freq_DH2))
print(head(freq_genitores))

# 9. Classificar os parentais e identificar os híbridos SS × NSS
# ==============================================================================
# Categoria        | Critério
# -----------------|------------------------------------------------------------
# Stiff_Stalk      | Classificação individual adotada para as referências SS
# Non_Stiff_Stalk  | Referências NSS e DHs da população GEM Derived NSS do G2F
# Tropical_Exotico | Ancestralidade GEM atribuída aos IDs GEMN-XXXX e BGEM-XXXX-N
# Indefinido       | Evidência insuficiente para atribuição
#
# Tropical_Exotico indica ancestralidade. Parte do germoplasma GEM, incluindo
# derivados BGEM, foi desenvolvida por introgressão e retrocruzamentos com
# linhagens elite temperadas.
#
# Referências
# G2F, Project Overview & Scope: https://www.genomes2fields.org/about/project-overview/
# McFarland et al. (2020): https://doi.org/10.1186/s13104-020-4922-8
# Mikel & Dudley (2006): https://doi.org/10.2135/cropsci2005.10-0371
# Tracy & Chandler (2006): https://www.researchgate.net/publication/230027777
# White et al. (2020): https://doi.org/10.1002/csc2.20050
# USDA–ARS GRIN, LH244: https://npgsweb.ars-grin.gov/gringlobal/accessiondetail?id=1594201
# USDA–ARS GRIN, exemplo de retrocruzamento GEM:
# https://npgsweb.ars-grin.gov/gringlobal/accessiondetail?id=1982870
# DK3IIH6, Tabela 1: https://pmc.ncbi.nlm.nih.gov/articles/PMC13153409/
# ==============================================================================

df_pedigree = tibble(
  Linhagem = unique(trimws(c(as.character(phen$DH1), as.character(phen$DH2))))
) %>%
  filter(!is.na(Linhagem), nzchar(Linhagem)) %>%
  mutate(
    origem_gem_nss = grepl(
      "^GEMN-(0096|0097|0192|0225)_(LH212HT|PHJ89|PHK76|PHN46|PHP02)_[0-9]+$",
      Linhagem,
      ignore.case = TRUE
    ),
    Grupo_Heterotico = case_when(
      Linhagem %in% c(
        "B14A", "B37", "B73", "LH195", "LH244",
        "PHB47", "PHJ40", "PHRE1", "PHW52"
      ) ~ "Stiff_Stalk",
      Linhagem %in% c(
        "B97", "DK3IIH6", "OH43", "LH212HT", "PHN46",
        "PHJ89", "PHK76", "PHP02", "LH82", "LH185",
        "PHG29", "PHG47", "PHK56", "PHM49", "PHN82",
        "PHR03", "PHW30", "PHZ51"
      ) | origem_gem_nss ~ "Non_Stiff_Stalk",
      grepl(
        "^(GEMN-[0-9]{4}|BGEM-[0-9]{4}-N)$",
        Linhagem,
        ignore.case = TRUE
      ) ~ "Tropical_Exotico",
      TRUE ~ "Indefinido"
    ),
    Evidencia = case_when(
      origem_gem_nss ~ "Origem na população GEM Derived NSS",
      Grupo_Heterotico == "Tropical_Exotico" ~ "Origem GEM atribuída pelo ID",
      Grupo_Heterotico == "Indefinido" ~ "Evidência insuficiente",
      TRUE ~ "Classificação individual adotada"
    ),
    Grupo_Heterotico = factor(
      Grupo_Heterotico,
      levels = c("Stiff_Stalk", "Non_Stiff_Stalk", "Tropical_Exotico", "Indefinido")
    )
  ) %>%
  select(Linhagem, Grupo_Heterotico, Evidencia) %>%
  arrange(Grupo_Heterotico, Linhagem)

resumo_grupos = df_pedigree %>%
  count(Grupo_Heterotico, name = "n_linhagens", .drop = FALSE)

print(resumo_grupos)
if (interactive()) View(df_pedigree)

grupos = setNames(
  as.character(df_pedigree$Grupo_Heterotico),
  df_pedigree$Linhagem
)

phen = phen %>%
  mutate(
    grupo_DH1 = unname(grupos[trimws(as.character(DH1))]),
    grupo_DH2 = unname(grupos[trimws(as.character(DH2))])
  )

hibridos_ss_nss = phen %>%
  filter(
    (grupo_DH1 == "Stiff_Stalk" & grupo_DH2 == "Non_Stiff_Stalk") |
      (grupo_DH1 == "Non_Stiff_Stalk" & grupo_DH2 == "Stiff_Stalk")
  ) %>%
  distinct(hib, DH1, DH2, grupo_DH1, grupo_DH2)

print(n_distinct(hibridos_ss_nss$hib))

# 10. Preparar os dados para o ASReml
dados_modelo = phen %>%
  mutate(
    env = factor(trimws(as.character(env))),
    hib = factor(hib, levels = rownames(matriz_aditiva))
  ) %>%
  arrange(env) %>%
  as.data.frame()

stopifnot(
  !anyNA(dados_modelo$env),
  all(nzchar(as.character(dados_modelo$env))),
  !anyNA(dados_modelo$hib),
  is.numeric(dados_modelo$BLUE)
)

print(nlevels(dados_modelo$env))
print(table(dados_modelo$env[!is.na(dados_modelo$BLUE)]))


# 7.1. Regularizar, inverter e salvar a matriz aditiva
A = matriz_aditiva
ajuste_diagonal = 1e-6 * mean(diag(A))

A_regularizada = A
diag(A_regularizada) = diag(A_regularizada) + ajuste_diagonal

matriz_aditiva_inversa = chol2inv(chol(A_regularizada))
dimnames(matriz_aditiva_inversa) = dimnames(A)
attr(matriz_aditiva_inversa, "INVERSE") = TRUE

arquivo_inversa = file.path(dados, "matriz_aditiva_inversa.rds")
saveRDS(matriz_aditiva_inversa, arquivo_inversa)

print(ajuste_diagonal)
print(arquivo_inversa)