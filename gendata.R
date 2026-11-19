
n_nucleos = 5L  # Ajuste conforme os núcleos e a memória disponíveis.

preparacao = local({
  entrada = NULL
  pasta_temp = NULL
  provis = NULL
  tryCatch({
    pasta = file.path("G2F", "GenomesToFields_G2F_genotypic_data_2014_to_2023")
    arquivo_vcf = file.path(pasta, "inbreds_G2F_2014-2023_437k.vcf")
    arquivo_blue = file.path("data", "blue_data.csv")
    arquivo_chave = file.path(pasta, "key_inbreds_G2F_2014-2023.txt")
    destino = file.path("data", "snp.clean")
    entradas = c(arquivo_vcf, arquivo_blue, arquivo_chave)
    if (!all(file.exists(entradas))) {
      stop("Arquivos ausentes: ", paste(entradas[!file.exists(entradas)], collapse = ", "))
    }
    # 1. Uma linha por híbrido basta para encontrar os parentais.
    # As médias por local continuam disponíveis no CSV original.
    colunas = c("Pedigree", "DH1", "DH2")
    medias = read.csv(arquivo_blue, colClasses = setNames(rep("character", 3L), colunas))
    hibridos = unique(medias[, colunas])
    stopifnot(nrow(hibridos) > 0L, !anyNA(hibridos$Pedigree),
              all(nzchar(hibridos$Pedigree)), !anyDuplicated(hibridos$Pedigree))
    hibridos$DH1 = trimws(hibridos$DH1)
    hibridos$DH2 = trimws(hibridos$DH2)

    # 2. Leia o cabeçalho para localizar as amostras sem carregar os marcadores.
    entrada = file(arquivo_vcf, "rt")
    metadados = character()
    repeat {
      linha = readLines(entrada, n = 1L, warn = FALSE)
      if (!length(linha)) stop("Cabeçalho #CHROM ausente.")
      if (startsWith(linha, "#CHROM\t")) break
      metadados = c(metadados, linha)
    }
    cabecalho = strsplit(linha, "\t", fixed = TRUE)[[1L]]
    amostras = cabecalho[-seq_len(9L)]
    stopifnot(length(amostras) > 0L, !anyDuplicated(amostras))

    # 3. Conte os pares exatos e trate somente a exceção revisada.
    # A chave oficial sugere nomes alternativos para revisão manual.
    exatos = (hibridos$DH1 %in% amostras) + (hibridos$DH2 %in% amostras)
    correcao = hibridos$Pedigree == "PHG80xPHP02" &
      (is.na(hibridos$DH1) | hibridos$DH1 == "") &
      (is.na(hibridos$DH2) | hibridos$DH2 == "")
    if (any(correcao) && all(c("PHG80", "PHP02") %in% amostras)) {
      hibridos$DH1[correcao] = "PHG80"
      hibridos$DH2[correcao] = "PHP02"
    }
    hibridos$Parentais_no_VCF = (hibridos$DH1 %in% amostras) +
      (hibridos$DH2 %in% amostras)
    pendentes = hibridos[hibridos$Parentais_no_VCF < 2L, , drop = FALSE]
    elegiveis = hibridos[hibridos$Parentais_no_VCF == 2L, , drop = FALSE]
    ids = sort(unique(c(elegiveis$DH1, elegiveis$DH2)))
    if (!length(ids)) stop("Nenhum híbrido tem os dois parentais no VCF.")

    chave = read.delim(arquivo_chave, check.names = FALSE, colClasses = "character",
                       na.strings = c("", "NA"))
    nomes = trimws(chave[["Cultivar"]])
    alternativos = trimws(chave[["Alternative name"]])
    stopifnot(length(nomes) == nrow(chave), length(alternativos) == nrow(chave))
    faltantes = unique(c(pendentes$DH1, pendentes$DH2))
    faltantes = sort(setdiff(faltantes[!is.na(faltantes) & nzchar(faltantes)], amostras))
    candidatos_vcf = character(length(faltantes))
    for (i in seq_along(faltantes)) {
      id = faltantes[i]
      candidatos = c(nomes[which(alternativos == id)], alternativos[which(nomes == id)])
      candidatos_vcf[i] = paste(intersect(candidatos, amostras), collapse = " | ")
    }
    revisao = data.frame(ID = faltantes, Candidatos_no_VCF = candidatos_vcf)
    print(data.frame(
      Parentais_no_VCF = 0:2,
      Hibridos_match_exato = tabulate(exatos + 1L, nbins = 3L),
      Hibridos_apos_revisao = tabulate(hibridos$Parentais_no_VCF + 1L, nbins = 3L)
    ), row.names = FALSE)
    print(pendentes, row.names = FALSE)
    print(revisao, row.names = FALSE)
    message("Amostras parentais distintas: ", length(ids))

    # 4. Faça o QC em lotes. Conte a frequência no painel de parentais.
    # MAF >= 0,05 elimina marcadores raros. Exigir GT completo evita imputação
    # adicional nesta análise. Os genótipos da publicação já passaram pelo PHG.

    pasta_temp = tempfile("g2f_qc_")
    stopifnot(dir.create(pasta_temp))
    col_amostras = 9L + match(ids, amostras)
    num_amostras = length(ids)
    cromossomos = as.character(1:10)
    arquivos = file.path(pasta_temp, paste0(cromossomos, ".tsv"))
    binarios = file.path(pasta_temp, paste0(cromossomos, ".bin"))
    codigos = c("0/0", "0/1", "1/0", "1/1")
    dosagens = c(0L, 1L, 1L, 2L)
    bases = c("A", "C", "G", "T")
    contagens = c(Variantes = 0L, SNPs_bialelicos = 0L, QC = 0L)
    contagem_chr = integer(length(cromossomos))
    repeat {
      linhas = readLines(entrada, n = 500L, warn = FALSE)
      if (!length(linhas)) break
      contagens["Variantes"] = contagens["Variantes"] + length(linhas)
      registros = strsplit(linhas, "\t", fixed = TRUE)
      saida_lote = character(length(registros))
      chr_lote = integer(length(registros))
      dose_lote = vector("list", length(registros))

      for (i in seq_along(registros)) {
        registro = registros[[i]]
        if (length(registro) != length(cabecalho)) stop("Registro VCF incompleto.")
        if (!(registro[4L] %in% bases && registro[5L] %in% bases &&
              registro[4L] != registro[5L])) next
        contagens["SNPs_bialelicos"] = contagens["SNPs_bialelicos"] + 1L
        formatos = strsplit(registro[9L], ":", fixed = TRUE)[[1L]]
        pos_gt = match("GT", formatos)
        if (is.na(pos_gt)) next

        chamadas = registro[col_amostras]
        gt = if (pos_gt == 1L) {
          sub(":.*$", "", chamadas)
        } else {
          vapply(strsplit(chamadas, ":", fixed = TRUE), function(x) {
            x[pos_gt]
          }, character(1L))
        }
        dose = dosagens[match(chartr("|", "/", gt), codigos)]
        if (anyNA(dose)) next
        freq_alt = sum(dose) / (2 * num_amostras)
        if (min(freq_alt, 1 - freq_alt) < 0.05) next

        chr = match(registro[1L], cromossomos)
        if (is.na(chr)) stop("Cromossomo inesperado: ", registro[1L])
        chr_lote[i] = chr
        dose_lote[[i]] = as.raw(dose)
        saida_lote[i] = paste(registro[c(seq_len(9L), col_amostras)], collapse = "\t")
      }
      for (grupo in split(which(chr_lote > 0L), chr_lote[chr_lote > 0L])) {
        chr = chr_lote[grupo[1L]]
        cat(saida_lote[grupo], file = arquivos[chr], sep = "\n", append = TRUE)
        con = file(binarios[chr], "ab")
        tryCatch(writeBin(do.call(c, dose_lote[grupo]), con), finally = close(con))
        contagem_chr[chr] = contagem_chr[chr] + length(grupo)
      }
    }
    contagens["QC"] = sum(contagem_chr)
    print(contagens)
    if (!contagens["QC"]) stop("Nenhum SNP passou no QC.")

    # 5. Reduza redundância por LD dentro de cada cromossomo.
    # Janela de 50 SNPs, avanço de 5 SNPs e limiar r² de 0,50.
    # Em um par acima do limiar, priorize o maior número de alelos menores.
    # Em empate, fica o SNP que aparece primeiro no VCF.
    # Cada processo recebe um cromossomo e executa os mesmos loops.
    
    podar_cromossomo = function(tarefa) {
      n = tarefa$n
      if (!n) return(logical())
      num_amostras = tarefa$num_amostras
      bruto = readBin(tarefa$arquivo, "raw", n = num_amostras * n)
      stopifnot(length(bruto) == num_amostras * n)
      dose = matrix(as.integer(bruto), nrow = num_amostras)
      rm(bruto)
      totais = colSums(dose)
      mac = pmin(totais, 2 * num_amostras - totais)
      ativos = rep(TRUE, n)

      for (inicio in seq.int(1L, n, by = 5L)) {
        fim = min(n, inicio + 49L)
        janela = seq.int(inicio, fim)
        alvos = janela[ativos[janela]]
        if (length(alvos)) {
          x = dose[, alvos, drop = FALSE]
          x = sweep(x, 2L, colMeans(x), "-")
          norma = sqrt(colSums(x^2))
          ativos[alvos[norma == 0]] = FALSE
          alvos = alvos[norma > 0]
          if (length(alvos) > 1L) {
            x = sweep(x[, norma > 0, drop = FALSE], 2L, norma[norma > 0], "/")
            r2 = crossprod(x)^2
            repeat {
              alterou = FALSE
              for (i in seq_len(length(alvos) - 1L)) {
                if (!ativos[alvos[i]]) next
                direita = seq.int(i + 1L, length(alvos))
                conflito = direita[ativos[alvos[direita]] &
                                    r2[i, direita] > 0.5 * (1 + 1e-12)]
                if (!length(conflito)) next
                j = conflito[1L]
                remover = if (mac[alvos[i]] < mac[alvos[j]]) i else j
                ativos[alvos[remover]] = FALSE
                alterou = TRUE
              }
              if (!alterou) break
            }
          }
        }
        if (fim == n) break
      }
      ativos
    }

    environment(podar_cromossomo) = baseenv()
    tarefas = vector("list", length(cromossomos))
    for (chr in seq_along(cromossomos)) {
      tarefas[[chr]] = list(n = contagem_chr[chr], arquivo = binarios[chr],
                            num_amostras = num_amostras)
    }
    stopifnot(length(n_nucleos) == 1L, is.finite(n_nucleos),
              n_nucleos >= 1, n_nucleos == floor(n_nucleos))
    processos = min(n_nucleos, sum(contagem_chr > 0L))
    message("Poda por LD com ", processos, " processo(s).")
    if (processos == 1L) {
      manter = lapply(tarefas, podar_cromossomo)
    } else {
      cluster = parallel::makePSOCKcluster(processos)
      manter = tryCatch(
        parallel::parLapplyLB(cluster, tarefas, podar_cromossomo),
        finally = parallel::stopCluster(cluster)
      )
    }
    apos_ld = vapply(manter, sum, integer(1L))
    print(data.frame(Cromossomo = cromossomos, Apos_QC = contagem_chr,
                     Apos_LD = apos_ld), row.names = FALSE)
    message("SNPs após LD: ", sum(apos_ld))
    if (!sum(apos_ld)) stop("Nenhum SNP passou na poda por LD.")

    # 6. Monte a matriz de linhagens com os SNPs retidos pela poda.
    # Cada linha é uma linhagem. As colunas guardam a dosagem do alelo ALT.
    # Aproveite as dosagens já lidas do VCF na etapa de controle de qualidade.
    n_snps = sum(apos_ld)
    matriz = matrix(NA_integer_, nrow = num_amostras, ncol = n_snps)
    marcadores = character(n_snps)
    deslocamento = 0L

    # Grave o VCF compactado e reúna os IDs na mesma ordem das dosagens.
    provis = tempfile(".g2f_", tmpdir = dirname(destino), fileext = ".vcf.gz")
    saida = gzfile(provis, "wt")
    tryCatch({
      writeLines(c(metadados, paste(c(cabecalho[seq_len(9L)], ids), collapse = "\t")), saida)
      for (i in which(apos_ld > 0L)) {
        bruto = readBin(binarios[i], "raw", n = num_amostras * contagem_chr[i])
        stopifnot(length(bruto) == num_amostras * contagem_chr[i])
        dose = matrix(as.integer(bruto), nrow = num_amostras)
        colunas_chr = deslocamento + seq_len(apos_ld[i])
        matriz[, colunas_chr] = dose[, manter[[i]], drop = FALSE]
        rm(bruto, dose)

        vcs = file(arquivos[i], "rt")
        tryCatch({
          inicio = 0L
          repeat {
            lote = readLines(vcs, n = 1000L, warn = FALSE)
            if (!length(lote)) break
            selecionar = manter[[i]][inicio + seq_along(lote)]
            retidos = lote[selecionar]
            writeLines(retidos, saida)
            if (length(retidos)) {
              # O ID do SNP é o terceiro campo de cada registro VCF.
              indices = deslocamento + seq_along(retidos)
              marcadores[indices] = sub("^[^\t]*\t[^\t]*\t([^\t]*).*", "\\1", retidos)
              deslocamento = deslocamento + length(retidos)
            }
            inicio = inicio + length(lote)
          }
          stopifnot(inicio == contagem_chr[i])
        }, finally = close(vcs))
      }
    }, finally = close(saida))
    stopifnot(deslocamento == n_snps, all(nzchar(marcadores) & marcadores != "."),
              !anyDuplicated(marcadores))
    dimnames(matriz) = list(ids, marcadores)
    if (!file.rename(provis, destino)) {
      stop("Não foi possível gravar o VCF final em: ", destino)
    }
    message("VCF salvo em: ", destino)

    invisible(list(hibridos = elegiveis, pendentes = pendentes, aliases = revisao,
                   parentais = ids, contagens = contagens, snps_apos_ld = sum(apos_ld),
                   arquivo = destino, M_linhagens = matriz))
  }, finally = {
    if (!is.null(entrada)) close(entrada)
    if (!is.null(provis)) unlink(provis)
    if (!is.null(pasta_temp)) unlink(pasta_temp, recursive = TRUE)
  })
})

# 7. Calcule a heterozigosidade nos SNPs retidos pela poda por LD.
# A dosagem 1 indica uma cópia de cada alelo, ou seja, um heterozigoto.
M_linhagens = preparacao$M_linhagens
preparacao$M_linhagens = NULL
heterozigotos = rowSums(M_linhagens == 1L)
n_snps = ncol(M_linhagens)
heterozigosidade = data.frame(
  Linhagem = rownames(M_linhagens),
  SNPs_validos = n_snps,
  SNPs_heterozigotos = heterozigotos,
  Heterozigosidade_pct = 100 * heterozigotos / n_snps,
  row.names = NULL
)
heterozigosidade = heterozigosidade[
  order(heterozigosidade$Heterozigosidade_pct, decreasing = TRUE), ,
  drop = FALSE
]
linhagens_revisar = heterozigosidade$Linhagem[heterozigosidade$SNPs_heterozigotos > 0L]
print(head(heterozigosidade, 20L), row.names = FALSE)

# 8. Retenha nas matrizes os parentais sem heterozigotos e seus cruzamentos.
# Um parental heterozigoto pode transmitir qualquer um dos seus dois alelos.
# O arquivo snp.clean registra os SNPs após a poda por LD, antes desta etapa.
cruzamentos = preparacao$hibridos
excluir = cruzamentos$DH1 %in% linhagens_revisar |
  cruzamentos$DH2 %in% linhagens_revisar
hibridos_revisar = cruzamentos[excluir, c("Pedigree", "DH1", "DH2"), drop = FALSE]
print(hibridos_revisar, row.names = FALSE)
message("Linhagens excluídas: ", length(linhagens_revisar))
message("Híbridos excluídos: ", nrow(hibridos_revisar))

cruzamentos = cruzamentos[!excluir, , drop = FALSE]
M_linhagens = M_linhagens[
  !rownames(M_linhagens) %in% linhagens_revisar, , drop = FALSE
]
if (!nrow(cruzamentos)) stop("Nenhum híbrido tem dois parentais homozigotos.")

# 9. Reconstrua a dosagem dos híbridos a partir dos dois parentais.
# Cada parental homozigoto transmite uma cópia do seu alelo.
# Por exemplo, um cruzamento 0/0 × 1/1 produz dosagem 1 no híbrido.
pai_1 = match(cruzamentos$DH1, rownames(M_linhagens))
pai_2 = match(cruzamentos$DH2, rownames(M_linhagens))
stopifnot(!anyNA(pai_1), !anyNA(pai_2))

M_hibridos = matrix(
  NA_integer_, nrow = nrow(cruzamentos), ncol = n_snps,
  dimnames = list(cruzamentos$Pedigree, colnames(M_linhagens))
)
for (j in seq_len(n_snps)) {
  M_hibridos[, j] = (M_linhagens[pai_1, j] + M_linhagens[pai_2, j]) %/% 2L
}

# Use Pedigree para relacionar os híbridos da matriz às médias ajustadas.
print(data.frame(
  Matriz = c("Linhagens", "Híbridos"),
  Genotipos = c(nrow(M_linhagens), nrow(M_hibridos)),
  SNPs = c(ncol(M_linhagens), ncol(M_hibridos))
), row.names = FALSE)

# Salve as matrizes com os IDs dos genótipos e dos SNPs.
saveRDS(M_hibridos, file.path("data", "marcadores_hibridos.rds"))
saveRDS(M_linhagens, file.path("data", "marcadores_linhagens.rds"))
cat("\nMarcadores dos híbridos\n")
print(M_hibridos[1:5, 1:8, drop = FALSE])

cat("\nMarcadores das linhagens\n")
print(M_linhagens[1:5, 1:8, drop = FALSE])

snp = "S1_201031"

c(
  B37 = M_linhagens["B37", snp],
  H95 = M_linhagens["H95", snp],
  "B37/H95" = M_hibridos["B37/H95", snp]
)

