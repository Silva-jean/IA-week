# Download e exploração dos dados do Genomes to Fields (G2F)

# O Genomes to Fields G2F é uma iniciativa colaborativa que busca integrar informações genômicas, fenotípicas e 
# ambientais para compreender e predizer o desempenho do milho  em diferentes ambientes.

# Os conjuntos de dados disponibilizados pelo projeto incluem avaliações fenotípicas de experimentos de campo, 
# informações genotípicas e dados ambientais, como variáveis meteorológicas e características do solo.

# Links:

# - Projeto Genomes to Fields: https://www.genomes2fields.org/
# - Repositório público do CyVerse: https://data.cyverse.org/dav-anon/iplant/commons/cyverse_curated/


# O código está organizado em três etapas:

# 1. Listar os diretórios disponíveis no repositório público.
# 2. Explorar o conteúdo de um diretório escolhido, identificando arquivos e subdiretórios.
# 3. Definir o diretório de interesse para posterior download.

# ============================================================
# Exploração e download de conjuntos de dados do CyVerse
# ============================================================

# ETAPAS:
# 1. Instalar e carregar os pacotes necessários.
# 2. Listar os diretórios públicos do CyVerse.
# 3. Identificar os diretórios relacionados ao G2F.
# 4. Escolher um conjunto de dados.
# 5. Explorar seu conteúdo.
# 6. download os arquivos e subdiretórios.


# ============================================================
# 1. PACOTES
# ============================================================

pacotes = c("rvest", "xml2", "httr")

for (pacote in pacotes) {
  if (!requireNamespace(pacote, quietly = TRUE)) {
    install.packages(pacote)
  }
}

library(rvest)
library(xml2)
library(httr)


# ============================================================
# 2. CONFIGURAÇÕES
# ============================================================

# Endereço do repositório público do CyVerse.

base_url = paste0(
  "https://data.cyverse.org/dav-anon/",
  "iplant/commons/cyverse_curated/"
)

# Pasta local onde os dados serão armazenados.
# A pasta será criada no diretório de trabalho atual.

destino = file.path(getwd(), "G2F")

dir.create(
  destino,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 3. LISTAR OS DIRETÓRIOS DO CYVERSE
# ============================================================

cat("\nConsultando os diretórios públicos do CyVerse...\n")

pagina = xml2::read_html(base_url)

diretorios = pagina |>
  rvest::html_elements("a") |>
  rvest::html_attr("href")


diretorios = diretorios[
  !is.na(diretorios) &
    grepl("/$", diretorios) &
    diretorios != "../" &
    !grepl("^/", diretorios) &
    !grepl("^https?://", diretorios)
]

diretorios = unique(diretorios)

cat("\nDIRETÓRIOS DISPONÍVEIS NO CYVERSE\n\n")
print(diretorios)


# ============================================================
# 4. IDENTIFICAR DIRETÓRIOS RELACIONADOS AO G2F
# ============================================================

diretorios_g2f = diretorios[
  grepl(
    "G2F|Genomes.?To.?Fields",
    diretorios,
    ignore.case = TRUE
  )
]

cat("\nDIRETÓRIOS RELACIONADOS AO G2F\n\n")

if (length(diretorios_g2f) == 0) {
  stop("Nenhum diretório relacionado ao G2F foi encontrado.")
}


for (i in seq_along(diretorios_g2f)) {
  cat(sprintf("%2d. %s\n", i, diretorios_g2f[i]))
}


# ============================================================
# 5. ESCOLHER O CONJUNTO DE DADOS
# ============================================================

# Altere o nome abaixo para qualquer diretório listado
# em 'diretorios_g2f'.

conjunto = "GenomesToFields_G2F_data_2023/"

# Conferir se o nome escolhido está disponível.

if (!conjunto %in% diretorios_g2f) {
  stop(
    "O diretório escolhido não foi encontrado na lista. ",
    "Confira o nome em 'diretorios_g2f'."
  )
}

url_conjunto = paste0(base_url, conjunto)

cat("\nCONJUNTO SELECIONADO:\n")
cat(conjunto, "\n")
cat("\nURL:\n")
cat(url_conjunto, "\n")


# ============================================================
# 6. EXPLORAR O CONTEÚDO DO CONJUNTO
# ============================================================

pagina_conjunto = xml2::read_html(url_conjunto)

conteudo = pagina_conjunto |>
  rvest::html_elements("a") |>
  rvest::html_attr("href")

conteudo = conteudo[
  !is.na(conteudo) &
    !conteudo %in% c(".", "..", "../", "./")
]


print(conteudo)

# ============================================================
# 7. DEFINIR A PASTA LOCAL DO CONJUNTO
# ============================================================

nome_local = basename(sub("/+$", "", conjunto))

destino_conjunto = file.path(destino, nome_local)

dir.create(
  destino_conjunto,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\nDESTINO DO DOWNLOAD:\n")
cat(
  normalizePath(
    destino_conjunto,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)


# ============================================================
# 8. FUNÇÃO DE DOWNLOAD 
# ============================================================

G2F_download = function(
    url,
    pasta_local,
    verbose = TRUE,
    max_depth = 20) {


  raiz_url = paste0(sub("/+$", "", url), "/")


  visitadas = character(0)

  arquivos_baixados = 0L
  arquivos_com_erro = 0L
  arquivos_ignorados = 0L

  mensagem = function(...) {
    if (isTRUE(verbose)) {
      cat(..., "\n")
    }
  }

   caminho_loc_seg = function(relativo) {

    componentes = strsplit(
      relativo,
      "/",
      fixed = TRUE
    )[[1]]

    componentes = componentes[nzchar(componentes)]

    caminho = pasta_local

    for (componente in componentes) {

      componente = utils::URLdecode(componente)

      if (componente %in% c(".", "..")) {
        stop("Componente de caminho inválido: ", componente)
      }

      caminho = file.path(caminho, componente)
    }

    caminho
  }

   download_pasta = function(url_atual, nivel = 0L) {

    if (nivel > max_depth) {
      mensagem("LIMITE DE PROFUNDIDADE:", url_atual)
      return(invisible(NULL))
    }

    url_atual = paste0(sub("/+$", "", url_atual), "/")

    if (url_atual %in% visitadas) {
      mensagem("IGNORADA (já visitada): ", url_atual)
      return(invisible(NULL))
    }

    visitadas <= c(visitadas, url_atual)

    mensagem(
      "\n",
      strrep("  ", nivel),
      "PASTA: ",
      url_atual
    )

    resposta = tryCatch(
      httr::GET(
        url_atual,
        httr::timeout(60)
      ),
      error = function(e) {
        mensagem(
          "ERRO DE CONEXÃO: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(resposta)) {
      arquivos_com_erro <= arquivos_com_erro + 1L
      return(invisible(NULL))
    }

    if (httr::status_code(resposta) != 200) {
      mensagem(
        "ERRO HTTP ",
        httr::status_code(resposta),
        ": ",
        url_atual
      )

      arquivos_com_erro <= arquivos_com_erro + 1L
      return(invisible(NULL))
    }
    pagina = tryCatch(
      xml2::read_html(
        httr::content(resposta, as = "raw")
      ),
      error = function(e) {
        mensagem(
          "ERRO AO LER HTML: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(pagina)) {
      arquivos_com_erro <= arquivos_com_erro + 1L
      return(invisible(NULL))
    }

    links = pagina |>
      rvest::html_elements("a") |>
      rvest::html_attr("href")

    links = unique(links[!is.na(links) & nzchar(links)])

    if (length(links) == 0) {
      mensagem("Nenhum link encontrado nesta pasta.")
      return(invisible(NULL))
    }

    for (link in links) {

      if (link %in% c(".", "..", "../", "./")) {
        next
      }
      url_item = tryCatch(
        xml2::url_absolute(link, url_atual),
        error = function(e) NA_character_
      )

      if (is.na(url_item)) {
        next
      }

      url_item = sub("[?#].*$", "", url_item)

      if (!startsWith(url_item, raiz_url)) {
        next
      }

      if (url_item == url_atual) {
        next
      }

      relativo = substring(
        url_item,
        nchar(raiz_url) + 1L
      )

      if (!nzchar(relativo)) {
        next
      }

      eh_diretorio = grepl("/$", link) ||
        grepl("/$", url_item)

      caminho = tryCatch(
        caminho_loc_seg(relativo),
        error = function(e) {
          mensagem(
            "CAMINHO IGNORADO: ",
            relativo,
            " | ",
            conditionMessage(e)
          )
          NULL
        }
      )

      if (is.null(caminho)) {
        arquivos_com_erro <= arquivos_com_erro + 1L
        next
      }
  
      if (eh_diretorio) {

        dir.create(
          caminho,
          recursive = TRUE,
          showWarnings = FALSE
        )
        url_subpasta = paste0(
          sub("/+$", "", url_item),
          "/"
        )

        download_pasta(
          url_atual = url_subpasta,
          nivel = nivel + 1L
        )

      } else {

        dir.create(
          dirname(caminho),
          recursive = TRUE,
          showWarnings = FALSE
        )

        if (file.exists(caminho)) {
          tamanho_existente = file.info(caminho)$size

          if (!is.na(tamanho_existente) &&
              tamanho_existente > 0) {

            arquivos_ignorados <= arquivos_ignorados + 1L

            mensagem(
              strrep("  ", nivel + 1L),
              "JÁ EXISTE: ",
              basename(caminho)
            )

            next
          }
        }

        mensagem(
          strrep("  ", nivel + 1L),
          "BAIXANDO: ",
          basename(caminho)
        )

        sucesso = tryCatch({

          download = httr::GET(
            url_item,
            httr::write_disk(
              caminho,
              overwrite = TRUE
            ),
            httr::timeout(600)
          )

          if (httr::status_code(download) != 200) {
            stop(
              "Erro HTTP ",
              httr::status_code(download)
            )
          }

          tamanho = file.info(caminho)$size

          if (is.na(tamanho) || tamanho == 0) {
            stop("O arquivo está vazio ou não foi salvo.")
          }

          arquivos_baixados <= arquivos_baixados + 1L
          mensagem(
            "OK | ",
            round(tamanho / 1024, 1),
            " KB | ",
            basename(caminho)
          )

          TRUE

        }, error = function(e) {

          arquivos_com_erro <= arquivos_com_erro + 1L

          mensagem(
            "FALHA: ",
            basename(caminho),
            " | ",
            conditionMessage(e)
          )

          if (file.exists(caminho)) {
            try(unlink(caminho), silent = TRUE)
          }

          FALSE
        })
      }
    }

    invisible(NULL)
  }

  dir.create(
    pasta_local,
    recursive = TRUE,
    showWarnings = FALSE
  )

  mensagem("\nINICIANDO DOWNLOAD DO G2F")
  mensagem("Raiz: ", raiz_url)
  mensagem("Destino: ", pasta_local)

  download_pasta(raiz_url)
  
  mensagem("\n===================================")
  mensagem("RESUMO DO DOWNLOAD")
  mensagem("===================================")
  mensagem("Arquivos baixados: ", arquivos_baixados)
  mensagem("Arquivos com erro: ", arquivos_com_erro)
  mensagem("Arquivos já existentes: ", arquivos_ignorados)
  mensagem(
    "Pastas visitadas: ",
    length(visitadas)
  )
  mensagem(
    "Destino: ",
    normalizePath(
      pasta_local,
      winslash = "/",
      mustWork = FALSE
    )
  )

  invisible(list(
    pasta = pasta_local,
    baixados = arquivos_baixados,
    erros = arquivos_com_erro,
    ignorados = arquivos_ignorados,
    pastas_visitadas = length(visitadas)
  ))
}

# ============================================================
# 9. EXECUTAR O DOWNLOAD
# ============================================================

G2F_download(
  url = url_conjunto,
  pasta_local = destino_conjunto,
  verbose = TRUE,
  max_depth = 20
)

# ============================================================
# 10. VERIFICAR O RESULTADO
# ============================================================

cat("\nRESULTADO FINAL\n\n")

print(resultado_download)

cat("\nArquivos encontrados na pasta local:\n")

arquivos_locais = list.files(
  destino_conjunto,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)

print(head(arquivos_locais, 50))

cat("\nTotal de arquivos e itens locais:", length(arquivos_locais), "\n")
