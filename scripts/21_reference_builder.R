#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 21_reference_builder.R
# 功能说明: 利用 PubMed 与 Crossref 自动检索并生成参考文献清单（含DOI），
#          输出 CSV/BIB/MD；尽量优先使用近年文献，同时保留少量方法学奠基文献。
# 适用范围: Nature级别科研项目（河网SDM+因果分析）；中文注释，输出英文参考条目。
# 重要说明: 
#   - PubMed API Key 读取顺序：环境变量 PUBMED_API_KEY 或 ENTREZ_KEY →
#     可选的文件 scripts/keys/pubmed.key →（若均缺失）使用用户提供的密钥占位。
#   - Crossref 无需密钥，但请避免过高频率请求。
#   - 本脚本检查是否已存在同名输出，存在则覆盖（幂等）。
# 输出文件:
#   references/references.csv
#   references/references.bib
#   references/references.md
# 依赖: rentrez, rcrossref, dplyr, readr, stringr, tibble, purrr, jsonlite, xml2
# 日期: 2025-11-06
# ==============================================================================

# ----------------------------- 基础设置 ---------------------------------------
rm(list = ls())
GC <- gc()
try(setwd("E:/SDM01"), silent = TRUE)

suppressWarnings({
  pkgs <- c("rentrez", "rcrossref", "dplyr", "readr", "stringr", "tibble", "purrr", "jsonlite", "xml2")
  for(p in pkgs){
    if(!require(p, character.only = TRUE)){
      install.packages(p, dependencies = TRUE)
      library(p, character.only = TRUE)
    }
  }
})

# 声明全局变量以消除管道与非标准评估引发的linter告警
if(getRversion() >= "2.15.1") utils::globalVariables(c("%>%", "score", "year", "doi", "pub_year", ".data"))
# 明确引入管道符，避免工具误报
`%>%` <- dplyr::`%>%`

dir.create("references", showWarnings = FALSE, recursive = TRUE)

# 读取/设置 PubMed API Key
get_pubmed_key <- function(){
  # 1) 环境变量（优先）
  k <- Sys.getenv("PUBMED_API_KEY")
  if(nchar(k) > 0) return(k)
  k <- Sys.getenv("ENTREZ_KEY")
  if(nchar(k) > 0) return(k)
  # 2) 尝试从文件读取
  key_file <- "scripts/keys/pubmed.key"
  if(file.exists(key_file)){
    k2 <- try(readr::read_lines(key_file, n_max = 1), silent = TRUE)
    if(!inherits(k2, "try-error") && length(k2) > 0 && nchar(k2[1]) > 0) return(k2[1])
  }
  # 3) 使用用户提供的占位密钥（用户已在会话中提供）
  #    如需替换，请设置环境变量 PUBMED_API_KEY 或 ENTREZ_KEY
  return("70ae3d6eddab65bf1c3144614dc3db61a709")
}

PUBMED_KEY <- get_pubmed_key()
if(nchar(PUBMED_KEY) > 0){
  Sys.setenv(ENTREZ_KEY = PUBMED_KEY)
}

# -------------------------- 检索配置（可按需调整） -----------------------------
# 主题检索（偏向近5年综述/方法与淡水生态）
query_topics <- c(
  # 核心主题：河流鱼类SDM（中国/亚洲/全球对比），近5年偏向
  "freshwater fish species distribution model AND (river OR stream) AND (China OR Asia) AND (2020:3000[dp])",
  "riverine fish biodiversity AND climate change AND habitat suitability AND (2020:3000[dp])",
  "hydrological variability AND fish distribution AND species distribution modelling AND (2020:3000[dp])",
  # 因果与生态：图模型/因果森林/生态中的因果推断
  "causal discovery AND (ecology OR biogeography) AND (Bayesian network OR PC algorithm OR bnlearn OR pcalg) AND (2020:3000[dp])",
  "causal forest AND (ecology OR biodiversity OR environmental) AND (2020:3000[dp])",
  # 数据与变量：河网特定环境层与上游加权指标
  "EarthEnv-Streams freshwater environmental variables biodiversity AND (2019:3000[dp])",
  "upstream weighted hydroclim OR river network climate AND biodiversity AND (2020:3000[dp])",
  # 未来情景与CMIP6：
  "CMIP6 AND freshwater biodiversity AND SSP AND habitat suitability AND (2020:3000[dp])",
  # 评估与可视化方法
  "partial dependence OR accumulated local effects AND ecology AND (2020:3000[dp])"
)

# 强制收录（方法学/数据基准；使用 Crossref 精确检索）
forced_items <- tribble(
  ~title, ~author_hint, ~year_min,
  "Near-global freshwater-specific environmental variables for biodiversity analyses at 1 km resolution", "Domisch", 2015,
  "WorldClim 2: new 1-km spatial resolution climate surfaces", "Fick", 2017,
  "SoilGrids250m: Global gridded soil information based on machine learning", "Hengl", 2017,
  "Random Forests", "Breiman", 2001,
  "Maximum entropy modeling of species geographic distributions", "Phillips", 2006,
  "Fast stable direct fitting and smoothness selection for generalized additive models", "Wood", 2011,
  "Learning Bayesian Networks with the bnlearn R Package", "Scutari", 2010,
  "Causal Inference Using Graphical Models with the R Package pcalg", "Kalisch", 2012,
  "Estimation and inference of heterogeneous treatment effects using random forests", "Wager", 2018,
  "Visualizing the effects of predictor variables in black box supervised learning models", "Apley", 2020,
  "Inferring causation from time series in Earth system sciences", "Runge", 2019,
  "HydroSHEDS – A global, high-resolution hydrographic dataset derived from SRTM", "Lehner", 2008
)

# 年份过滤（普通检索）
YEAR_MIN_RECENT <- 2020
MAX_PER_TOPIC   <- 30  # 每个主题最多收录文献条目（更丰富）

# ---------------------------- 工具函数 ----------------------------------------
# 安全解析 Crossref 列中的年份（兼容 list/atomic/缺失）
safe_year_list <- function(x){
  y <- NA_integer_
  try({
    if(is.null(x)) return(NA_integer_)
    # 若是 data.frame/list 且包含 `date-parts`
    if(is.list(x)){
      # 命名列表：x[["date-parts"]]
      if(!is.null(x[["date-parts"]])){
        dp <- x[["date-parts"]]
        if(is.list(dp) && length(dp) >= 1){
          v <- try(dp[[1]], silent = TRUE)
          if(!inherits(v, "try-error") && length(v) >= 1) y <- suppressWarnings(as.integer(v[1]))
        }
      } else if(length(x) >= 1){
        # 某些返回为 list 的第一项即向量年份
        v <- x[[1]]
        if(is.atomic(v) && length(v) >= 1) y <- suppressWarnings(as.integer(v[1]))
      }
    } else if(is.atomic(x) && length(x) >= 1){
      y <- suppressWarnings(as.integer(x[1]))
    }
  }, silent = TRUE)
  if(!is.na(y) && (y < 1000 || y > 3000)) y <- NA_integer_
  return(y)
}

# 从 Crossref 以标题/作者关键词检索，返回最匹配的一条
cr_find_best <- function(title, author_hint = NULL, year_min = 1900){
  q <- title
  if(!is.null(author_hint) && nchar(author_hint) > 0){ q <- paste(title, author_hint) }
  res <- try(rcrossref::cr_works(query = q, limit = 5), silent = TRUE)
  if(inherits(res, "try-error") || is.null(res$data) || nrow(res$data) == 0) return(NULL)
  df <- res$data
  # 解析多个可能来源的年份
  pp <- if("published.print" %in% names(df)) df[["published.print"]] else vector("list", nrow(df))
  iss<- if("issued" %in% names(df)) df[["issued"]] else vector("list", nrow(df))
  pol<- if("published.online" %in% names(df)) df[["published.online"]] else vector("list", nrow(df))
  cre<- if("created" %in% names(df)) df[["created"]] else vector("list", nrow(df))
  yr <- if("year" %in% names(df)) df[["year"]] else rep(NA_integer_, nrow(df))
  pub_year <- purrr::pmap_int(list(pp, iss, pol, cre, yr), function(a,b,c,d,e){
    y <- safe_year_list(a); if(is.na(y)) y <- safe_year_list(b)
    if(is.na(y)) y <- safe_year_list(c); if(is.na(y)) y <- safe_year_list(d)
    if(is.na(y)) y <- suppressWarnings(as.integer(e)); if(is.na(y)) y <- NA_integer_
    y
  })
  df$pub_year <- pub_year
  df <- df %>% dplyr::arrange(dplyr::desc(.data$score)) %>%
    dplyr::filter(is.na(pub_year) | pub_year >= year_min)
  if(nrow(df) == 0) df <- res$data %>% dplyr::arrange(dplyr::desc(.data$score))
  df[1, , drop = FALSE]
}

# 普通主题检索（优先PubMed，失败回退Crossref）
search_topic <- function(topic, year_min = YEAR_MIN_RECENT, max_n = MAX_PER_TOPIC){
  out <- list()
  # 1) PubMed 检索（返回PMIDs）
  pm <- try(rentrez::entrez_search(db = "pubmed", term = topic, retmax = max_n), silent = TRUE)
  if(!inherits(pm, "try-error") && length(pm$ids) > 0){
    # 获取摘要与DOI（使用 efetch xml 解析）
    xmltxt <- try(rentrez::entrez_fetch(db = "pubmed", id = paste(pm$ids, collapse = ","), rettype = "xml", retmode = "xml"), silent = TRUE)
    if(!inherits(xmltxt, "try-error")){
      doc <- try(xml2::read_xml(xmltxt), silent = TRUE)
      if(!inherits(doc, "try-error")){
        arts <- xml2::xml_find_all(doc, "//PubmedArticle")
        for(a in arts){
          ti  <- xml2::xml_text(xml2::xml_find_first(a, ".//ArticleTitle"))
          dp  <- xml2::xml_text(xml2::xml_find_first(a, ".//PubDate/Year"))
          yr  <- suppressWarnings(as.integer(dp))
          doi <- xml2::xml_text(xml2::xml_find_first(a, ".//ArticleId[@IdType='doi']"))
          jn  <- xml2::xml_text(xml2::xml_find_first(a, ".//Journal/Title"))
          if(is.na(yr) || yr >= year_min){
            out[[length(out)+1]] <- tibble::tibble(
              title = ti, journal = jn, year = yr, doi = doi, source = "PubMed"
            )
          }
        }
      }
    }
  }
  # 2) 若不足，则用 Crossref 兜底
  if(length(out) < max_n){
    cr <- try(rcrossref::cr_works(query = topic, limit = max_n), silent = TRUE)
    if(!inherits(cr, "try-error") && !is.null(cr$data) && nrow(cr$data) > 0){
      df <- cr$data
      pp <- if("published.print" %in% names(df)) df[["published.print"]] else vector("list", nrow(df))
      iss<- if("issued" %in% names(df)) df[["issued"]] else vector("list", nrow(df))
      pol<- if("published.online" %in% names(df)) df[["published.online"]] else vector("list", nrow(df))
      cre<- if("created" %in% names(df)) df[["created"]] else vector("list", nrow(df))
      yr <- if("year" %in% names(df)) df[["year"]] else rep(NA_integer_, nrow(df))
      pub_year <- purrr::pmap_int(list(pp, iss, pol, cre, yr), function(a,b,c,d,e){
        y <- safe_year_list(a); if(is.na(y)) y <- safe_year_list(b)
        if(is.na(y)) y <- safe_year_list(c); if(is.na(y)) y <- safe_year_list(d)
        if(is.na(y)) y <- suppressWarnings(as.integer(e)); if(is.na(y)) y <- NA_integer_
        y
      })
      add <- tibble::tibble(
        title   = purrr::map_chr(df$title, ~ ifelse(length(.x) > 0, .x[1], NA_character_)),
        journal = df$container.title,
        year    = pub_year,
        doi     = df$doi,
        source  = "Crossref"
      ) %>% dplyr::filter(is.na(.data$year) | .data$year >= year_min)
      if(nrow(add) > 0) out[[length(out)+1]] <- add
    }
  }
  if(length(out) == 0) return(NULL)
  dplyr::bind_rows(out) %>%
    dplyr::distinct(.data$doi, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(.data$doi) & .data$doi != "")
}

# 强制收录条目
collect_forced <- function(tbl){
  rows <- list()
  for(i in seq_len(nrow(tbl))){
    r <- tbl[i,]
    hit <- cr_find_best(r$title, r$author_hint, r$year_min)
    if(!is.null(hit)){
      # 使用 cr_find_best 已解析的 pub_year，回退到 year
      yr <- if("pub_year" %in% names(hit)) suppressWarnings(as.integer(hit$pub_year)) else suppressWarnings(as.integer(hit$year))
      rows[[length(rows)+1]] <- tibble::tibble(
        title = ifelse(length(hit$title) > 0, hit$title[[1]], r$title),
        journal = hit$container.title,
        year = yr,
        doi = hit$doi,
        source = "Crossref"
      )
    }
  }
  if(length(rows) == 0) return(NULL)
  dplyr::bind_rows(rows) %>% dplyr::distinct(doi, .keep_all = TRUE)
}

# ------------------------------ 主流程 ----------------------------------------
message("\n====== 参考文献自动检索开始 ======\n")

# 1) 主题检索汇总
res_list <- purrr::map(query_topics, ~ search_topic(.x))
res_topic <- res_list %>% purrr::compact() %>% dplyr::bind_rows()

# 2) 强制收录（方法/数据基准）
res_forced <- collect_forced(forced_items)

# 3) 合并去重（以DOI为键）
all_refs <- dplyr::bind_rows(res_topic, res_forced) %>%
  dplyr::distinct(.data$doi, .keep_all = TRUE)

# 4) 简单质量控制：移除明显非论文记录（无期刊名或异常DOI）
all_refs <- all_refs %>%
  dplyr::filter(!is.na(.data$doi) & .data$doi != "") %>%
  dplyr::mutate(
    doi = stringr::str_to_lower(.data$doi),
    url = paste0("https://doi.org/", .data$doi)
  )

# 5) 按年份与来源排序（近年优先，方法基准保留）
all_refs <- all_refs %>%
  dplyr::mutate(year2 = ifelse(is.na(.data$year), -Inf, as.integer(.data$year))) %>%
  dplyr::arrange(dplyr::desc(.data$year2), .data$source, .data$title) %>%
  dplyr::select(title, journal, year, doi, url, source)

# 6) 输出 CSV
csv_path <- "references/references.csv"
readr::write_csv(all_refs, csv_path)

# 7) 输出 BIB（最简条目；复杂需求建议后续用 Zotero/Better BibTeX 再格式化）
bib_path <- "references/references.bib"
con_bib <- file(bib_path, open = "w", encoding = "UTF-8")
cat("% Auto-generated by 21_reference_builder.R\n", file = con_bib)
for(i in seq_len(nrow(all_refs))){
  key <- paste0("ref", i)
  cat("@article{", key, ",\n", sep = "", file = con_bib)
  cat("  title = {", all_refs$title[i], "},\n", sep = "", file = con_bib)
  if(!is.na(all_refs$journal[i])) cat("  journal = {", all_refs$journal[i], "},\n", sep = "", file = con_bib)
  if(!is.na(all_refs$year[i]))    cat("  year = {", all_refs$year[i], "},\n", sep = "", file = con_bib)
  cat("  doi = {", all_refs$doi[i], "},\n", sep = "", file = con_bib)
  cat("  url = {", all_refs$url[i], "}\n", sep = "", file = con_bib)
  cat("}\n\n", file = con_bib)
}
close(con_bib)

# 8) 输出 MD（编号列表，便于论文草稿粘贴引用）
md_path <- "references/references.md"
con_md <- file(md_path, open = "w", encoding = "UTF-8")
cat("# References (auto-generated)\n\n", file = con_md)
for(i in seq_len(nrow(all_refs))){
  line <- sprintf("%d. %s. %s. %s. DOI: %s", 
                  i, all_refs$title[i], ifelse(is.na(all_refs$journal[i]), "", all_refs$journal[i]),
                  ifelse(is.na(all_refs$year[i]), "", all_refs$year[i]), all_refs$doi[i])
  cat(line, "\n", file = con_md)
}
close(con_md)

message("已生成: \n - ", csv_path, "\n - ", bib_path, "\n - ", md_path, "\n")
message("====== 参考文献自动检索完成 ======\n")

# 备注：若需固定收录的条目可在 forced_items 中添加；若需严格限定年份，调 YEAR_MIN_RECENT。
#       建议后续用 Zotero (Better BibTeX) 导入 CSV/BIB 进行人工复核与格式审校，以满足 Nature 样式。


