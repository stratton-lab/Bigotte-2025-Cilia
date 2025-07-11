suppressPackageStartupMessages({
  library(clusterProfiler)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(forcats)
  library(ggplot2)
  library(SeuratObject)
  library(Seurat)
  library(SeuratWrappers)
  library(patchwork)
  library(tximeta)
  library(stringr)
  library(purrr)
  library(data.table)
  library(ggrepel)
})

GSE199460 <- function(data_dir) {
  so <- list(
    list(c("CTL1_total", "EAE1_total", "EAE2_total"), "GSE199460.1"),
    list(c("CTL2_total", "CTL3_total", "EAE3_total"), "GSE199460.2")
  ) |>
    purrr::map(
      purrr::partial(
        do.call,
        function(samples, batch) {
          purrr::map(
            samples,
            function(sample) {
              sname <- stringr::str_replace(
                sample,
                stringr::regex("_total$"),
                ""
              )

              file.path(data_dir, sample, "filtered_feature_bc_matrix") |>
                Read10X() |>
                CreateSeuratObject(project = "GSE199460") |>
                AddMetaData(NULL, "orig.ident") |>
                AddMetaData(sname, "sample") |>
                RenameCells(
                  add.cell.id = sname
                )
            }
          ) |>
            purrr::reduce(merge) |>
            JoinLayers() |>
            AddMetaData(batch, "batch")
        }
      )
    ) |>
    purrr::reduce(merge) |>
    JoinLayers()

  # fix broken cell identifiers in provided metadata
  meta <- data.table::fread(file.path(data_dir, "meta.csv.gz"))[, .(
    ident = V1,
    cell_type,
    condition
  )]

  # initialize cells hashset
  cells <- r2r::hashset()
  purrr::walk(Cells(so), purrr::partial(r2r::insert, cells))

  # id prefixes in cell identifiers were improperly set
  # specifically, CTL2 and CTL3 prefixes were
  # swapped for a small portion of cells.
  # we can detect these incorrect identifiers
  # and swap the sample prefix accordingly,
  # then verify that the new identifier is still valid.
  meta[, `:=`(
    ident = Vectorize(function(ident) {
      if (!r2r::has_key(cells, ident)) {
        split <- stringr::str_split_1(ident, stringr::fixed("_"))
        new <- if (split[[1]] == "CTL2") {
          sprintf("CTL3_%s", split[[2]])
        } else if (split[[1]] == "CTL3") {
          sprintf("CTL2_%s", split[[2]])
        } else {
          # broken ident is not prefixed with CTL2/CTL3, should never happen
          stop("error: unrecoverable ident")
        }
        if (!r2r::has_key(cells, new)) {
          # sanity check to make sure new ident is still valid
          stop("error: unrecoverable ident")
        }
        new
      } else {
        ident
      }
    })(ident)
  )]

  # drop cells not given annotations
  subset(so, cells = meta[, ident]) |>
    AddMetaData(data.frame(meta, row.names = "ident")) |>
    SetIdent(value = "cell_type") |>
    AddMetaData("GSE199460", "source")
}

GSE254863 <- function(data_dir) {
  c("EAE4p", "EAE4n") |>
    purrr::map(
      function(sample) {
        file.path(data_dir, sample) |>
          Read10X() |>
          CreateSeuratObject(project = "GSE254863") |>
          AddMetaData(NULL, "orig.ident") |>
          AddMetaData(sample, "sample") |>
          RenameCells(add.cell.id = sample)
      }
    ) |>
    purrr::reduce(merge) |>
    JoinLayers() |>
    AddMetaData("GSE254863.1", "batch") |>
    AddMetaData("GSE254863", "source") |>
    AddMetaData("E", "condition") |>
    AddMetaData(
      file.path(data_dir, "meta.csv.gz") |>
        data.table::fread(header = TRUE) |>
        _[, `:=`(
          V1 = stringr::str_replace(V1, "Neg_", "EAE4n_") |>
            stringr::str_replace("Pos_", "EAE4p_")
        )][] |>
        data.frame(row.names = "V1")
    )
}


term_is_offspring <- function(terms, filter.by) {
  setDT(
    AnnotationDbi::select(
      GO.db::GO.db,
      keys = filter.by,
      columns = "ONTOLOGY"
    )
  )[
    # NA ontology indicates GO term not in GO.db, ignore such cases
    # (most likely due to desync between GO.db and clusterProfiler)
    !is.na(ONTOLOGY),
    .(
      ID = switch(ONTOLOGY, "BP" = GO.db::GOBPOFFSPRING, "CC" = GO.db::GOCCOFFSPRING, "MF" = GO.db::GOMFOFFSPRING, )[GOID] |>
        as.list() |>
        reduce(c)
    ),
    by = .(ONTOLOGY),
  ][, terms %chin% ID]
}

gene_within_goset <- function(genes, filter_by, key.type, org.db) {
  setDT(
    AnnotationDbi::select(
      org.db,
      keys = genes,
      keytype = key.type,
      columns = c("GOALL")
    )
  )[
    GOALL %chin% filter_by,
    genes %chin% get(key.type)
  ]
}

mk_cilia_plots <- function(m.dt, log2FC.col, pval.col, org.db, n.genes = 50) {
  m.dt <- copy(m.dt)[, `:=`(
    is.ciliary = gene_within_goset(
      SYMBOL,
      "GO:0044782",
      "SYMBOL",
      org.db
    )
  )][]

  list(
    "pl.1" = ggplot(m.dt[order(is.ciliary)]) +
      aes(
        x = get(log2FC.col),
        y = -log10(get(pval.col)),
        colour = fct_relevel(
          fcase(
            is.ciliary & get(pval.col) <= 0.05 & get(log2FC.col) > 0,
            "up - ciliary",
            !is.ciliary & get(pval.col) <= 0.05 & get(log2FC.col) > 0,
            "up",
            is.ciliary & get(pval.col) <= 0.05 & get(log2FC.col) < 0,
            "down - ciliary",
            !is.ciliary & get(pval.col) <= 0.05 & get(log2FC.col) < 0,
            "down",
            default = "n.s."
          ),
          "down - ciliary",
          "down",
          "up - ciliary",
          "up",
          "n.s."
        )
      ) +
      geom_point(size = 1) +
      scale_alpha_identity() +
      guides(
        colour = guide_legend(
          override.aes = list(shape = 16, size = 3),
          nrow = 2
        )
      ) +
      labs(colour = NULL, x = "log2 fold change", y = "-log10(pval)") +
      scale_color_brewer(palette = "Set1") +
      theme(
        legend.position = "bottom",
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 16)
      ),
    "pl.2" = m.dt[is.ciliary & get(pval.col) < 0.05][
      data.table::last(order(abs(get(log2FC.col))), n.genes)
    ] |>
      ggplot() +
      aes(
        y = fct_reorder(SYMBOL, order(abs(get(log2FC.col)))),
        x = get(log2FC.col)
      ) +
      geom_col() +
      labs(y = NULL, x = "log2 fold change") +
      theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank()
      ),
    "pl.3" = rbind(
      as.data.table(
        enrichGO(
          m.dt[get(pval.col) < 0.05 & get(log2FC.col) > 0, SYMBOL],
          OrgDb = org.db,
          keyType = "SYMBOL",
          ont = "BP",
          minGSSize = 2,
          maxGSSize = 5000
        )@result
      )[, `:=`(regulated = "up")],
      as.data.table(
        enrichGO(
          m.dt[get(pval.col) < 0.05 & get(log2FC.col) < 0, SYMBOL],
          OrgDb = org.db,
          keyType = "SYMBOL",
          ont = "BP",
          minGSSize = 2,
          maxGSSize = 5000
        )@result
      )[, `:=`(regulated = "down")]
    )[
      qvalue <= 0.05
    ][
      term_is_offspring(ID, "GO:0044782")
    ][
      order(qvalue)
    ] |>
      ggplot() +
      aes(x = Count, y = fct_reorder(Description, Count), fill = qvalue) +
      geom_col() +
      scale_fill_viridis_c(option = "turbo") +
      labs(x = "count", y = NULL) +
      facet_wrap(vars(regulated), nrow = 2, scales = "free_y") +
      theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        axis.text.y = element_text(size = 16)
      ),
    "dt" = m.dt
  )
}

mk_dr_plot <- function(pm.dt) {
  wrap_plots(
    ggplot(pm.dt) +
      aes(x = i.PaCMAP_1, y = i.PaCMAP_2) +
      labs(shape = NULL, colour = NULL) +
      geom_point(
        aes(
          colour = cell_type,
          shape = fifelse(is.na(PaCMAP_1), "dropped", "kept")
        ),
        size = 0.4
      ) +
      geom_text_repel(
        aes(label = cell_type),
        pm.dt[,
          lapply(.SD, median),
          by = .(cell_type),
          .SDcols = c("i.PaCMAP_1", "i.PaCMAP_2")
        ]
      ) +
      theme_void() +
      theme(legend.position = "left") +
      guides(
        colour = "none",
        shape = guide_legend(override.aes = list(size = 3))
      ),
    ggplot(pm.dt[!is.na(PaCMAP_1)]) +
      aes(
        colour = fifelse(condition == "E", "EAE", "CTL"),
        x = PaCMAP_1,
        y = PaCMAP_2
      ) +
      geom_point(size = 0.4) +
      scale_colour_manual(values = c("CTL" = "blue", "EAE" = "red")) +
      guides(colour = guide_legend(override.aes = list(size = 3))) +
      labs(colour = NULL) +
      theme_void(),
    design = "
      AAAA###
      AAAA#BB
      AAAA#BB
      AAAA###
    "
  )
}

if (!interactive()) {
  # ===================================================
  # generate mouse EAE sc and rat in vitro bulkseq data
  # ===================================================

  so.eae <- merge(
    GSE199460("data/GSE199460"),
    GSE254863("data/GSE254863")
  ) |>
    JoinLayers()

  so.ec.eae <- suppressWarnings({
    subset(so.eae, cell_type == "Ependymal Cells")
  }) |>
    SetIdent(value = "source") |>
    subset(downsample = 500) |>
    NormalizeData()
  makeLinkedTxome(
    indexDir = "data/bulkseq/index",
    source = "RefSeq",
    organism = "Rattus norvegicus",
    release = "latest",
    genome = "GRCr8",
    fasta = c(
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_rna.fna.gz",
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.fna.gz"
    ),
    gtf = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.gff.gz"
  )

  m.dt.eae <- suppressMessages({
    FindMarkers(so.ec.eae, "E", "C", group.by = "condition") |>
      as.data.table(keep.rownames = "SYMBOL")
  })

  m.dt.bulk <- fread("data/bulkseq/meta.csv")[, `:=`(
    files = map_vec(names, function(n) {
      str_subset(
        list.files("data/bulkseq/salmon", "quant\\.sf", recursive = TRUE, full.names = TRUE),
        pattern = fixed(n)
      )
    }),
    group = fct_relevel(group, "nt")
  )][group %in% c("ms", "nt")] |>
    tximeta(type = "salmon") |>
    summarizeToGene() |>
    DESeqDataSet(~group) |>
    DESeq() |>
    results() |>
    as.data.table(keep.rownames = "SYMBOL")

  pacmap.dt <- as.data.table(
    so.ec.eae |>
      FindVariableFeatures() |>
      ScaleData() |>
      RunPCA() |>
      RunPaCMAP(dims = 1:30) |>
      FetchData(c("PaCMAP_1", "PaCMAP_2")),
    keep.rownames = "cell"
  )[
    as.data.table(
      so.eae |>
        NormalizeData() |>
        FindVariableFeatures() |>
        ScaleData() |>
        RunPCA() |>
        RunPaCMAP(dims = 1:30) |>
        FetchData(
          c("PaCMAP_1", "PaCMAP_2", "cell_type", "condition")
        ),
      keep.rownames = "cell"
    )[, `:=`(
      cell_type = fcoalesce(cell_type, "Unknown")
    )],
    on = .(cell)
  ]

  # ===================
  # generate/save plots
  # ===================

  eae.cilia.plots <- mk_cilia_plots(m.dt.eae, "avg_log2FC", "p_val_adj", org.Mm.eg.db)
  ggsave(
    "out/eae-out1.tiff",
    eae.cilia.plots[["pl.1"]],
    w = 1500,
    h = 1500,
    units = "px"
  )
  ggsave(
    "out/eae-out2.tiff",
    eae.cilia.plots[["pl.2"]],
    w = 2000,
    h = 2000,
    units = "px"
  )
  ggsave(
    "out/eae-out3.tiff",
    eae.cilia.plots[["pl.3"]],
    w = 2000,
    h = 1500,
    units = "px"
  )

  bulk.cilia.plots <- mk_cilia_plots(m.dt.bulk[!is.na(padj)], "log2FoldChange", "padj", org.Rn.eg.db)
  ggsave(
    "out/bulk-out1.tiff",
    bulk.cilia.plots[["pl.1"]],
    w = 1500,
    h = 1500,
    units = "px"
  )
  ggsave(
    "out/bulk-out2.tiff",
    bulk.cilia.plots[["pl.2"]],
    w = 2000,
    h = 2000,
    units = "px"
  )
  ggsave(
    "out/bulk-out3.tiff",
    bulk.cilia.plots[["pl.3"]],
    w = 2000,
    h = 1500,
    units = "px"
  )

  maid.cilia.plots <- mk_cilia_plots(
    setnames(fread("out/ependyma_ciliary_cmp.csv"), "V1", "SYMBOL"),
    "avg_log2FC",
    "p_val_adj",
    org.Hs.eg.db
  )
  ggsave(
    "out/human-out1.tiff",
    maid.cilia.plots[["pl.1"]],
    w = 1500,
    h = 1500,
    units = "px"
  )
  ggsave(
    "out/human-out2.tiff",
    maid.cilia.plots[["pl.2"]],
    w = 2000,
    h = 2000,
    units = "px"
  )
  ggsave(
    "out/human-out3.tiff",
    maid.cilia.plots[["pl.3"]],
    w = 2000,
    h = 1500,
    units = "px"
  )

  pl.4 <- mk_dr_plot(pacmap.dt)
  ggsave("out/sc_pacmap.tiff", pl.4, w = 4000, h = 2000, units = "px")
}
