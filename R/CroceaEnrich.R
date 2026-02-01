# 函数文档注释格式
#' Title: CroceaEnrich Function by liyin(liyin59375@gmail.com)
#' Description: This function performs GO and KEGG enrichment analysis.
#' @param file A data frame containing gene IDs.
#' @param prefix A character string used as the prefix for output files.
#' @export

CroceaEnrich <- function(file, orgdb, prefix = "CroceaEnrich_Result") {

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    install.packages("clusterProfiler", dependencies = TRUE)}
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse", dependencies = TRUE)}
  if (!requireNamespace("egg", quietly = TRUE)) {
    install.packages("egg", dependencies = TRUE)}
  
  library(clusterProfiler)
  library(tidyverse) 



  #01 读取基因列表文件并设置列名
  gene <- file
  colnames(gene) <- 'id'

  #02 获取当前时间
  current_datetime <- format(Sys.time(), "%Y%m%d%H%M")

  #03 创建输出目录
  output_dir <- paste0(prefix, "_Enrichment_", current_datetime)
  dir.create(output_dir, showWarnings = FALSE)

  #04 GO 富集分析
  go_results <- enrichGO(gene = gene$id,
                         keyType = "GID",
                         OrgDb = orgdb,  # 使用传递的 OrgDb 参数
                         ont = "all",
                         pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         readable = FALSE)

  go_results_df <- as.data.frame(go_results)
  go_results_df <- go_results_df[order(go_results_df$p.adjust),]
  go_filtered <- go_results_df[go_results_df$p.adjust < 0.05 & go_results_df$qvalue < 0.2,]
  go_top20 <- go_results_df[1:20,]

  ## 计算富集因子
  go_top20$BgRatio <- as.numeric(sub("/.*", "", go_top20$BgRatio))
  go_top20$GeneRatio <- as.numeric(sub("/.*", "", go_top20$GeneRatio))
  go_top20$Rich.factor <- go_top20$GeneRatio / go_top20$BgRatio

  ## 输出 GO 富集结果
  write.csv(go_results_df, file = paste0(output_dir, "/", prefix, "_GO_all_result.csv"), row.names = FALSE, quote = TRUE)
  write.csv(go_filtered, file = paste0(output_dir, "/", prefix, "_GO_filtered_result.csv"), row.names = FALSE, quote = TRUE)
  write.csv(go_top20, file = paste0(output_dir, "/", prefix, "_GO_top20.csv"), row.names = FALSE, quote = TRUE)

  #05 KEGG 富集分析
  pathway2name <- clusterProfiler:::kegg_list("pathway") %>%
    as.data.frame() %>%
    mutate(Pathway = gsub("map", "", from), Name = to) %>%
    mutate(ko = 'ko', Pathway = str_c(ko, Pathway)) %>%
    dplyr::select(Pathway, Name)
  pathway2gene <- AnnotationDbi::select(orgdb,
                                        keys = keys(orgdb),
                                        columns = c("Pathway", "Ko")) %>%
    na.omit() %>%
    dplyr::select(Pathway, GID)
  gene2des <- left_join(pathway2gene, pathway2name) %>%
    unique() %>%
    na.omit()

  kegg_results <- enricher(gene = gene$id,
                           TERM2GENE = gene2des[c("Pathway", "GID")],
                           TERM2NAME = gene2des[c("Pathway", "Name")],
                           pvalueCutoff = 1,
                           pAdjustMethod = 'BH',
                           qvalueCutoff = 1,
                           maxGSSize = 500)

  kegg_results_df <- as.data.frame(kegg_results)
  kegg_filtered <- kegg_results_df[kegg_results_df$p.adjust < 0.05 & kegg_results_df$qvalue < 0.2,]
  kegg_top20 <- kegg_results_df[1:20,]

  ## 计算富集因子
  kegg_top20$BgRatio <- as.numeric(sub("/.*", "", kegg_top20$BgRatio))
  kegg_top20$GeneRatio <- as.numeric(sub("/.*", "", kegg_top20$GeneRatio))
  kegg_top20$Rich.factor <- kegg_top20$GeneRatio / kegg_top20$BgRatio

  ## 输出 KEGG 富集结果
  write.csv(kegg_results_df, file = paste0(output_dir, "/", prefix, "_KEGG_all_result.csv"), row.names = FALSE, quote = TRUE)
  write.csv(kegg_filtered, file = paste0(output_dir, "/", prefix, "_KEGG_filtered_result.csv"), row.names = FALSE, quote = TRUE)
  write.csv(kegg_top20, file = paste0(output_dir, "/", prefix, "_KEGG_top20.csv"), row.names = FALSE, quote = TRUE)

  #06 GO 富集结果可视化 
  labels_go <- (levels(factor(go_top20$Description))[as.factor(go_top20$Description)])
  go_top20$number <- factor(rev(1:nrow(go_top20)))
  names(labels_go) <- rev(1:20)
  go_top20$shortname <- labels_go

  p_go <- ggplot(data = go_top20, aes(x = number, y = Rich.factor)) +
    geom_point(mapping = aes(size = Count, colour = -log10(qvalue), shape = ONTOLOGY)) +
    coord_flip() +
    scale_color_gradient(low = "#6BAED6",high = "#EE6A50")+
    scale_x_discrete(labels = labels_go) +
    labs(title = "Top 20 of GO enrichment", x = " ", y = "Rich factor",
         colour = "-log10(qvalue)", size = "Gene number") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(lineheight=.4, size=13, face = "bold", hjust = 0.5),
        legend.position = "right")
    
    go_fixed <- set_panel_size(
    p_go,
    width = unit(4, "cm"),
    height = unit(15, "cm"))

  ggsave(filename = paste0(output_dir, "/", prefix, "_GO_enrichment.pdf"),
         plot = go_fixed, height = 8)

  #07 KEGG 富集结果可视化 
  labels_kegg <- (levels(factor(kegg_top20$Description))[as.factor(kegg_top20$Description)])
  kegg_top20$number <- factor(rev(1:nrow(kegg_top20)))
  names(labels_kegg) <- rev(1:20)
  kegg_top20$shortname <- labels_kegg

  p_kegg <- ggplot(data = kegg_top20, aes(x = number, y = Rich.factor)) +
    geom_point(mapping = aes(size = Count, colour = -log10(qvalue))) +
    coord_flip() +
    scale_color_gradient(low = "#6BAED6",high = "#EE6A50")+
    scale_x_discrete(labels = labels_kegg) +
    labs(title = "Top 20 of KEGG enrichment", x = " ", y = "Rich factor",
         colour = "-log10(qvalue)", size = "Gene number") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(lineheight=.4, size=13, face = "bold", hjust = 0.5),
        legend.position = "right")
  
  kegg_fixed <- set_panel_size(
  p_kegg,
  width = unit(4, "cm"),
  height = unit(15, "cm"))

  ggsave(filename = paste0(output_dir, "/", prefix, "_KEGG_enrichment.pdf"),
         plot = kegg_fixed, height = 8)
}
