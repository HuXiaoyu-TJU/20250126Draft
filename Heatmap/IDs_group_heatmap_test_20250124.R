library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(circlize)

data <- read.csv("IDs_group_heatmap_test.csv")
data_sorted <- data


expr_matrix <- as.matrix(data_sorted[, 14:31])
rownames(expr_matrix) <- data_sorted$ID

expr_matrix_norm <- sweep(expr_matrix, 
                          MARGIN = 2,  
                          STATS = colSums(expr_matrix, na.rm = TRUE),  
                          FUN = "/") * 100
expr_matrix_norm[expr_matrix_norm >= 3] <- 3

create_row_order <- function(data, expr_matrix) {

  row_na_counts <- rowSums(!is.na(expr_matrix))
  groups <- cut(row_na_counts, 
                breaks = c(5, 6, 7, 8, 9, 13, 18),
                labels = c("6 organs", "7 organs", "8 organs",
                           "9 organs", "10-13 organs", "14-18 organs"),
                include.lowest = TRUE)

  ordering_df <- data.frame(
    ID = rownames(expr_matrix),
    Group = groups,
    Glycan_subtype = factor(data$Glycan.subtype, 
                            levels = c("High-mannose", "Complex", "Hybrid")),
    Core_type = factor(data$Core.type,
                       levels = c("Common", "Core Fuc", "Bisect Core", "Bisect Core Fuc")),
    Antenna = factor(data$Antenna.count, 
                     levels = c("0", "1", "2", "3", "4")) 
  )
  
  for(g in unique(groups)) {
    cat("\n\n=== Group:", g, "===\n")
    group_data <- ordering_df[ordering_df$Group == g, ]
    ordered_group <- group_data[order(group_data$Glycan_subtype, 
                                      group_data$Core_type,
                                      group_data$Antenna
                                      ), ]
    
    cat("\nGlycan_subtype分布：\n")
    print(table(ordered_group$Glycan_subtype))
    
    cat("\nCore type分布：\n")
    print(table(ordered_group$Core_type))
    
    cat("\nAntenna count分布：\n")
    print(table(ordered_group$Antenna))

    cat("\n详细排序检查（前10行）：\n")
    print(head(data.frame(
      Glycan = ordered_group$Glycan_subtype,
      Core = ordered_group$Core_type,
      Antenna = ordered_group$Antenna
    ), 10))
  }

  ordered_df <- ordering_df %>%
    split(.$Group) %>%
    lapply(function(x) {
      x[order(x$Glycan_subtype, 
              x$Core_type, 
              x$Antenna
              ), ]
    }) %>%
    bind_rows()
  
  return(list(
    order = ordered_df$ID,
    groups = ordered_df$Group
  ))
}

row_order_result <- create_row_order(data_sorted, expr_matrix_norm)
row_order_new <- row_order_result$order
row_groups_new <- row_order_result$groups


expr_matrix_norm <- expr_matrix_norm[row_order_new, ]
data_sorted <- data_sorted[match(row_order_new, data_sorted$ID), ]

s_colors <- c("TRUE" = "purple", "FALSE" = "white")
f_colors <- c("TRUE" = "orange", "FALSE" = "white")
antenna_colors <- c(
  "0" = "#FFFFFF", "1" = "#BDD7E7", "2" = "#6BAED6",
  "3" = "#2171B5", "4" = "#08306B"
)

core_type_colors <- c(
  "Common" = "#4DAF4A",        
  "Core Fuc" = "#FF7F00",      
  "Bisect Core" = "#377EB8",   
  "Bisect Core Fuc" = "#E41A1C" 
)

glycan_type_colors <- c(
  'High-mannose' = '#B2DD8B',
  'Complex' = '#FDC06F',
  'Hybrid' = '#A6CEE5'
)

g_linkage_colors <- structure(
  colorRampPalette(c("#74C476", "#31A354", "#74ADD1", "#ABD9E9", 
                     "#E0F3F8", "#FEE090", "#FDBE85", "#A1D99B")
  )(length(unique(data_sorted$g.Linkage.))),
  names = unique(data_sorted$g.Linkage.)
)

accession_colors <- structure(
  colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                     "#FF7F00", "#FFFF33"))(length(unique(data_sorted$Accession))),
  names = unique(data_sorted$Accession)
)


value_breaks <- c(
  min(expr_matrix_norm, na.rm = TRUE),
  quantile(expr_matrix_norm, 0.25, na.rm = TRUE),
  median(expr_matrix_norm, na.rm = TRUE),
  quantile(expr_matrix_norm, 0.75, na.rm = TRUE),
  3
)

col_fun = colorRamp2(
  value_breaks,
  c("#4DA6FF", "#B3D9FF", "white", "#FFB3B3", "#FF6B6B")
)

ha = rowAnnotation(
  Accession = data_sorted$Accession,
  g_Linkage = data_sorted$g.Linkage.,
  glycan_type = data_sorted$Glycan.subtype,
  core_type = data_sorted$Core.type,
  Antenna = as.character(data_sorted$Antenna.count),
  F_status = data_sorted$F > 0,
  S_status = data_sorted$S > 0,
  col = list(
    Accession = accession_colors,
    g_Linkage = g_linkage_colors,
    glycan_type = glycan_type_colors,
    core_type = core_type_colors,
    Antenna = antenna_colors,
    F_status = f_colors,
    S_status = s_colors
  ),
  annotation_width = unit(rep(2, 7), "cm"),
  gap = unit(rep(2, 6), "mm"),
  annotation_legend_param = list(
    core_type = list(
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 16)
    ),
    glycan_type = list(
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 16)
    ),
    Antenna = list(
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 16)
    )
  ),
  show_annotation_name = TRUE,
  show_legend = c(Accession = FALSE,
                  g_Linkage = FALSE,
                  glycan_type = TRUE,
                  core_type = TRUE,
                  Antenna = TRUE,
                  F_status = FALSE,
                  S_status = FALSE),
  annotation_name_gp = gpar(fontsize = 16),
  gp = gpar(col = "grey", lwd = 0.2),
  border = FALSE
)

tissue_matrix <- t(expr_matrix_norm)
tissue_dist <- dist(tissue_matrix)
tissue_clust <- hclust(tissue_dist)
tissue_groups <- cutree(tissue_clust, k = 3)

ht = Heatmap(expr_matrix_norm,
             name = "Relative\nAbundance(%)",
             col = col_fun,
             na_col = "grey96",
             cluster_rows = FALSE,
             cluster_columns = TRUE,
             show_column_names = TRUE,
             show_row_names = FALSE,
             column_names_gp = gpar(fontsize = 16),
             row_names_gp = gpar(fontsize = 16),
             row_title_gp = gpar(fontsize = 16),
             column_title_gp = gpar(fontsize = 16),
             left_annotation = ha,
             column_title_rot = 90,
             row_title_rot = 0,
             border = FALSE,
             rect_gp = gpar(col = "grey", lwd = 0.2),
             column_split = tissue_groups,
             row_split = row_groups_new,
             row_gap = unit(2, "mm"),
             column_gap = unit(1, "mm"),
             heatmap_legend_param = list(
               legend_height = unit(4, "cm"),
               title_gp = gpar(fontsize = 16),
               labels_gp = gpar(fontsize = 16),
               at = c(0, 1, 2, 3),
               labels = c("0%", "1%", "2%", "≥3%")
             )
)

svg("all_groups_map.svg", 
    width = 297/25.4,
    height = 297/25.4)

draw(ht, 
     heatmap_legend_side = "right",
     annotation_legend_side = "bottom",
     show_heatmap_legend = TRUE,
     padding = unit(c(2, 2, 2, 2), "cm"))

dev.off()


grouped_data <- data.frame(
  ID = row_order_new,
  Group = row_groups_new,
  Glycan_subtype = data_sorted$Glycan.subtype,
  Core_type = data_sorted$Core.type,
  Antenna = data_sorted$Antenna.count
)
write.csv(grouped_data,
          file = "all_groups_data.csv",
          row.names = FALSE)
