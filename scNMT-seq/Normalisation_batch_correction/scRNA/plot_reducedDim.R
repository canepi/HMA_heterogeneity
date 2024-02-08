plot_reducedDim <- function(sce, reducedDim = 'PCA', comps = c(1,2), colBy = 'Batch', facet_by = NULL, cols, axis.title = 'PC', ...)
{
    df <- sce@int_colData$reducedDims[[reducedDim]]
    df <- df[,comps]
    df <- as.data.frame(df)
    colnames(df) <- paste0(axis.title, comps)
    pcs <- colnames(df)
    df$group <- colData(sce)[,colBy]
    df$facet <- colData(sce)[,facet_by]
    cols <- cols[unique(as.character(df$group))]

    p <- ggplot(df, aes_string(pcs[1], pcs[2])) + geom_point(aes(col = group),...) +
        scale_colour_manual(values = cols) +
        theme_classic() + labs(col = colBy)
    if (!is.null(facet_by))
    {
        p <- p + facet_wrap(.~facet, scales = 'free')
    }
    p
}
