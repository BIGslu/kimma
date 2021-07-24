
#' Title
#'
#' @param dat DGEList or EList object with expression data (counts, E), sample metadata (samples, targets), and gene annotation (genes)
#' @param counts If dat not provided. Data frame of gene expression. Genes are rows, samples are columns. geneID must be rownames or first column
#' @param meta If dat not provided. Data frame of sample meta data with samples as rows.
#' @param genes Optional if dat not provided. Data frame of gene annotation with genes as rows.
#' @param fdr Optional. Model results output by kmFit( )
#' @param model.name Optional. Character string of model name to include in fdr results
#' @param libraryID Character string of variable name to use when combining expression and sample data
#' @param geneID Character string of variable name to use when combining expression and gene data
#' @param subset.genes Optional. Character vector of genes to plot. Must match names in geneID column. If not provided, all genes are plotted.
#' @param variables Character vector of variable names to include in plot. Variables can be character, factor, or numeric
#' @param interaction Logical if plot interaction effect of first two variables
#' @param colorID Character string of variable to color data points
#' @param colors Optional character vector of colors for colorID. If not provided, default ggplot is used
#' @param ylab Character string of y-axis lab
#' @param width Numeric figure width in inches
#' @param height Numeric figure height in inches
#' @param outdir Character string of output directory
#' @param cores Numeric cores to run in parallel
#'
#' @return Save individual pdf plots for each gene to outdir
#' @export
#'
#' @examples
#' subset.genes <- c("ENSG00000250479","ENSG00000250510","ENSG00000255823")
#' fdr <- kmFit(dat = dat.voom.example,
#'       patientID = "donorID",
#'       kin = kin.example,
#'       subset.genes = subset.genes,
#'       model = "~ virus + (1|donorID)")
#' gene_plots(dat = dat.voom.example, fdr = fdr, subset.genes = subset.genes,
#'      variables = "virus", colorID = "virus")

gene_plots <- function(
  dat=NULL, counts=NULL, meta=NULL, genes=NULL,
  fdr, model.name=NULL,
  libraryID="libID", geneID="geneName",
  subset.genes=NULL,
  variables, interaction=NULL,
  colorID=NULL, colors=NULL,
  ylab="Normalized log2 expression",
  width=5, height=5, outdir="figs/",
  cores=2){

  #Set seed
  set.seed(4389)
  model <- gene <- variable <- FDR <- i <- E <- NULL
  ########## Check data ##########
  if(!is.null(dat) & !is.null(counts)){ stop("Please provide only one of dat or counts.")}
  if(!is.null(counts) & is.null(meta)){ stop("When using counts, meta must also be provided.")}

  ########## Load data ##########

  if(!is.null(dat) & !is.null(counts)){ stop("Please provide only one of dat or counts.")}
  if(!is.null(counts) & is.null(meta)){ stop("When using counts, meta must also be provided.")}

  if(!is.null(dat)){
    if(class(dat) == "DGEList"){
      dat.counts <- as.data.frame(dat$counts) %>%
        tibble::rownames_to_column(geneID)
      dat.meta <- as.data.frame(dat$samples)
      dat.genes <- as.data.frame(dat$genes)
      } else if(class(dat) == "EList"){
        dat.counts <- as.data.frame(dat$E) %>%
          tibble::rownames_to_column(geneID)
        dat.meta <- as.data.frame(dat$targets)
        dat.genes <- as.data.frame(dat$genes)
      } else {
        stop("dat must be a DGEList or EList object.")
        }
    }
  if(!is.null(counts)){
    #Move rownames, if exist
    if(is.numeric(as.matrix(counts))){
        dat.counts <- as.data.frame(counts) %>%
          tibble::rownames_to_column(geneID)
        } else {
          dat.counts <- as.data.frame(counts)
          colnames(dat.counts)[1] <- geneID
        }}
  if(!is.null(meta)){
    dat.meta <- as.data.frame(meta)
  }
  if(!is.null(genes)){
    dat.genes <- as.data.frame(genes)
  }

  if(!is.null(fdr)){
    if(!is.null(model.name)){
    dat.fdr <- dplyr::select(fdr, model, gene, variable, FDR) %>%
      dplyr::filter(model==model.name)
    } else {
      dat.fdr <- dplyr::select(fdr, gene, variable, FDR)
    }
    } else { dat.fdr <- fdr}

  ########## Combine data ##########
  plot.dat <- dat.counts %>%
    tidyr::pivot_longer(-geneID, names_to = libraryID,
                        values_to = "E") %>%
    dplyr::left_join(dat.meta, by=libraryID) %>%
    dplyr::left_join(dat.genes, by=geneID)

  ########## Plots ##########
  #List all genes/modules
  if(!is.null(subset.genes)){
    to_plot <- subset.genes
  } else{
    #list all genes
    to_plot <- sort(unique(unlist(plot.dat[,geneID])))
  }

  # Setup parallel computing
  doParallel::registerDoParallel(cores=cores)

  ##########  Loop through genes ##########
  foreach::foreach(i = 1:length(to_plot), .verbose = TRUE) %dopar% {

    #Subset data to gene/module of interest
    plot.dat.sub <- plot.dat %>%
      dplyr::filter(get(geneID) == to_plot[i])

    ########## Loop through variables ##########
    plot_list = list()

    #fdr table
    if(!is.null(fdr)){
      dat.fdr.sub <- dplyr::filter(dat.fdr, gene==to_plot[i])
      fdr.plot <- ggpubr::ggtexttable(dat.fdr.sub, rows = NULL)
    } else { fdr.plot = NULL}

    for(j in 1:length(variables)){

      #Type = factor or character variables
      if(is.factor(plot.dat.sub[[variables[j]]]) |
         is.character(plot.dat.sub[[variables[j]]])) {

        #Force first level if exist
        plot.dat.sub.fct <- plot.dat.sub %>%
          dplyr::mutate_at(dplyr::vars(variables[j]), ~forcats::fct_relevel(as.factor(.),
                                                "none", "media","control",
                                                after = 0))

        plot1 <- plot.dat.sub.fct %>%
          ggplot2::ggplot(ggplot2::aes_string(x=variables[j], y="E")) +
          ggplot2::geom_boxplot(outlier.shape = NA) +
          ggplot2::geom_jitter(ggplot2::aes_string(color=colorID), height=0, width=0.2) +
          ggplot2::theme_classic() +
          ggplot2::labs(y=ylab) +
          ggplot2::theme(legend.position = "none")

        if(is.null(colors)){
          plot_list[[variables[j]]] <- plot1
        } else{
          plot1 <- plot1 + ggplot2::scale_color_manual(values=colors)
          plot_list[[variables[j]]] <- plot1
        }

      } else
        #Type = numeric
        if(is.numeric(plot.dat.sub[[variables[j]]])){
          plot1 <- plot.dat.sub %>%
            ggplot2::ggplot(ggplot2::aes_string(x=variables[j], y="E")) +
            ggplot2::geom_point(ggplot2::aes_string(color=colorID)) +
            ggplot2::theme_classic() +
            ggplot2::labs(y=ylab) +
            ggplot2::theme(legend.position = "none") +
            ggplot2::geom_smooth(method='lm', formula= y~x, color="black",
                        se=FALSE)

          if(is.null(colors)){
            plot_list[[variables[j]]] <- plot1
          } else{
            plot1 <- plot1 + ggplot2::scale_color_manual(values=colors)
            plot_list[[variables[j]]] <- plot1
          }
        } else{
          stop("Variables of interest must be numeric, character, or factor.")
        }

    }

    #Interaction plot
    if(!is.null(interaction)){
      plot3 <- plot.dat.sub %>%
        ggplot2::ggplot(ggplot2::aes(x=paste(get(interaction[1]),get(interaction[1]), sep=":"),
                   y=E)) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_jitter(ggplot2::aes_string(color=colorID), height=0, width=0.2) +
        ggplot2::theme_classic() +
        ggplot2::labs(y=ylab, x="", color=colorID) +
        ggplot2::theme(legend.position = "right")

      ##### Color points #####
      if(!is.null(colors)){
        plot3 <- plot3 + ggplot2::scale_color_manual(values=colors)
      }
      } else{
      plot3 <- NULL
    }

    #### Combine plots ####
    #Main plot title
    if (!is.null(dat.genes)){
      hgnc <- dat.genes %>%
        dplyr::filter(get(geneID) == to_plot[i]) %>%
        dplyr::select(dplyr::contains("hgnc"), dplyr::contains("HGNC")) %>%
        unlist()

      title <- paste(to_plot[i], unique(hgnc), sep=" ", collapse=" ")
      plot_row0 <- cowplot::ggdraw() +
        cowplot::draw_label(title, fontface='bold', x=0, hjust=0, vjust=4) +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 7))
    } else{
      titel <- to_plot[i]
      plot_row0 <- cowplot::ggdraw() +
        cowplot::draw_label(title,
                   fontface='bold', x=0, hjust=0, vjust=4) +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 7))
    }

    plot_row2 <- cowplot::plot_grid(plotlist = plot_list,
                           align="hv", nrow=1)

    if(!is.null(plot3) & !is.null(fdr.plot)){
      plot_final <- cowplot::plot_grid(plot_row0, fdr.plot, plot_row2, plot3,
                              nrow=4,
                              rel_heights = c(0.4,0.4,1,1))
    } else if(!is.null(plot3)){
      plot_final <- cowplot::plot_grid(plot_row0, plot_row2, plot3,
                                       nrow=4,
                                       rel_heights = c(0.4,1,1))
    } else if(!is.null(fdr.plot)){
      plot_final <- cowplot::plot_grid(plot_row0, fdr.plot, plot_row2,
                                       nrow=4,
                                       rel_heights = c(0.4,0.4,1))
      } else {
      plot_final <- cowplot::plot_grid(plot_row0, fdr.plot, plot_row2,
                              nrow=3,
                              rel_heights = c(0.4,0.4,1))
    }

    #### Save to disk
    dir.create(path=outdir, showWarnings = FALSE)

    if(!grepl("/$", outdir)) { outdir <- paste(outdir, "/", sep="")}
    filename <- paste(outdir, gsub(" ","_",title),".pdf", sep="")
    ggplot2::ggsave(filename, plot_final, width=width, height=height)
  }

  print("All plots complete.")
}
