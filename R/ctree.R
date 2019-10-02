#' Construct a clone tree.
#'
#' @description 
#'
#' This constructor creates an object of class `'ctree'`, which represents a clone tree for a 
#' patient. The tree must be created from a set of clusters computed for a patient, usually
#' with methods that carry out tumour subclonal deconvolution from bulk sequencing data.
#' To create a tree a list of drivers should be annotated for the input clusters. 
#' 
#' 
#' @param CCF_clusters 
#' @param drivers 
#' @param samples 
#' @param patient 
#' @param M 
#' @param score 
#' @param annotation 


#' @param x A REVOLVER cohort.
#' @param patient The id of this patient in the cohort.
#' @param M The tree adjacency matrix.
#' @param score A score for this model, that we seek to maximize.
#' @param annotation A string to annotate this tree.
#'
#' @return An object of class \code{"rev_phylo"} that represents this tree.
#' @export
#'
#' @import crayon
#' @import tidygraph
#'
#' @examples
#' 
#'
#' 
#'
#' 
#' # We re-build the trees for one patient
#' data(TRACERx_cohort)
#' 
#' # Because the function will return a new cohort, we can just
#' M = ...
#' 
#' # Create a new tree, the score is made up
#' new_tree = revolver_phylogeny(TRACERx_cohort, patient = "CRUK0002", M, score = 17)
#' 
#' # S3 functions for this object
#' print(new_tree)
#' plot(new_tree)
ctree = 
  function(
    CCF_clusters,
    drivers,
    samples,
    patient,
    M,
    score,
    annotation = paste0("Clone tree for patient ", patient))
  {
    # This function will create this output object
    obj <-
      structure(
        list(
          adj_mat = NA,
          # Adjacency matrix for the tree
          tb_adj_mat = NA,
          # Tidygraph object for this tree
          score = NA,
          # Tree score
          patient = NA,
          # Patient ID
          samples = NA,
          # Samples name
          drivers = NA,
          # Driver mapping to clusters
          CCF = NA,
          # CCF (aggregated per cluster)
          transfer = NA,
          # Information Transfer
          annotation = NA,
          # Some custom annotation
          tree_type = "Clone tree"  # Clone tree 
        ),
        class = "ctree",
        call = match.call()
      )
    
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    # The information that we need for each tree are
    # CCF clusters, data and driver events (mapped to clusters)
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    obj$CCF = CCF_clusters
    obj$samples = samples
    obj$drivers = drivers
    
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    # We begin to create a representation of the tree
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    
    # Input can should be an adjacency matrix (which works also for monoclonal tumours)
    if (class(M) != 'matrix')
      stop("Input `M` should be an adjacency matrix, aborting.")
    
    adj_mat =  M
    df_mat = MatrixToDataFrame(M)
    
    # We check for that to be a tree - empty is OK in this case if monoclonal
    if (!is_tree(adj_mat, allow.empty = nrow(obj$CCF) == 1))
      stop("The input adjacency matrix is not a valid tree, aborting.")
    
    # Add GL node, beware of the special case of empty adj_mat (might happen for a monoclonal tumour)
    M_root = ifelse(sum(adj_mat) == 0, rownames(adj_mat), root(adj_mat))
    df_mat = rbind(df_mat,
                   data.frame(
                     from = 'GL',
                     to = M_root,
                     stringsAsFactors = FALSE
                   ))
    
    # Update the adj matrix with GL, which can now go into obj
    obj$adj_mat = DataFrameToMatrix(df_mat)
    
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    # Create a tidygraph object for this tree
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    tb_adj_mat = as_tbl_graph(df_mat) %>%
      activate(nodes) %>%
      rename(cluster = name)
    
    # We can add information specific for this tree to the tidygraph
    tb_adj_mat = tb_adj_mat %>%
      left_join(obj$CCF, by = 'cluster')
    
    # Sample attachment for input data
    attachment = obj$CCF %>%
      select(cluster, obj$samples) %>%
      reshape2::melt(id = 'cluster') %>%
      mutate(value = ifelse(value > 0, 1, 0)) %>%
      group_by(variable) %>%
      filter(sum(value) == 1) %>%
      ungroup() %>%
      filter(value == 1) %>%
      select(-value) %>%
      rename(attachment = variable) %>%
      mutate(attachment = paste(attachment)) %>%
      group_by(cluster) %>%
      summarise(attachment = paste(attachment, collapse = ', '))
    
    tb_adj_mat = tb_adj_mat %>%
      left_join(attachment, by = 'cluster')
    
    # Drivers per node
    tb_adj_mat = tb_adj_mat %>%
      left_join(obj$drivers %>%
                  group_by(cluster) %>%
                  summarise(driver = paste(variantID, collapse = ', ')),
                by = 'cluster')
    
    # Store it in obj
    obj$tb_adj_mat = tb_adj_mat
    
    # Compute the information transfer
    obj$transfer = information_transfer(obj)
    
    # Extras
    obj$score = score
    obj$patient = patient
    
    obj$annotation = annotation
    
    obj$tree_type = ifelse(
      all(obj$CCF %>% select(obj$samples) %in% c(0, 1)),
      "Mutation tree (binary data)",
      "Phylogenetic tree (CCF)"
    )
    
    return(obj)
  }



