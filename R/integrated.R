#' @rdname doc_integrated
#' @export
nonselftalk <- function(mat_type, lrp) {
    # make sure function takes in a Matrix type
    mat_type <- Matrix::Matrix(mat_type)

    lrp_index <- match_lr_pairs(mat_type, lrp)
    errorifnot(0 < nrow(lrp_index), "no ligand-receptor pairs found")

    index <- which(Matrix::rowSums(mat_type != 0) == 0)

    if (length(index) > 0) {
        mat_type[index, ] <- add_noise(mat_type[index, ])
    }
    mat_disc <- discretize_sparse(
        Matrix::t(mat_type), "equalwidth", max(2, ncol(mat_type)^(1/2)))

    mi <- NULL
    for (i in seq_len(nrow(lrp_index))) {
        row <- unlist(lrp_index[i, ])
        x <- mat_disc[, row[1]]
        y <- mat_disc[, row[2]]
        mi <- c(mi, mutinfo_xy(x, y, "mm", TRUE))
    }

    -log10(ifelse(mi < 0, 0, mi))
}

#' @rdname doc_integrated
#' @export
gene_relevance <- function(mat_intra, lrp) {
    # create a new matrix, same size
    mat_nsq <- as.matrix(mat_intra) * 0
    # set the diagonal as the network rowsums
    diag(mat_nsq) <- Matrix::rowSums(mat_intra)

    # negative square root
    mat_nsq <- corpcor::mpower(mat_nsq, -0.5)
    # symmetric matrix
    mat_wnorm <- mat_nsq %*% as.matrix(mat_intra) %*% mat_nsq
    # match genes to LR pair
    lrp_index <- !is.na(cmatch(rownames(mat_intra), unlist(lrp)))

    # compute gene relevance value,
    # random walk with restart
    n <- 50
    alpha <- 0.9
    offset <- (1 - alpha) * lrp_index
    mat_relev <- as.matrix(lrp_index)

    for (i in seq_len(n)) {
        mat_relev <- alpha * (mat_wnorm %*% mat_relev) + offset
    }

    vec_relev <- as.numeric(mat_relev)
    names(vec_relev) <- rownames(mat_intra)

    vec_relev
}

#' @rdname doc_integrated
#' @export
node_prize <- function(mat_pem, cell_type, vec_relev) {
    # match PEM cell type
    vec_pem <- mat_pem[, cell_type == colnames(mat_pem)]

    # match gene relevance names to cell type vector
    index <- match(toupper(names(vec_relev)), toupper(names(vec_pem)))
    vec_pem_match <- vec_pem[index]

    # relevance times cell specific,
    # if negative then zero out
    vec_relev * ifelse(vec_pem_match < 0, 0, vec_pem_match)
}

#' @rdname doc_integrated
#' @export
crosstalk <- function(
    mat_pem, cell_type_a, cell_type_b, vec_nst_a, vec_nst_b, mat_type, lrp) {

    # grab relevant PEM scores, zero out negatives and NaNs
    vec_pem_names <- rownames(mat_pem)
    vec_pem_a <- zero_na_neg(mat_pem[, cell_type_a])
    vec_pem_b <- zero_na_neg(mat_pem[, cell_type_b])

    # zero out bad NST scores, merge the them together
    df_nst <- match_lr_pairs(mat_type, lrp)
    df_nst[, "mi_a"] <- zero_na_neg(vec_nst_a)
    df_nst[, "mi_b"] <- zero_na_neg(vec_nst_b)

    # initialize variables
    df_ct <- data.frame()

    for (i in seq_len(nrow(df_nst))) {
        # find the LR pair
        lig <- df_nst[i, "ligand"]
        rec <- df_nst[i, "receptor"]

        # non-self talk score
        scr_nst <- sum(df_nst[i, c("mi_a", "mi_b")]) / 2

        # expressed score
        scr_expr <- (vec_pem_a[lig] + vec_pem_b[rec]) / 2

        df_ct <- rbind(df_ct, data.frame(
            ligand = lig, receptor = rec, ligand_type = cell_type_a,
            receptor_type = cell_type_b, nst = scr_nst, expr = scr_expr
        ))

        # from type B to A (sometimes skipped)
        if (!identical(lig, rec)) {
            # expressed score (notice the difference!)
            scr_expr <- (vec_pem_b[lig] + vec_pem_a[rec]) / 2

            df_ct <- rbind(df_ct, data.frame(
                ligand = lig, receptor = rec, ligand_type = cell_type_b,
                receptor_type = cell_type_a, nst = scr_nst, expr = scr_expr
            ))
        }
    }

    # compute crosstalk score
    df_ct[, "crosstalk"] <- minmax(df_ct[, "expr"]) * minmax(df_ct[, "nst"])
    rownames(df_ct) <- NULL
    df_ct
}

#' @noRd
extract_lower_nonzero <- function(mat) {
    index <- Matrix::which(
        lower.tri(mat, diag = TRUE) & mat != 0, arr.ind = TRUE
    )
    df <- data.frame(apply(index, 2, function(x) rownames(mat)[x]))
    df[, "val"] <- mat[index]
    df
}

#' @rdname doc_integrated
#' @export
integrate_network <- function(
    vec_nst_a, vec_nst_b, mat_intra_a, mat_intra_b,
    cell_type_a, cell_type_b, mat_pem, mat_type, lrp) {

    # gene relevance
    vec_gr_a <- gene_relevance(mat_intra_a, lrp)
    vec_gr_b <- gene_relevance(mat_intra_b, lrp)

    # node prize
    vec_np_a <- node_prize(mat_pem, cell_type_a, vec_gr_a)
    vec_np_b <- node_prize(mat_pem, cell_type_b, vec_gr_b)

    # crosstalk
    df_edge_ct <- crosstalk(
        mat_pem, cell_type_a, cell_type_b, vec_nst_a, vec_nst_b, mat_type, lrp)

    # extract MI scores from sparse matrix
    df_edge_a <- extract_lower_nonzero(mat_intra_a)
    df_edge_b <- extract_lower_nonzero(mat_intra_b)

    # casefold names, determine node and edge types
    node_names_a <- add_suffix(names(vec_np_a), cell_type_a)
    node_names_b <- add_suffix(names(vec_np_b), cell_type_b)
    edge_names_a <- apply(df_edge_a[, c(1, 2)], 2, add_suffix, cell_type_a)
    edge_names_b <- apply(df_edge_b[, c(1, 2)], 2, add_suffix, cell_type_b)
    edge_names_ct <- cbind(
        add_suffix(df_edge_ct[, "ligand"], df_edge_ct[, "ligand_type"]),
        add_suffix(df_edge_ct[, "receptor"], df_edge_ct[, "receptor_type"])
    )

    # compile full node names
    node_names_full <- unique(c(node_names_a, node_names_b))

    # validate crosstalk edges
    index_ct <- rowSums(apply(edge_names_ct, 2, "%in%", node_names_full)) == 2

    # compile full edges and costs
    edge_full <- rbind(edge_names_a, edge_names_b, edge_names_ct[index_ct, ])
    cost <- list(
        df_edge_a$val, df_edge_b$val, df_edge_ct[index_ct, "crosstalk"])

    # normalize edge costs
    cost_norm <- unlist(lapply(cost, scale))
    cost_norm[is.na(cost_norm)] <- 0
    cost_full <- 1 - minmax(cost_norm)

    # grab relevant PEM scores,
    # zero out negatives and NaNs
    vec_pem_names <- rownames(mat_pem)
    vec_pem_a <- zero_na_neg(mat_pem[, cell_type_a])
    vec_pem_b <- zero_na_neg(mat_pem[, cell_type_b])

    # match PEM
    vec_pem_match_a <- vec_pem_a[cmatch(names(vec_np_a), names(vec_pem_a))]
    vec_pem_match_b <- vec_pem_b[cmatch(names(vec_np_b), names(vec_pem_b))]

    # node dataframe
    df_net_node <- data.frame(
        node = c(node_names_a, node_names_b),
        prize = c(vec_np_a, vec_np_b),
        pem = c(vec_pem_match_a, vec_pem_match_b),
        gene_relevance = c(vec_gr_a, vec_gr_b)
    )

    # edge dataframe
    df_net_edge <- data.frame(
        node1 = edge_full[, 1],
        node2 = edge_full[, 2],
        cost = cost_full
    )

    # return out
    lst_net <- list(nodes = df_net_node, edges = df_net_edge)
    lst_net
}
