#' @noRd
filter_arti <- function(x) {
    index <- startsWith(names(x), "node")
    x[rowSums(x[index] == "ARTI") == 0 & 0.5 <= x["omega"], ]
}

#' @rdname doc_pcst
#' @export
run_pcst <- function(lst_net, beta_max, omega_min, omega_max) {
    # extract dataframes
    df_nodes <- lst_net$nodes[, c(1, 2)]
    df_edges <- lst_net$edges[, c(1, 2, 3)]

    # generate template network
    df_nodes_tmp <- data.frame(node = "ARTI", prize = 1)
    df_edges_tmp <- data.frame(
        node1 = "ARTI", node2 = df_nodes$node, cost = NA)

    # bind to old network
    df_nodes_tmp <- rbind(df_nodes_tmp, df_nodes)
    df_edges_tmp <- rbind(df_edges, df_edges_tmp)

    # find index-0 node indices
    edge_inds <- apply(df_edges_tmp[, c(1, 2)], 2, function(x) {
        as.integer(match(toupper(x), toupper(df_nodes_tmp$node)) - 1)
    })

    # reticulate import (could fail)
    reticulate::use_condaenv(condaenv = "PCSF", required = TRUE)  #use a specific conda env.
    pcst_fast <- reticulate::import("pcst_fast")

    # for each beta
    lst_beta <- list()
    for (beta in seq_len(beta_max)) {
        # copy of template network
        df_nodes_spc <- df_nodes_tmp
        df_edges_spc <- df_edges_tmp
        index <- is.na(df_edges_tmp$cost)

        # beta as prize multiplier
        df_nodes_spc$prize <- (df_nodes_spc$prize * beta)

        lst_omega <- list()
        for (omega in seq(omega_min, omega_max, 0.1)) {
            # omega as ARTI edge cost
            df_edges_spc$cost[index] <- omega

            # call to Python pcst_fast function
            lst <- pcst_fast$pcst_fast(
                edge_inds,
                df_nodes_spc$prize,
                df_edges_spc$cost,
                as.integer(0),
                as.integer(1),
                "strong",
                as.integer(0)
            )

            lst$beta <- beta
            lst$omega <- round(omega, 1)

            lst_omega <- append(lst_omega, list(lst))
        }
        lst_beta <- append(lst_beta, lst_omega)
    }

    # collect the resulting lists
    df_nodes_all <- NULL
    df_edges_all <- NULL
    for (lst in lst_beta) {
        # grab run parameters
        beta <- lst$beta
        omega <- lst$omega

        # find node and edge names
        node_names <- df_nodes_tmp[lst[[1]] + 1, "node", drop=FALSE]
        edge_names <- df_edges_tmp[lst[[2]] + 1, c("node1", "node2")]

        # if there is data to add, then add it
        if (nrow(node_names) != 0) {
            df_nodes_new <- data.frame(beta, omega, node_names)
            df_nodes_all <- rbind(df_nodes_all, df_nodes_new)
        }
        if (nrow(edge_names) != 0) {
            df_edges_new <- data.frame(beta, omega, edge_names)
            df_edges_all <- rbind(df_edges_all, df_edges_new)
        }
    }

    # return out
    list(nodes = df_nodes_all, edges = df_edges_all)
}

#' @rdname doc_pcst
#' @export
summarize_pcst <- function(lst_pcst) {
    # extract dataframe,
    # remove artificial nodes
    df_edge <- filter_arti(lst_pcst[["edges"]])

    # count occurances, filter to less than 2000
    tab <- table(df_edge$beta, df_edge$omega)
    tab_which <- which(tab < 2000, arr.ind = TRUE, useNames = FALSE)

    # summary dataframe
    data.frame(
        beta = rownames(tab)[tab_which[, 1]],
        omega = colnames(tab)[tab_which[, 2]],
        n_edges = tab[tab_which]
    )
}

#' @rdname doc_pcst
#' @export
ks_test_pcst <- function(lst_pcst) {
    # extract dataframe,
    # remove artificial nodes
    df_edge <- filter_arti(lst_pcst[["edges"]])

    vec_param <- paste(df_edge$beta, df_edge$omega)
    vec_edges <- paste(df_edge$node1, df_edge$node2)
    tab_edges <- table(vec_edges)

    vec_counts <- tab_edges[cmatch(vec_edges, names(tab_edges))]
    lst_counts <- tapply(vec_counts, vec_param, c)

    # parallel loop for Kolmogorov-Smirnov test
    i <- NULL
    vec_pval <- foreach::`%dopar%`(
        foreach::foreach(i = seq_len(length(lst_counts)), .combine = c), {
        suppressWarnings(stats::ks.test(
            lst_counts[[i]], unlist(lst_counts[-i]),
            alternative = "less"
        )$p.value)
    })

    # order by low to hight p-values
    index <- order(as.numeric(vec_pval))
    mat_param <- do.call(rbind, strsplit(names(lst_counts)[index], " "))

    # test dataframe
    data.frame(
        beta = as.numeric(mat_param[, 1]),
        omega = as.numeric(mat_param[, 2]),
        pval = vec_pval[index],
        row.names = NULL
    )
}

#' @rdname doc_pcst
#' @export
extract_network <- function(lst_net, lst_pcst, mat_pem, beta, omega) {
    df_param <- data.frame(beta, omega)
    df_select <- filter_arti(merge(df_param, lst_pcst$edges))[, -c(1, 2)]

    node_names <- lst_net$nodes$node
    node_prizes <- lst_net$nodes$prize

    df_net <- data.frame(
        do.call(rbind, strsplit(df_select$node1, "__")),
        do.call(rbind, strsplit(df_select$node2, "__")),
        node_prizes[cmatch(df_select$node1, node_names)],
        node_prizes[cmatch(df_select$node2, node_names)],
        merge(df_select, lst_net$edges)$cost
    )

    names(df_net) <- c(
        "node1", "node1_type", "node2", "node2_type",
        "node1_prize", "node2_prize", "cost"
    )

    # new columns
    rownames(mat_pem) <- totitle(rownames(mat_pem))
    df_net$is_ct_edge <- (df_net$node1_type != df_net$node2_type)
    df_net$node1_pem <- apply(df_net, 1, function(row) {
        mat_pem[row["node1"], row["node1_type"]]
    })
    df_net$node2_pem <- apply(df_net, 1, function(row) {
        mat_pem[row["node2"], row["node2_type"]]
    })

    # reorder columns
    index <- c(
        "node1", "node2", "node1_type", "node2_type",
        "node1_prize", "node2_prize", "node1_pem", "node2_pem",
        "is_ct_edge", "cost"
    )

    # return out
    df_net[, index]
}

#' @rdname doc_pcst
#' @export
extract_pathways <- function(df_net, cell_type_a, depth) {
    # copy of dataframe
    df_cpy <- df_net

    # make names unique
    df_cpy$node1 <- paste0(
        df_cpy$node1, ifelse(df_cpy$node1_type == cell_type_a, "", "_"))
    df_cpy$node2 <- paste0(
        df_cpy$node2, ifelse(df_cpy$node2_type == cell_type_a, "", "_"))

    # subset to crosstalk edges
    df_lig <- df_cpy[df_cpy$is_ct_edge, ]

    if (nrow(df_lig) == 0) {
        return()
    }

    # loop through each
    df_path <- do.call(rbind, apply(df_lig, 1, function(row) {
        nodes_all <- unlist(df_cpy[, c(1, 2)])
        nodes_mat <- matrix(nodes_all, ncol = 2)

        nodes_sel <- c()
        nodes_new <- unlist(row[c(1, 2)])

        for (d in seq_len(depth)) {
            nodes_sel <- unique(c(nodes_sel, nodes_new))
            index <- rowSums(matrix(nodes_all %in% nodes_sel, ncol = 2)) != 0
            nodes_new <- as.vector(nodes_mat[index, ])
        }

        cbind(pathway = paste(row[c(1, 2)], collapse = "--"), df_cpy[index, ])
    }))

    keys <- unique(df_path$pathway)
    lst_path <- lapply(keys, function(x) {
        df_path[df_path$pathway == x, -1]
    })
    names(lst_path) <- keys

    # return out
    lst_path
}
