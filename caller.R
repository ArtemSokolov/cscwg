## Makes cell type calls based marker expression and marker -> cell type mapping
##
## by Artem Sokolov

library( tidyverse )

main <- function()
{
    ## Name of the column that contains cluster assignments
    iCl <- "Cluster"

    ## Load clustering results
    Xraw <- read_csv( "segResultsRF-pam3.csv.gz", col_types=cols() ) %>% mutate_at( iCl, as.factor )

    ## Load the marker -> cell type mapping
    ## Determine the overlap between the two marker spaces
    M <- read_csv( "markers-v1.0.csv", col_types=cols() ) %>% filter(!is.na(Class)) %>%
        filter( Marker %in% colnames(Xraw) ) %>% select( -Note ) %>% group_by(Class) %>%
        mutate( Norm = n() ) %>% ungroup

    ## Display the name of markers used for making calls
    if( nrow(M) == 0 )
        stop( "No mapping available for supplied marker names" )
    cat( "Using the following", nrow(M), "markers to infer cell type:", str_flatten(M$Marker, ", "), "\n" )

    ## Isolate the markers of interest
    X <- Xraw %>% select( iCl, one_of(M$Marker) ) %>% gather( Marker, Value, -Cluster ) %>%
        nest( Value, .key=In ) %>% mutate_at( "In", map, pull, "Value" )

    ## Compare the expression of each marker between in-cluster and out-of-cluster
    R <- X %>% group_by( Marker ) %>% mutate( Out = map( 1:n(), ~flatten_dbl(In[-.x]) ) ) %>% ungroup %>%
        mutate( wcx = map2_dbl(In, Out, ~wilcox.test(.x,.y,"g")$p.value), Ave = map_dbl(In,mean) ) %>%
        mutate( Vote = 1 - wcx ) %>% select( -In, -Out, -wcx )

    ## Plot the votes
    v <- c(8, 10, 4, 2, 5, 6, 9, 11, 3, 12, 1, 7)
    RP <- R %>% mutate_at( "Cluster", ~str_c("Cluster", as.character(.x)) ) %>% select( -Ave ) %>%
        spread( Marker, Vote ) %>% as.data.frame() %>% column_to_rownames("Cluster") %>% as.matrix %>% round
    pheatmap::pheatmap( RP[,v], legend=FALSE, cluster_cols=FALSE, cluster_rows=FALSE,
                       filename="votes.png", width=7, height=3.5 )

    ## Map markers to cell types and tally the votes
    inner_join( R, M, by="Marker" ) %>% group_by(Class) %>% mutate( Vote = Vote / Norm ) %>%
        group_by( Cluster, Class ) %>% summarize( Votes = sum(Vote) ) %>% mutate( Votes = Votes / sum(Votes) ) %>%
        ungroup %>% spread( Class, Votes )
}

