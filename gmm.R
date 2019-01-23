## Callers based on Gaussian Mixture Models and Naive Bayes
##
## by Artem Sokolov

library( tidyverse )

## Clamps qq^th quantile to -1 and (1-qq)^th quantile to 1
quantnorm <- function( v, qq = 0.001 )
{
    stopifnot( qq < .5 )
    lo <- quantile( v, qq )
    hi <- quantile( v, 1 - qq )
    2 * (v - lo) / (hi - lo) - 1
}

## Fits a two-component mixture model to the provided 1-D vector
mixEM <- function( v, tag )
{
    cat( "Fitting a mixture for", tag, "\n" )
    mix <- mixtools::normalmixEM( v, mu=c(-0.8,0.8), maxit=2000 )
    mix[c("lambda","mu","sigma")]
}

## Given a model returned by mixEM(), computes posterior probabilities
##   for column vals in data frame .df
mutate_probs <- function( .df, vals, mix )
{
    ## Ensure the means are ordered
    stopifnot( mix$mu[1] < mix$mu[2] )
    s <- ensym(vals)

    ## Compute posterior probabilities
    .df %>% mutate( Neg = ifelse(!!s < -1, 1, mix$lambda[1]*dnorm(!!s, mix$mu[1], mix$sigma[1])),
                   Pos = ifelse(!!s < -1, 0, mix$lambda[2]*dnorm(!!s, mix$mu[2], mix$sigma[2])) ) %>%
        mutate( Pos = ifelse( !!s > 1, 1, Pos ), Neg = ifelse( !!s > 1, 0, Neg ) ) %>%
        mutate( .tmp.sum = Neg + Pos, Pos = Pos / .tmp.sum, Neg = Neg / .tmp.sum ) %>%
        select( -.tmp.sum )
}

## Plots values from a single cell against background distribution of each marker
plotCell <- function( X, CL )
{
    CLn <- CL %>%
        mutate( Value = ifelse(Value < -1, -1, Value) ) %>%
        mutate( Value = ifelse(Value > 1, 1, Value) )
    ggplot( X, aes(x=Value) ) + theme_bw() + xlim( c(-1.1, 1.1) ) +
        geom_density() + facet_wrap( ~Marker, nrow=4 ) +
        geom_vline( aes(xintercept=Value), data = CLn, color="red" )
}

main <- function()
{
    ## Load the segmented matrix
    Xraw <- read_csv( "segResultsRF.csv.gz", col_types=cols() ) %>%
        mutate( CellID = str_c("Cell", 1:n()) )

    ## Load markers -> cell types mapping
    M <- read_csv( "markers-v1.0.csv", col_types=cols() ) %>%
        filter( !is.na(Class), Marker %in% colnames(Xraw) ) %>%
        select( -Note )

    ## Consider markers for which there's a mapping
    ## Apply a log transform and quantile normalization
    X <- Xraw %>% select( CellID, one_of(M$Marker) ) %>% mutate_at( vars(-CellID), ~log10(.x+1) ) %>%
        mutate_at( vars(-CellID), quantnorm, 0.01 ) %>% gather( Marker, Value, -CellID )

    ## Fit mixture models to the outlier-free regions
    ## MX <- X %>% filter( Value >= -1, Value <= 1 ) %>% nest( -Marker, .key=Values ) %>%
    ##     mutate_at( "Values", map, pull, "Value" ) %>%
    ##     mutate( Fit = map2( Values, Marker, mixEM ) ) %>% select( -Values )
    ## save( MX, file="em-fits.RData" )
    load( "em-fits.RData" )

    ## Retrieve and plot the means of each fit
    MM <- MX %>% mutate( Value = map(Fit, "mu") ) %>% unnest(Value)
    plotCell( X, MM )
    ## + ggsave( "em-fits.pdf", width=8, height=6 )

    ## Compute posterior probabilities for each sample / model pair
    P <- X %>% nest( -Marker, .key=Values ) %>% inner_join( MX, by="Marker" ) %>%
        mutate( PostProb = map2(Values, Fit, ~mutate_probs(.x, "Value", .y)) ) %>%
        select( -Values, -Fit ) %>% unnest()

    ## For any given marker, choose the positive Gaussian probability for the corresponding
    ##   class and the negative Gaussian probability otherwise
    ee <- unique(M$Class) %>% set_names %>% map( ~expr(ifelse(Class == !!.x, Pos, Neg)) )
    Y <- inner_join( P, M, by="Marker" ) %>% mutate( !!!ee ) %>% select( -Neg, -Pos, -Class )

    ## Compute the (unnormalized) posterior probability of for each cell/class combo
    R <- Y %>% group_by( CellID ) %>% summarize_at( c("Tumor","Immune","Stroma"), prod )

    ## Match cell probabilities against their positions and save everything to a file
    Xraw %>% select( CellID, Pos_X = CellPosition_X, Pos_Y = CellPosition_Y ) %>%
        inner_join( R, by="CellID" ) %>% write_csv( "segResultsRF-calls.csv.gz" )
    
    ## Identify the top 1000 cells in each class
    R1k <- R %>% gather( Class, Prob, -CellID ) %>% group_by( Class ) %>%
        top_n( 1000, Prob ) %>% ungroup()

    ## Match them up against cell positions
    F <- Xraw %>% select( CellID, Pos_X = CellPosition_X, Pos_Y = CellPosition_Y ) %>%
        inner_join( R1k, by="CellID" ) %>% select( -Prob )
    write_csv( F, "Top1000-labels.csv" )
}
