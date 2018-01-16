remain = function(x){
	1-sum(x)
}

prop_tax_row = function(taxrow){
	## taxrow: a row from a tax_table of a phyloseq object

	ranks = c('Phylum','Class','Order','Family','Genus')

	hasNA = FALSE
	for (i in 1:length(ranks)){

	    if (!is.na(taxrow[,ranks[i]])){
	        tax = taxrow[,ranks[i]]
	    } else {
	        hasNA = TRUE
	        break
	    }
	}

	if (!hasNA){
		return (taxrow)
	} else {
		init = tolower(substring(ranks[i-1],1,1))
		assn = paste(init,tax,sep = '_')
		taxrow[,ranks[i:length(ranks)]] = assn
		return(taxrow)
	}
}

prop_tax_down = function(taxtab){
    # taxtab: a taxonomy table to propogate taxa down in

    ## I don't know why I can't use apply for this, but I can't.
    for (r in 1:nrow(taxtab)){
        taxtab[r,] = prop_tax_row(taxtab[r,])
    }

    return(taxtab)
}

taxa_other_df = function(phyl_rel, rank, cutoff){
	## phyl_rel:relative abundance phyloseq object
	## rank:	taxonomic rank to glom to
	## cutoff:	proportion below which things get called 'other'

	require(phyloseq)
	require(dplyr)

	# Propogate taxonomic assignments down the tree
    tax_table(phyl_rel) = prop_tax_down(tax_table(phyl_rel))

	# Glom to the correct taxonomic rank
	phyl_glommed = tax_glom(phyl_rel, taxrank = rank)

	# Set all counts < 2% to zero
	otu_table(phyl_glommed)[otu_table(phyl_glommed) < cutoff] = 0

	# Filter out all taxa that are zero (<cutoff) everywhere, melt, and sort
	phyl_glommed %>%
		filter_taxa(function(x) sum(x) > 0, prune = TRUE) %>%
		psmelt() %>%
		arrange_(rank) %>%
		data.frame() -> abunds

	# Order the rank by mean abundance

	abunds %>%
		group_by(UQ(sym(rank))) %>%
		summarize(Mean = mean(Abundance)) %>%
		data.frame() -> mean_abunds

	lev_ord = levels(mean_abunds[,rank])
	lev_ord = lev_ord[order(mean_abunds$Mean)]

	abunds[,rank] = factor(abunds[,rank], levels = lev_ord)

	# List all the metadata columns so that they are included in the data frame
	metacols = names(abunds)[4:(match('Rank1',names(abunds))-1)]
	# Make an 'Other' row for each sample
	abunds %>%
		group_by_(.dots = metacols) %>%
		summarize(Abundance=remain(Abundance)) -> others

	# Add in the taxonomic data columns
	taxcols = names(abunds)[match('Rank1',names(abunds)):ncol(abunds)]
	others[taxcols] = 'Other'

	# Combine the 'Other' data frame with the original
	newdf = abunds[,metacols]
	newdf$Abundance = abunds$Abundance
	newdf[,taxcols] = abunds[,taxcols]
	newdf = rbind(as.data.frame(others),
				as.data.frame(newdf))
	newdf[,rank] = factor(newdf[,rank], levels = c('Other',lev_ord))

	return(newdf)

}


taxa_plot = function(taxa_df,rank,colours = NULL,
					 sample = 'X.SampleID', abund = 'Abundance'){
	## taxa_df:	The data frame produced by taxa_other_df()
	## rank: The taxonomic rank to plot by
	## colours:	A character vector with the right number of colours. If you
	## don't provide one it uses my 21-colour vector, which might not be enough
	## left NULL, the bars will be ordered alphabetically by sample ID.
	## sample: the name of the sample ID column.
	## abund: the name of the abundance column

	require(ggplot2)

	# Pick colours
	if (is.null(colours)){
		colours = c('grey69',"#4f8579","#783fcc","#69d24d","#cb4bbd","#c6dc46",
					"#542871","#78d792","#cc4472","#83d7d0","#d44d33","#676fcd",
					"#ceb854","#403d57","#b97839","#84a4cb","#588038","#c68ac4",
					"#48472a","#c9c39c","#6e2b34","#c78889")
	} else {
		colours = c('grey69',colours)
	}

	# Make sure the x axis is categorical
	taxa_df[,sample] = factor(taxa_df[,sample])

	indiv = ggplot(taxa_df, aes_string(x = sample, y = abund, fill = rank)) +
	geom_bar(stat = "identity") +
	theme(axis.title.x = element_blank(),
			axis.text.x = element_text(size = 10,
									   angle = 90,
									   hjust = 1,
									   vjust = 0.5)) +
	scale_fill_manual(values = colours, guide = guide_legend(reverse = TRUE)) +
	ylab(paste("Relative Abundance (",rank,")\n",sep=''))

	return(indiv)
}
