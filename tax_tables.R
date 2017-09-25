remain = function(x){
	1-sum(x)
}

tax_bar = function(phyl_rel, rank, cutoff, colours=NULL,){
	## phyl_rel:relative abundance phyloseq object
	## rank:	taxonomic rank to glom to
	## cutoff:	proportion below which things get called 'other'
	## colours: a character vector with the right number of hex values. If null,
	## uses a 21-value vector selected by k-means clustering.
	## assumes the sample ID column is SampleID

	require(phyloseq)
	require(dplyr)

	# Glom to the correct taxonomic rank
	phyl_glommed = tax_glom(phyl_rel, taxrank = rank)

	# Set all counts < 2% to zero
	otu_table(phyl_glommed)[otu_table(phyl_glommed) < cutoff] = 0

	# Filter out all taxa that are zero (<cutoff) everywhere, melt, and sort
	phyl_glommed %>%
		filter_taxa(function(x) sum(x) > 0, prune = TRUE) %>%
		psmelt() %>%
		arrange(Family) -> abunds

## Remove the unneeded columns (This you're going to want to play with
## The idea is to separate out the taxonomy columns (everything after (and
## including) Rank1) so that you can set them all to 'Other'

## I just didn't want these columns any more
# lose = c('OTU','X.SampleID','BarcodeSequence',
#	'LinkerPrimerSequence','Description')

## Keep all the ones I still want
#keep = colnames(abunds)[!(colnames(abunds) %in% lose)]
#abunds = abunds[,keep]

## Divide the columns into sample metadata (up to Rank1)
#metacols = keep[c(1,3:11)]

## and taxonomy metadata (Rank1, Phylum, Class, etc)
#taxcols = keep[12:16]

	# Make an 'Other' row for each sample
	abunds %>% 
		group_by(SampleID) %>%
		summarize(Abundance=remain(Abundance)) -> others
	others[taxcols] = 'Other'		HANDLE TAXCOLS HERE

	# Combine the 'Other' data frame with the original
	abunds = rbind(as.data.frame(others),
				as.data.frame(abunds))
	abunds[,rank] = factor(abunds[,rank], levels = unique(abunds[,rank]))


	# Plot

	# Pick colours
	if (is.null(colours)){
		colours = c("#4f8579","#783fcc","#69d24d","#cb4bbd","#c6dc46","#542871",
					"#78d792","#cc4472","#83d7d0","#d44d33","#676fcd","#ceb854",
					"#403d57","#b97839","#84a4cb","#588038","#c68ac4","#48472a",
					"#c9c39c","#6e2b34","#c78889")
	}

	# Plot individual samples
	abunds = arrange(abunds, PatientID, Family)
indiv = ggplot(abunds, aes(x = Sample, y = Abundance, fill = Family)) +
	geom_bar(stat = "identity") +
	theme(axis.title.x = element_blank(), 
			axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) + 
	scale_fill_manual(values = colours) +
	facet_grid(~GroupID, scales = 'free', space = 'free')+
	ylab("Relative Abundance (Families > 2%) \n")
indiv

# Plot averages
avg = ggplot(abund_fam_means, aes(x = GroupID, y = Average, fill = Family)) +
	geom_bar(stat = "identity") +
	theme(axis.title.x = element_blank(), 
			axis.text.x = element_text(size = 10)) +
	scale_fill_manual(values = colours) +
	ylab("Relative Abundance (Families > 2%) \n")
avg

fecal_rel %>% 
	tax_glom(taxrank = "Family") %>%
	psmelt() %>%
	filter(Abundance > 0.02) %>%
	arrange(Family)-> abunds

