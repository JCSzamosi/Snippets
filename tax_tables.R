# Glom to family
fecal_fam = tax_glom(fecal_rel, taxrank = 'Family')

# Set all counts < 2% to zero
otu_table(fecal_fam)[otu_table(fecal_fam) < 0.02] = 0

# Filter out all taxa that are zero (<2%) everywhere, melt, and sort by family
fecal_fam %>%
	filter_taxa(function(x) sum(x) > 0, prune = TRUE) %>%
	psmelt() %>%
	arrange(Family) -> abund_fams

# Remove the unneeded columns (This you're going to want to play with
# The idea is to separate out the taxonomy columns (everything after (and
# including) Rank1) so that you can set them all to 'Other'

# I just didn't want these columns any more
lose = c('OTU','X.SampleID','BarcodeSequence',
	'LinkerPrimerSequence','Description')

# Keep all the ones I still want
keep = colnames(abund_fams)[!(colnames(abund_fams) %in% lose)]
abund_fams = abund_fams[,keep]

# Divide the columns into sample metadata (up to Rank1)
metacols = keep[c(1,3:11)]

# and taxonomy metadata (Rank1, Phylum, Class, etc)
taxcols = keep[12:16]

# Make an 'Other' row for each sample
remain = function(x){
	1-sum(x)
}
abund_fams %>% 
	group_by_(.dots=metacols) %>%
	summarize(Abundance=remain(Abundance)) -> others
others[taxcols] = 'Other'

# Combine the 'Other' data frame with the original
abund_fams = rbind(as.data.frame(others),
			as.data.frame(abund_fams))
abund_fams$Family = factor(abund_fams$Family, levels = unique(abund_fams$Family))


# Plot

# Pick colours
clrs = c("#4f8579","#783fcc","#69d24d","#cb4bbd","#c6dc46","#542871","#78d792",
		"#cc4472","#83d7d0","#d44d33","#676fcd","#ceb854","#403d57","#b97839",
		"#84a4cb","#588038","#c68ac4","#48472a","#c9c39c","#6e2b34","#c78889")

# Plot individual samples
abund_fams = arrange(abund_fams, PatientID, Family)
indiv = ggplot(abund_fams, aes(x = Sample, y = Abundance, fill = Family)) +
	geom_bar(stat = "identity") +
	theme(axis.title.x = element_blank(), 
			axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) + 
	scale_fill_manual(values = clrs) +
	facet_grid(~GroupID, scales = 'free', space = 'free')+
	ylab("Relative Abundance (Families > 2%) \n")
indiv

# Plot averages
avg = ggplot(abund_fam_means, aes(x = GroupID, y = Average, fill = Family)) +
	geom_bar(stat = "identity") +
	theme(axis.title.x = element_blank(), 
			axis.text.x = element_text(size = 10)) +
	scale_fill_manual(values = clrs) +
	ylab("Relative Abundance (Families > 2%) \n")
avg

fecal_rel %>% 
	tax_glom(taxrank = "Family") %>%
	psmelt() %>%
	filter(Abundance > 0.02) %>%
	arrange(Family)-> abund_fams

