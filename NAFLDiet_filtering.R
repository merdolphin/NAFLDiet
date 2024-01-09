#### Prevalence Filtering
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ex1b),
               MARGIN = ifelse(taxa_are_rows(ex1b), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ex1b),
                    tax_table(ex1b))

plyr::ddply(prevdf, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Subset to the remaining phyla
prevdf1 = subset(prevdf, phylum %in% get_taxa_unique(ex1b, "phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ex1b),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  
  geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +  
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none") +
  clean_background

table(tax_table(ent.trim)[,"phylum"], exclude=NULL)