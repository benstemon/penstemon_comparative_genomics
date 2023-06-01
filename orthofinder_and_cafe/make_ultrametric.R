library(ape)


intree <- read.tree(file = "SpeciesTree_rooted.tre")

#generate tree calibration:
#from Wolfe et al 2022
#P. davidsonii -> any other non-Dasanthera: 1.419-6.416
treecal <- makeChronosCalib(intree,
                            node = getMRCA(intree, c("single_isoform_davidsonii_FUNCTIONAL-INCLUDED-nucleotide-CMKEVH",
                                                     "single_isoform_petiolatus-nucleotide-CMKEVH")),
                            age.min = 1.419, age.max = 6.416)


#make calibrated tree
calibrated_tree <- chronos(intree, lambda = 1, model = "correlated", calibration = treecal,
                           control = chronos.control(iter.max = 1e6, eval.max = 1e6))
#plot(calibrated_tree)
#is.ultrametric(calibrated_tree)
write.tree(calibrated_tree, file = "calibrated_species_tree.tre")
