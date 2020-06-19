1. Motivating questions
    * When we work with a set of bifurcating trees, but the true history was a network, can our downstream inferences about lineage diversification and trait evolution be misled?
        - Does a burst of hybridization get interpreted as a burst of speciation?
        - Are rates of trait evolution over-estimated?
        - Are reversals inferred for traits that are evolving irreversibly?
        - Do results depend on the genetic architecture of the trait (e.g., single locus vs additive effect of many loci)?
        - Are problems worse when reticulation is recent vs deep in the clade?

2. Obtain a true, network history of relationships among species
    * Largely tree-shaped, but with some reticulation
    * Various possible sources
        - Invent one by hand
        - Use an empirical result
        - Generate by simulation

3. Convert the true network to a collection of bifurcating trees
    * "Geometric" option
        - Break the network in all possible ways
        - Get a set of trees and their probabilities/weights
    * "Inference" option
        - Use `ms` to generate gene trees based on the network
        - Use `seq-gen` to generate DNA sequences for each gene
        - Use `mrbayes` to generate a posterior set of trees from the concatenated gene sequences

4. Other needs
    * If the focal question is about speciation and extinction rates only, don't need anything but the trees
    * If the focal question is about trait evolution, need trait values for the extant species
        - Simulate a discrete or continuous-valued trait on the network's dominant species tree
        - Compute the trait value based on info from `ms/seq-gen`

5. Downstream analyses
    * Use various models of speciation, extinction, trait evolution
    * Fit the relevant model to each tree in the set, and then combine to get overall answer and its uncertainty
