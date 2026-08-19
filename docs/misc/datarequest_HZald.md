
Written for H. Zald
https://github.com/vvandermeersch/climategrowthshifts/issues/56

Estimating neighborhood effects on tree growth to test the role of shared evolutionary history

Requested by:
EM Wolkovich
D Deng
V Van der Meersch

Effects of the neighborhood on tree growth are well established across most forests, with a higher basal area of surrounding trees correlated with lower growth, and suggesting competition contributes to determining tree growth rates. Research often also finds that competitive effects are stronger when trees are surrounded more by conspecific trees, which fits with the idea that intra-specific competition is generally stronger than inter-specific competition. Given the shared evolutionary history of trees, this effect likely extends to congeneric species and beyond but is harder to test. 

We have been developing a model that estimates the effect of tree neighborhood on its growth by decomposing the effect of the neighborhood into a component that is explained by shared evolutionary history (measured as a phylogenetic correlation matrix) and the remaining componen that is unexplained by evolutionary history. The model requires growth over time data (we have been using tree rings, but repeated DBH also works) and the neighborhood data for each measured tree defined by the basal area (or similar) of each species (in some set radius, e.g., 10-15 m). Thus far the model estimates these effects as 'static' -- we assume the neighborhood is unchanged over time, but the model could be updated to estimated the change in growth in response to a change in the neighborhood environment. 

The current model formulation is here: https://github.com/vvandermeersch/climategrowthshifts/wiki/Competition-models#model-specification-of-delta-model 

Our request is for tree growth data and surrounding neighborhood data. For each measured tree we need both its growth data, which species it is, and its neighborhood including basal area (or similar) given by species. If the neighborhood changes over time, we could work to take that into account (though to start we would likely model each time-point separately). 