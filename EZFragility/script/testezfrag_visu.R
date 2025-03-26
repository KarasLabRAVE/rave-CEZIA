

data("pt01EcoG")

## Visualize a subject of electrodes
sozIndex <- attr(pt01EcoG, "sozIndex")
display <- c(sozIndex, 77:80)

epoch <- Epoch(pt01EcoG)
visuIEEGData(epoch = epoch[display, ])

data("pt01EcoG")

## sozIndex is the index of the electrodes we assume are in the SOZ
sozIndex <- attr(pt01EcoG, "sozIndex")

## index of the electrodes to display
display <- c(sozIndex, 77:80)

## precomputed fragility object
data("pt01Frag")

stat <- fragStat(pt01Frag, sozIndex)

## plot the fragility heatmap
plotFragHeatmap(frag = pt01Frag, sozIndex = sozIndex)

## plot the fragility quantiles
plotFragQuantile(frag = pt01Frag, sozIndex = sozIndex)

## plot the fragility distribution
plotFragDistribution(frag = pt01Frag, sozIndex = sozIndex)
