This repository contains code for D. Vial, V. Subramanian, "On the role of clustering in Personalized PageRank estimation," which considers the relationship between Personalized PageRank (PPR) estimation complexity and clustering in the underlying graph. A [preprint](https://arxiv.org/abs/1706.01091) is available on arXiv.

As a disclaimer, the code will be confusing without the preprint. However, comments reference relevant preprint sections, so the preprint need not be read in its entirety. If questions persist, my email can be found in the preprint.

## Overview
First, a brief explanation of the sub-directories:
- *algo* contains the underlying algorithms used in our experiments
- *data* contains code that pre-processes real datasets and generates synthetic datasets; once datasets are processed/generated, they are also saved in this directory
- *results* will contain experimental output

The remaining files in this directory are MATLAB scripts used for experiments in the preprint. Generally, these scripts proceed as follows: first, a graph (or set of graphs) is loaded from the *data* directory; next, experiments are run using code from the *algo* directory; finally, the output is saved as a *.mat* file in the *results* directory.

See below for further notes on datasets and experiments.

## Notes on datasets
### Real graphs
Our experiments use graphs from the [Stanford Network Analysis Project](http://snap.stanford.edu) (SNAP); all are available [here](http://snap.stanford.edu/data/index.html). SNAP publishes these as *.txt* edge-list files (i.e. each line contains two integers, corresponding to two nodes that share an edge). Our code requires that these are first processed as a particular *.mat* file (details below). Among other things, the pre-processing extracts the largest strongly-connected component from the graph, for which we use [David Gleich](https://www.cs.purdue.edu/homes/dgleich/)'s Graphs Algorithms in MATLAB Code (gaimc). The pre-processing described below assumes the *data* directory contains the *gaimc* directory, available [here](https://www.mathworks.com/matlabcentral/fileexchange/24134-gaimc-graph-algorithms-in-matlab-code) on MATLAB's File Exchange.

To explain the pre-processing, we use SNAP's *Slashdot-0902* graph (the smallest/most manageable graph considered in the paper). The basic workflow is as follows:
1. Download *Slashdot0902.txt* from SNAP (bottom of [this](http://snap.stanford.edu/data/soc-Slashdot0902.html) page), add to *data* directory
2. From *data* directory, run `python relabelSnap.py Slashdot0902` (this relabels nodes with consecutive integer labels in *data/Slashdot0902-r.txt*)
3.  In MATLAB, and again from *data* directory, run `processGraph('Slashdot0902')` (this yields *data/Slashdot0902.mat*, which is used in our experiments)

More specifically, *data/Slashdot0902.mat* contains a MATLAB struct *G* composed of the following:
- *G.n* and *G.m* are the number of nodes and edges in the graph
- *G.Dout* and *G.Din* are vectors containing out- and in-degrees
- *G.Nout* and *G.Nin* are cell arrays containing out- and in-neighbors

This custom format is useful because the only graph accesses required by our algorithms are degree and neighbor lookups.

### Synthetic graphs
The paper also contains empirical results for [Erdos-Renyi graphs](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model) and the [stochastic block model](https://en.wikipedia.org/wiki/Stochastic_block_model). The *data* directory contains code to generate these graphs and save them in the appropriate format. In particular, running `directConnectEr(n,p)` from the *data* directory will generate the graph format described above, for a directed/strongly-connected E-R  graph of *n* nodes with edge probability *p*. Running `directConnectSbm(n,c,p,q)` will do the same, but for an SBM with *n* nodes, *c* communities, intra-community edge probability *p*, and inter-community edge probability *q*. See *data/directConnectEr.m*, *data/directConnectSbm.m*, and preprint Appendix H for more details.

## Notes on experiments

### Single pair algorithms
Preprint Section 4 contains algorithms to estimate PPR for a single source/target node pair. These algorithms can be tested in two stages (assuming graph pre-processing is complete):
1. Run *highPrecisionEstimate.m* with desired parameters. This estimates PPR values to high precision (most graphs on SNAP are too large to actually compute PPR) for accuracy evaluation. Results are saved in *results/{graph name}_hpe.mat*.
2. Run *singlePairComparison.m* with desired parameters. This runs the algorithms from preprint Section 4, using the results of Step 1 to evaluate accuracy. Results are saved in *results/{graph name}_spc.mat*.

For example, assuming *Slashdot0902* has been pre-processed as above, running
```
highPrecisionEstimate({'Slashdot0902'},50,0.2,1,10);
singlePairComparison({'Slashdot0902'},0.2,2e-3,12.2e-3,4.2e-3,7,17);
```
from the main directory will test the Section 4 algorithms for 50 source/target pairs. Here the algorithmic parameters are those used in our experiments (see preprint Table 2). Runtime/accuracy results can be viewed by typing
```
load results/Slashdot0902_spc.mat;
disp(mean(error(:,:,1),2)');
disp(sum(mean(times(:,:,:,1),3)));
```
which should yield something like
```
0.1138  0.1127
0.0086  0.0097
```
indicating that both algorithms had (on average) ~11% relative error and ~9ms runtime (similar to preprint Table 2). See *highPrecisionEstimate.m*, *singlePairComparison.m* code comments for more details.

### Multiple pair algorithms (scalar viewpoint)
Preprint Section 5.1 considers PPR estimation for multiple source/target pairs (with scalar accuracy guarantees). These can be tested via *multiPairComparison.m*. For example, assuming *Slashdot0902* has been pre-processed as above, running
```
multiPairComparison({'Slashdot0902'},0.2,2e-3,12.2e-3,4.2e-3,7,17,5,50,'uni',1);
multiPairComparison({'Slashdot0902'},0.2,2e-3,12.2e-3,4.2e-3,7,17,5,50,'clust',1);
```
will test the algorithms for 5 trials of 50 source/target pairs each, again using preprint Table 2 parameters. Loading/viewing the results via
```
load results/Slashdot0902_uni_1_mpc.mat;
disp(mean(sum(times(:,2,:),1)./sum(times(:,1,:),1),3))
load results/Slashdot0902_clust_1_mpc.mat;
disp(mean(sum(times(:,2,:),1)./sum(times(:,1,:),1),3))
```
should yield something like
```
1.3041
3.4730
```
indicating our scheme is ~1.3x faster than the baseline when pairs are sampled uniformly, and ~3.5x faster when sampled clustered (similar to preprint Figure 7). See *multiPairComparison.m* code comments for further details.

### Multiple pair algorithms (matrix viewpoint)
Preprint Section 5.2 is similar to Section 5.1 but adopts a matrix approximation viewpoint. These algorithms can be tested via *matrixApprox.m*. For example, running
```
matrixApprox({'Slashdot0902'},0.2,2e-3,12.2e-3,4.2e-3,7,17,5,5);
load results/Slashdot0902_ma.mat;
disp(mean(times(:,1,:)./times(:,2:3,:),3));
```
should yield something like
```
    2.2562    2.5881
    1.9076    2.2416
```
indicating our algorithms are ~2-2.5x faster than the baseline (similar to preprint Figure 9). See *matrixApprox.m* code comments for further details.

### Distributed setting
Preprint Section 7 considers a distributed estimation setting; the corresponding script is *distSetting.m*.  For example, running
```
distSetting({'Slashdot0902'},0.2,7,17,2e-3,12.2e-3,5,50,5,5);
load results/Slashdot0902_ds.mat;
disp(mean(sum(times(1,:,:),2)./sum(times(2:3,:,:),2),3));
```
should yield something like
```
2.4649
2.4296
```
indicating the oracle and heuristic methods in Section 7 are ~2-2.5x faster than the baseline method (similar to preprint Figure 10). See *distSetting.m* code comments for more details.

### Final note
The tests above use much smaller source/target sets and far fewer trials than the experiments in the preprint (our goal was simply to validate the code with reasonably low runtime). Consequently, results may differ from those shown here, and from those in the preprint. However, if these tests are scaled up to the level of the preprint, results should be similar to those in the preprint (at the cost of additional runtime), as the preprint shows our empirical findings are statistically significant.
