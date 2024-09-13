[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4553387.svg)](https://doi.org/10.5281/zenodo.4553387)

# RSS-NET: Regression with Summary Statistics exploiting NEtwork Topology

The present repository contains source codes and documentations of RSS-NET,
a novel Bayesian framework for simultaneous enrichment and prioritization
analysis of complex trait GWAS and gene regulatory networks. 

## Getting started

1. [Install](https://suwonglab.github.io/rss-net/setup.html) the RSS-NET software.

2. Try RSS-NET on [two synthetic datasets](https://suwonglab.github.io/rss-net/wtccc_bcell.html).

3. Try RSS-NET on [a real-world dataset](https://suwonglab.github.io/rss-net/ibd2015_nkcell.html). 

## Citing this work

If you find any part of this repository useful for your work,
please kindly cite the following research article:

> Zhu, X., Duren, Z. & Wong, W.H.
> Modeling regulatory network topology improves genome-wide analyses of complex human traits.
> *Nat Commun* 12, 2851 (2021). <https://doi.org/10.1038/s41467-021-22588-0>

We originally developed RSS-NET to integrate GWAS with gene regulatory networks,
as implemented in [`rss_net.m`](src/rss_net.m).
We recently extended RSS-NET to integrate GWAS with other genomic annotations such as
[sequence-conserved enhancers](https://github.com/SUwonglab/m2h-ele),
and this extension is available as [`rss_gset.m`](src/rss_gset.m).
If you find this extension useful for your work,
please kindly cite the following research article,
in addition to the original [RSS-NET publication](https://doi.org/10.1038/s41467-021-22588-0).

> Zhu, X., Ma, S. & Wong, W.H.
> Genetic effects of sequence-conserved enhancer-like elements on human complex traits.
> *Genome Biol* 25, 1 (2024). <https://doi.org/10.1186/s13059-023-03142-1>

Correspondence should be addressed to X.Z. and W.H.W.

## Support

1. Refer to [RSS-NET wiki](https://SUwonglab.github.io/rss-net/)
for tutorials and documentations.

2. Create a new [GitHub issue](https://github.com/SUwonglab/rss-net/issues)
to report bugs and/or request features.

## Contact

[Xiang Zhu, Ph.D.](https://github.com/xiangzhu)<br>
[Wing Hung Wong Lab](https://statistics.stanford.edu/people/wing-hung-wong)<br>
[Department of Statistics](https://statistics.stanford.edu/)<br>
[Stanford University](https://www.stanford.edu/)


