FastAAI v2 updates the concept of FastAAI developed originally by Kenji Gerhardt, Carlos Ruiz-Perez, and Miguel Rodriguez-Rojas. The original FastAAI can be fount at https://github.com/cruizperez/FastAAI.

FastAAI is an ultrafast estimator of whole-genome aggregate amino acid identity, abbreviated as "AAI." AAI is a metric for assessing genome relatedness for moderately to deeply divergent organisms, resolving relationships betwen the phylum and genus levels. However, the calculation of AAI using alignment-based methods such as BLAST or DIAMOND is now too slow to reasonably calculate AAI for comparisons involving modern collections of genomes. FastAAI and FastAAI v2 accurately estimate AAI orders of magnitude more quickly than alignment-based approaches, eliminating this bottleneck. FastAAI achieves this through identifying universally shared, single-copy, protein coding genes (abbreviated as single-copy proteins or SCPs) and comparing genomes over only these shared proteins using the Jaccard similarity index of shared amino acid tetramers. The Jaccard indices of each shared SCP between two genomes are averaged to produce a final point estimate of AAI. This comparsion, while computationally simple, accurately estimates AAI.

The original FastAAI accurately estimated AAI and was fast, but suffered from limitations. Namely:

* The SCP set used by FastAAI v1 was completely fixed and could not be changed
* The chosen set of SCPs prevented the use of FastAAI with non-prokaryotic genomes
* FastAAI v1 used a single mathematical formula for translating Jaccard similarity into AAI; this means that the same model was used for comparing bacteria to bacteria and bacteria to archaea, for example
* FastAAI v1 was prone to occasional large misestimations of AAI, usually as the result of two genomes sharing a set of SCPs which had similar biases towards over or underestimating AAI

FastAAI v2 rectifies these issues:

* Any set of SCPs can be used, so long as they are in HMM (for prokaryotes, see HMMER 3.0 documentation) or block HMM (for eukaryotes, see AUGUSTUS documentation) formats
* Any microbial organism, e.g. prokaryotes or fungi
* Per-SCP models for translating Jaccard similarity to AAI are allowed; these models include attenuation for genome length, SCP match quality
* User-defined labelling of collections of genomes within a database, e.g. bacteria and archaea, with support for post-processing to further increase AAI estimation accuracy
* Taxonomic labelling if deisred

FastAAI v2 consists of 5 steps: database initialization, preprocessing, database construction, database finalization, and database searching
