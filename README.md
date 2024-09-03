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

FastAAI v2 consists of 5 steps: database initialization, preprocessing, database construction, database finalization, and querying. See the example_prokaryote_files zip for an example of the instructions for building and querying a complete FastAAI v2 database.

## Initialization

fastaai_main init

The initialization step creates an empty databases, flags it as either prokaryotic or eukaryotic and adds HMM models to the database. Optionally, a user can also provide a ruleset for identifying genomes as belonging to user-defined genome groups, post-processing for each pair of groups, and include per-SCP AAI models. 

## Preprocessing

fastaai_main preproc

Uses an initialized database as a repository for SCPs and predicts proteins for each genome in a collection, identifies SCP representatives within those proteins, and extracts all of the information FastAAI uses to construct a database into a file FastAAI calls a "crystal." Crystals are portable, lightweight representations of a genome ready to be added to a FastAAI database in your project or elsewhere. Optionally, add taxonomic information for each genome; this information will be stored in the crystal as well.

## Database construction

fastaai_main consume

Add a collection of crystals to a FastAAI database. The SCPs within each crystal must be the same or a subset of the SCPs within the database; FastAAI will skip any SCPs which are not.

## Database finalization

fastaai_main finalize

Convert the records in a FastAAI database into a faster, query-friendly representation prior to a search.

## Querying

fastaai_main query

Search one FastAAI database against another. At least some SCPs must be shared in common between the two databases for a search to proceed. FastAAI will use the target database, specificed with -tdb in the fastaai_main query command, as the source of per-SCP AAI models and genome labelling rules.
