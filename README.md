# KrakMeOpen

A [Kraken 2] downstream analysis toolkit. More specifically, calculate a series of quality metrics for Kraken 2 classifications.

## Installation
KrakMeOpen is available to install through conda. Simply run the following command to install it:

`conda install -c conda-forge -c bioconda krakmeopen`

## Usage

A good start is to run `krakmeopen --help`.

To calculate quality metrics for a Kraken 2 classification, run:

`krakmeopen --input classifications.kraken2 --output metrics_out.tsv --names names.dmp --nodes nodes.dmp [--tax_id INT | --tax_id_file tax_ids.txt] --output_kmer_tally kmer_tally_out.tsv`

Where:
* _input_ is the read-by-read classifications output by Kraken 2 (or [StringMeUp]). **Required**.
* _output_ is the file to write the output to. **Required**.
* _name_ and _nodes_ are the same NCBI style taxonomy files used for the building of the database that was used to produce the original Kraken 2 classifications. **Required**.
* _tax_id_ takes a single taxonomic ID (taxID), while _tax_id_file_ is a file containing multiple taxIDs (one per line) that you wish to get quality metrics for.
Must specify one of the two. **Required**.
* _output_kmer_tally_ is a file to output the complete kmer tally for each tax ID to (human readable). Optional.

KrakMeOpen consists of two steps: (1) tally kmers, and (2) calculate metrics. It is possible to do only the first step and save the result in a pickle for later use. This is
for example when you want to calculate metrics from multiple classifications. The kmer tallies from multiple pickles are added together and metrics are calculated on that sum of kmers.

In the below example, calculate combined metrics for the clade rooted at human (9606) from two kraken 2 classifications (sample_1, sample_2). First tally the kmers from sample_1 and sample_2 by specifying --output_pickle. Then, calculate metrics and output the final result to combined_metrics.tsv. Note that the tax_id, names.dmp, and nodes.dmp must be the same in all calls.

`krakmeopen --input sample_1.kraken2 --output_pickle sample_1_metrics.pickle --names names.dmp --nodes nodes.dmp --tax_id 9606`

`krakmeopen --input sample_2.kraken2 --output_pickle sample_2_metrics.pickle --names names.dmp --nodes nodes.dmp --tax_id 9606`

`krakmeopen --input_file_list pickle_file_list.txt --output combined_metrics.tsv --names names.dmp --nodes nodes.dmp --tax_id 9606`

## Quality metrics

The metrics are calculated on the clade-level. All kmers from all reads that are classified to any of the nodes in the
clades rooted at the supplied tax IDs are aggregated, and metrics are calculated on those aggregations. Input is
Kraken2 read-by-read classification files (can be gzipped). Output is a tab separated file containing the metrics.

The following metrics are calculated:

| Metric name | Description |
|-------------|-------------|
| nkmers_total | Total number of kmers |
| nkmers_classified | Total number of classified kmers |
| nkmers_unclassified | Total number of unclassified kmers |
| nkmers_clade | Total number of kmers classified to any tax ID within the clade |
| nkmers_lineage | Total number of kmers classified to any tax ID directly above the clade root tax ID |
| confidence_origin | The confidence score for the clade, calculated as described by Kraken2 |
| confidence_classified | An alternative confidence score where the unclassified kmers are removed from the denominator |
| other_kmers_lineage_ratio | Ratio of nkmers_lineage / (nkmers_total - nkmers_clade) |
| other_kmers_root_ratio | Ratio of "kmers classified to root" / (nkmers_total - nkmers_clade) |
| other_kmers_classified_ratio | Ratio of (nkmers_total - nkmers_clade - nkmers_unclassified) / (nkmers_total - nkmers_clade) |
| other_kmers_distance | Average distance between the clade root tax ID and the tax IDs which kmers are classified to |
| other_kmers_distance_lineage_excluded | Like other_kmers_distance but kmers classified to tax IDs above the clade are excluded |

[Kraken 2]: https://github.com/DerrickWood/kraken2
[StringMeUp]: https://github.com/danisven/stringmeup
