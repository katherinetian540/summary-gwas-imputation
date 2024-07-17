__author__ = "alvaro barbeira"

import logging
import os
import re
import sqlite3
import pandas as pd
import numpy as np
import gzip
import pyarrow.parquet as pq

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.miscellaneous import matrices, PandasHelpers
from genomic_tools_lib.file_formats import Parquet

def get_file_map(args):
    logging.log(9, "Loading parquet files")
    r = re.compile(args.parquet_genotype_pattern)
    files = os.listdir(args.parquet_genotype_folder)
    print("Files in directory:", files)  # Debugging line
    files = {r.search(f).group(1): os.path.join(args.parquet_genotype_folder, f) for f in files if r.search(f)}
    print("Files found:", files)  # Debugging line
    p = {}
    for k, v in files.items():
        logging.log(9, "Loading %s: %s", k, v)
        dataset = pq.ParquetDataset(v)
        p[k] = dataset
    return p

n_ = re.compile("^(\d+)$")

def process_gene_data(g_, w, file_map, individuals, output_file, args):
    chroms = set(w.varID.str.split("_").str[0].str.replace("chr", ""))
    chrom_files = {chrom: file_map[chrom] for chrom in chroms if chrom in file_map}
    
    # Combine data from all chromosomes
    d = {}
    for chrom, dataset in chrom_files.items():
        columns = w[w.varID.str.contains(f"chr{chrom}_")].varID.values.tolist()  # Convert NumPy array to list
        d.update(Parquet._read(dataset, columns=columns))

    if not d:
        logging.log(9, "No genotype available for %s, skipping", g_)
        return

    var_ids = list(d.keys())
    if args.output_rsids:
        ids = [x for x in pd.DataFrame({"varID": var_ids}).merge(w[["varID", "rsid"]], on="varID").rsid.values]
    else:
        ids = var_ids

    c = np.cov([d[x] for x in var_ids])
    c = matrices._flatten_matrix_data([(w.gene.values[0], ids, c)])
    result = []
    for entry in c:
        l = "{} {} {} {}\n".format(entry[0], entry[1], entry[2], entry[3])
        result.append(l)
    return result

def run(args):
    if os.path.exists(args.output):
        logging.info("Output already exists, either delete it or move it")
        return

    logging.info("Getting parquet genotypes")
    file_map = get_file_map(args)
    print("File map keys:", list(file_map.keys()))  # Debugging line

    logging.info("Getting genes")
    with sqlite3.connect(args.model_db) as connection:
        extra = pd.read_sql("SELECT * FROM EXTRA ORDER BY gene", connection)
        extra = extra[extra["n_snps_in_model"] > 0]
    print("Genes loaded:", extra.head())  # Debugging line

    individuals = TextFileTools.load_list(args.individuals) if args.individuals else None

    logging.info("Processing")
    Utilities.ensure_requisite_folders(args.output)

    results = []
    with sqlite3.connect(args.model_db) as connection:
        for i, t in enumerate(extra.itertuples()):
            g_ = t.gene
            logging.log(9, "Processing %i/%i: %s", i+1, extra.shape[0], g_)
            w = pd.read_sql("SELECT * FROM weights WHERE gene = '{}';".format(g_), connection)
            print("Weights for gene {}: {}".format(g_, w.head()))  # Debugging line
            result = process_gene_data(g_, w, file_map, individuals, args.output, args)
            results.append(result)

    # Compute all results and write to the output file
    with gzip.open(args.output, "w") as f:
        f.write("GENE RSID1 RSID2 VALUE\n".encode())
        for result in results:
            for line in result:
                f.write(line.encode())

    logging.info("Finished building covariance.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate BSLMM runs on study")
    parser.add_argument("-parquet_genotype_folder", help="Parquet Genotype folder")
    parser.add_argument("-parquet_genotype_pattern", help="Pattern to detect parquet genotypes by chromosome")
    parser.add_argument("-model_db", help="Where to save stuff")
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("--output_rsids", action="store_true")
    parser.add_argument("--individuals")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default=10)
    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)
