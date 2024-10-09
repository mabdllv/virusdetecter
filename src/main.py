from argparse import ArgumentParser
from subprocess import call
import os

from analyze_kraken_output import *


def run_kraken(db, classified_data_dir, kraken_output_dir, input_file1, input_file2, num_threads):
    if not (os.path.exists(classified_data_dir)):
        os.makedirs(classified_data_dir)

    if not (os.path.exists(kraken_output_dir)):
        os.makedirs(kraken_output_dir)

    out_fq_file = os.path.join(classified_data_dir, "cseqs#.fq")
    out_csv_file = os.path.join(kraken_output_dir, "kraken_output.csv")
    out_full_csv_file = os.path.join(kraken_output_dir, "kraken_full_output.csv")

    call(["kraken2",
          "--threads", num_threads,
          "--db", db,
          "--paired", "--use-names",
          "--classified-out", out_fq_file,
          input_file1, input_file2,
          "--output", out_csv_file,
          "--report", out_full_csv_file])


def run_spades(input_dir, output_dir, threads):
    for taxon_dir in os.listdir(input_dir):
        taxon_id = taxon_dir[taxon_dir.find("_") + 1:]

        taxon_output_dir = os.path.join(output_dir, f"taxon_{taxon_id}")
        if not (os.path.exists(taxon_output_dir)):
            os.makedirs(taxon_output_dir)

        fq_1 = os.path.join(input_dir, taxon_dir, f"taxon_{taxon_id}_1.fq")
        fq_2 = os.path.join(input_dir, taxon_dir, f"taxon_{taxon_id}_2.fq")

        call(["spades.py",
              "-1", fq_1,
              "-2", fq_2,
              "-o", taxon_output_dir,
              "-t", threads,
              "-k", 21, 33, 55, 77])


def run_blast(input_dir, output_dir):
    for taxon_dir in os.listdir(input_dir):
        scaffolds_file_path = os.path.join(input_dir, taxon_dir, "scaffolds.fasta")
        if not (os.path.exists(scaffolds_file_path)):
            continue

        taxon_id = taxon_dir[taxon_dir.find("_") + 1:]

        taxon_output_file = os.path.join(output_dir, f"taxon_{taxon_id}.out")

        call(["blastn",
              "-db","nt",
              "-query", scaffolds_file_path,
              "-out", taxon_output_file,
              "-remote"])


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-p", "--min_percent", default=0.0005,
                        help="minimum valuable percent of concrete virus sequences among all virus sequences in genome")
    parser.add_argument("-t", "--threads", default=4,
                        help="number of threads")
    parser.add_argument("-db", "--kraken_db",
                        help="absolute path to kraken viral_db")
    parser.add_argument("-out", "--output", default=os.getcwd(),
                        help="absolute path to virusdetecter output directory")
    parser.add_argument('filename1')
    parser.add_argument('filename2')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    classified_data_dir = os.path.join(args.output, "work_data","classified_data")
    kraken_output_dir = os.path.join(args.output, "work_data", "kraken_output")
    taxon_fq_dir = os.path.join(args.output, "work_data", "taxon_fq")
    spades_output_dir = os.path.join(args.output, "work_data", "spades_output")
    blast_output_dir = os.path.join(args.output, "work_data", "blast_output")

    run_kraken(args.kraken_db, classified_data_dir, kraken_output_dir, args.filename1, args.filename2, args.threads)

    analyze_kraken_output(os.path.join(kraken_output_dir, "kraken_full_output.csv"),
                          os.path.join(classified_data_dir, "kraken_plot.png"),
                          os.path.join(classified_data_dir, "kraken_taxon_ids_filtered.txt"),
                          args.min_percent if args.min_percent > 0 else 0.0005)

    parse_fq_by_taxon(os.path.join(classified_data_dir, "cseqs1.fq"),
                      os.path.join(classified_data_dir, "cseqs2.fq"),
                      os.path.join(classified_data_dir, "kraken_taxon_ids_filtered.txt"),
                      taxon_fq_dir)

    run_spades(taxon_fq_dir, spades_output_dir, args.threads)

    run_blast(spades_output_dir, blast_output_dir)
