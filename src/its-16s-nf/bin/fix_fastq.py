#!/usr/bin/env python3

#
# Fix Errors in a FASTQ file
# -
#

import gzip
from itertools import islice

import click


@click.command()
@click.option("--in-fastq-path", "-i", default="/dev/stdin", help="Input fastq file")
@click.option("--out-fastq-path", "-o", default="/dev/stdout", help="Output fastq file")
def fix_fastq(in_fastq_path, out_fastq_path):
    with gzip.open(in_fastq_path, "rt") as in_file, gzip.open(
        out_fastq_path, "wt"
    ) as out_file:
        good_qual_char = "I"
        while True:
            read = list(islice(in_file, 4))

            if len(read) == 0:
                break

            name, seq, plus, qual = read

            plus = "+"
            name = name.strip()
            qual = qual.strip()
            seq = seq.strip()

            if len(seq) < len(qual):
                # remove surplus qual
                qual = qual[0 : len(seq)]

            if len(qual) < len(seq):
                # add good qual
                qual = qual + good_qual_char * (len(seq) - len(qual))

            fixed_read = "\n".join([name, seq, plus, qual])

            out_file.write(fixed_read)
            out_file.write("\n")


if __name__ == "__main__":
    fix_fastq()
