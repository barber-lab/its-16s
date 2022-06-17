#!/usr/bin/env python3

import sys
import msea
from msea import SetLibrary

db_gmt_path = sys.argv[1]
test_taxa_path = sys.argv[2]
out_csv_path = sys.argv[3]
universe = int(sys.argv[4])

db_gmt = msea.read_gmt(db_gmt_path)

with open(test_taxa_path) as f:
  test_taxa = [x.strip() for x in f.readlines()]

msea_result = msea.enrich(test_taxa, d_gmt=db_gmt, universe=universe)
msea_result.to_csv(out_csv_path)
