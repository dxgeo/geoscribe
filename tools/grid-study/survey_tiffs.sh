#!/bin/sh
# What D-T2's reader must handle and what D-T7's licences say: for every grid
# in PROJ_DATA, its compression, predictor, interleave, declared interpolation,
# and copyright statement, tallied; then PROJ-data's own licence table.
# Usage: PROJ_DATA=/path/to/proj sh survey_tiffs.sh
set -eu
P=${PROJ_DATA:?set PROJ_DATA}
for f in "$P"/*.tif; do gdalinfo -json "$f"; echo; done | python3 -I -c '
import json, sys, collections
fmt = collections.Counter(); lic = collections.Counter(); buf = ""
for line in sys.stdin:
    if line.strip(): buf += line; continue
    if not buf: continue
    j = json.loads(buf); buf = ""
    m = j.get("metadata", {}); s = m.get("IMAGE_STRUCTURE", {}); d = m.get("", {})
    fmt[(s.get("COMPRESSION", "none"), s.get("PREDICTOR", "-"), s.get("INTERLEAVE", "-"), d.get("interpolation_method", "bilinear"))] += 1
    lic[d.get("TIFFTAG_COPYRIGHT", "(none)")] += 1
for k, n in fmt.most_common(): print(n, *k)
print()
for k, n in lic.most_common(): print(n, k)
'
echo; cut -d, -f3 "$P/copyright_and_licenses.csv" | sort | uniq -c | sort -rn
