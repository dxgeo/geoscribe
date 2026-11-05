#!/bin/sh
# Rerun NGS's own gridding command (gmtbat02's tension-0.4 line) on NGS's
# thinned points, one band at a time, and export each grid as float64 for
# cmp_surface.py and nadcon5_distance.py.
# Usage: sh nadcon5_surface.sh WORK      (needs GMT; about a second a band)
set -eu
W=$1; D=$W/nadcon5-download; G=$W/gmt-work; T=nad83_2007.nad83_2011.conus
mkdir -p "$G"
for b in lat lon eht; do
    gmt surface "$D/smtcd$b.$T.60" -R235/294/24/50 -I1m -T0.4 -A0.7986 -C0.01 -G"$G/$b.nc"
    gdal_translate -q -of ENVI -ot Float64 "$G/$b.nc" "$G/$b.bin"
done
