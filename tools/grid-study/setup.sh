#!/bin/sh
# Lay out a work directory and export the sample grids the per-file studies
# read: grids.txt (every current grid transformation and its PROJ file), and
# each sample grid as raw float32 *band-sequential* (.bin, forced: NADCON5's
# files are pixel-interleaved, and reading them as band-sequential scrambles
# them) beside its gdalinfo JSON. NTv2 Canada's 114 subgrids go to grids/ca/.
# Usage: PROJ_DATA=/path/to/proj sh setup.sh WORK
set -eu
W=$1; P=${PROJ_DATA:?set PROJ_DATA to the PROJ data directory}
mkdir -p "$W/grids/ca"
sqlite3 -separator '|' "$P/proj.db" "select distinct g.method_name, ga.proj_grid_name from grid_transformation g join grid_alternatives ga on ga.original_grid_name = g.grid_name where g.deprecated = 0;" > "$W/grids.txt"
for g in us_noaa_conus ca_nrc_ntv2_0 uk_os_OSTN15_NTv2_OSGBtoETRS us_nga_egm96_15 us_nga_egm08_25 \
         de_lgl_bw_BWTA2017 us_noaa_nadcon5_nad83_2007_nad83_2011_conus \
         us_noaa_g2018u0 de_bkg_gcg2016 es_ign_egm08-rednap au_ga_AUSGeoid2020_20180201; do
    gdal_translate -q -of ENVI -co INTERLEAVE=BSQ -ot Float32 "$P/$g.tif" "$W/grids/$g.bin"
    gdalinfo -json "$P/$g.tif" > "$W/grids/$g.json"
done
i=0
gdalinfo "$P/ca_nrc_ntv2_0.tif" | sed -n 's/.*SUBDATASET_[0-9]*_NAME=//p' | while read -r sd; do
    i=$((i + 1))
    gdal_translate -q -of ENVI -co INTERLEAVE=BSQ -ot Float32 "$sd" "$W/grids/ca/s$i.bin"
    gdalinfo -json "$sd" > "$W/grids/ca/s$i.json"
done
