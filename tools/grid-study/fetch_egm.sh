#!/bin/sh
# Download EGM96 and EGM2008's spherical-harmonic coefficients from ICGEM
# (GFZ; NGA's earth-info.nga.mil timed out on 2026-11-05) into a directory of
# their own, then parse them to .npz for egm96_test.py and egm2008_*.py.
# Downloaded files are data: they go in WORK/egm-download, the scripts stay here.
# Usage: sh fetch_egm.sh WORK
set -eu
W=$1; D=$W/egm-download; mkdir -p "$D" "$W/egm-work"
B=https://icgem.gfz.de/getmodel/zip
curl -sfL -o "$D/EGM96.zip" "$B/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.zip"
curl -sfL -o "$D/EGM2008.zip" "$B/c50128797a9cb62e936337c890e4425f03f0461d7329b09a8cc8561504465340/EGM2008.zip"
(cd "$D" && unzip -o -q EGM96.zip && unzip -o -q EGM2008.zip)
H=$(dirname "$0")
python3 -I "$H/egm_parse.py" "$D/EGM96.gfc" "$W/egm-work/egm96.npz"
python3 -I "$H/egm_parse.py" "$D/EGM2008.gfc" "$W/egm-work/egm2008.npz"
