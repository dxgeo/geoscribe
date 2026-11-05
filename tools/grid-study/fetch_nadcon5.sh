#!/bin/sh
# Download the parts of NGS's NADCON5 build (20160901 release) the source
# study reads, for NAD83(2007) -> NAD83(2011), CONUS: the control file, the
# grid parameters, the thinned surface inputs (smtcd*), one all-points file,
# the statistics, the log and the gridding batch files (gmtbat02/05).
# Usage: sh fetch_nadcon5.sh WORK
set -eu
W=$1; D=$W/nadcon5-download; mkdir -p "$D"
B=https://geodesy.noaa.gov/pub/nadcon5/20160901release; T=nad83_2007.nad83_2011.conus
cd "$D"
for f in RunNADCON5/Control/control.$T RunNADCON5/Data/grid.parameters \
         Builds/$T/dvstats.$T.60 Builds/$T/final.$T.log Builds/$T/cvacdlat.$T \
         Builds/$T/gmtbat02.$T.60 Builds/$T/gmtbat05.$T.60 \
         Builds/$T/smtcdlat.$T.60 Builds/$T/smtcdlon.$T.60 Builds/$T/smtcdeht.$T.60; do
    curl -sfL -O "$B/$f.gz"
done
gunzip -f ./*.gz
