#/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

# create output folder for initial run
mkdir -p Data
mkdir -p Figure
mkdir -p out
# start the initial run
 docker run -it --rm -v $PWD:/run partmc_with_mosaicc bash -c 'cd /run;/build/partmc urban_plume.spec'


