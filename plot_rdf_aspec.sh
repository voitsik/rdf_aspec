#!/bin/bash

### Setup

# Number of channels
ch=2048

# Number of spectra to integrate
int=100000

###

if [ -z "$1" ]; then
    echo "Usage: $(basename $0) <FILE.rdf> [Band1] [Band2]"
    exit 66
fi

rdf_file=$1

# Just to show on the plot
band1=$2
band2=$3

name=$(basename $rdf_file .rdf)

echo "Integration time: $(echo "scale=3; 2 * 31.25 * $ch * $int / 10^9" | bc -lq) s"
rdf_aspec $rdf_file $ch $int > ${name}.sp || exit 1


cat > "${name}.gnuplot" << EOF
set term postscript color solid
set output "${name}.ps"

set xlabel "Frequency (MHz)"
set title "$(sed -n '1 s/#//p' $name.sp)"
set key box center bottom

plot [-16:16][0:2] \
    '${name}.sp' u 1:2      w l lc 1 notitle, 0/0 w l lc 1 lw 5 title '$band1 (USB)', \
    '${name}.sp' u (-\$1):3 w l lc 2 notitle, 0/0 w l lc 2 lw 5 title '$band1 (LSB)', \
    '${name}.sp' u 1:4      w l lc 3 notitle, 0/0 w l lc 3 lw 5 title '$band2 (USB)', \
    '${name}.sp' u (-\$1):5 w l lc 4 notitle, 0/0 w l lc 4 lw 5 title '$band2 (LSB)'
EOF

gnuplot < ${name}.gnuplot

#ps2pdf "${name}.ps"
#rm -f "${name}.ps"
