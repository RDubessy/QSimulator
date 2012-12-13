#!/bin/sh
bin=qsimu
if [ $# -eq 0 ]; then
    echo "Usage: init.sh config_file output_dir optional_arg optional_arg"
    exit 0
fi
cfg=$1
if [ $# -gt 1 ]; then
    dir=`pwd`/$2
    if [ ! -d $dir ]; then
        mkdir $dir
    fi
else
    dir=`pwd`
fi
optarg=$3
optarg2=$4
$bin $cfg --groundstate --out=$dir/psi0.dat --general::dttol=0.9 $optarg $optarg2
for i in `seq 2 6`
do
    j=`echo 'scale='$i';1.0-10^(-'$i')' | bc`
    $bin $cfg --groundstate --in=$dir/psi0.dat --out=$dir/psi0.dat --general::dttol=$j $optarg2
done
# init.sh
