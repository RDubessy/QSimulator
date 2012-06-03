#!/bin/sh
bin=qsimu
if [ ! $# -eq 2 ]; then
    echo "Usage: movie.sh input_dir z_max"
    exit 0
fi
echo "Converting profiles to png images"
count=0
for i in `ls $1 | grep dat_t`
do
    $bin example.cfg --in=$1/$i --plot
    count=$((count+1))
    filename=$(printf "%04d.png" $count)
    echo "set view map;unset surface;set pm3d;set size square;
        set term png small;set cbrange [0:$2]
        set output \"$1/$filename\";
        splot \"/tmp/psi.txt\" using 1:2:(\$3*\$3+\$4*\$4)" | gnuplot
done
echo "Converting png to jpg"
for i in `ls $1 | grep png`
do
    convert $1/$i $1/$i.jpg
done
echo "Creating the movie"
mencoder "mf://$1/*.jpg" -mf fps=10 -o movie.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800
echo "Cleaning auxiliary files"
rm -f $1/*.png*
