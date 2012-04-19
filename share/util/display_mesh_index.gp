# !/usr/bin/gnuplot -persist
# $Id: display_mesh_index.gp 4945 2009-02-12 05:19:40Z dstrubbe $
#
# replace 'mesh_index' by whatever filename you'd like to use (assumes mesh_index output mode of octopus)
# delete space between hash and ! above and execute standalone
# or load this file inside gnuplot with: load "display_mesh_index.gp"

unset key
unset label
!tail -n +2 mesh_index | awk '(NF > 0) {print "set label \"" $1"\" at "$2+0.07","$3","$4}' > mesh_index.tmp
set xlabel "x"
set ylabel "y"
set zlabel "z"
load "mesh_index.tmp"
set style line 1 lt 3 lw 1 pt 13 ps 1.7
sp 'mesh_index' u 2:3:4 w p ls 1
!rm -f mesh_index.tmp
