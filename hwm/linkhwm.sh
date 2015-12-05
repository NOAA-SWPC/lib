#!/bin/sh

hwmdir="../hwm"

for i in $hwmdir/*.dat $hwmdir/*.bin; do
  f=$(basename $i)
  ln -sf $i $f
done
