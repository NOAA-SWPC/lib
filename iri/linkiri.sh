#!/bin/sh

iridir="../iri"

for i in $iridir/*.dat $iridir/*.asc; do
  f=$(basename $i)
  ln -sf $i $f
done
