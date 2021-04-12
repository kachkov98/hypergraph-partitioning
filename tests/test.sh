#! /bin/bash
for f in tests/ISPD/*
do
  echo "$f:"
  ./$1 $f $2 | grep "Final cut cost: "
done
