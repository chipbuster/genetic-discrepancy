#!/bin/bash

if [ ! -f knl-discrep.x ]; then
  echo "No KNL binary! Please run script to build KNL binaries"
  exit 1
fi
if [ ! -f sb-discrep.x ]; then
  echo "No SB binary! Please run script to build SB binaries"
  exit 1
fi

find knl-tests -type d | while read dir; do
  cp knl-discrep.x "$dir"
done

find sb-tests -type d | while read dir; do
  cp sb-discrep.x "$dir"
done
