#!/bin/sh
#this script removes ^M chars from a file.
#Used to process framefinder out file
tmpfile=$(mktemp)
tr -d '\r' <"$1" >"$tmpfile"
mv "$tmpfile" "$1"
