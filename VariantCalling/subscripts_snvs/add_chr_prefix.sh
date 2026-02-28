#!/bin/bash

input=$1
output=$2

# check if non-comment lines (comment lines start with #) start with "chr"
# and prepend "chr" if they don't

if [[ "$input" == *.gz ]]; then
    gunzip -c "$input"
else
    cat "$input"
fi | awk 'BEGIN { OFS="\t" }
     /^#/ { print; next }
     $1 !~ /^chr/ { $1 = "chr"$1 }
     { print }
' > "$output"

# if [[ "$input" == *.gz ]]; then
#     gunzip -c "$input"
# else
#     cat "$input"
# fi | awk 'BEGIN { OFS="\t" }
#      /^#/ { print; next }
#      $1 !~ /^chr/ { $1 = "chr"$1 }
#      { print }
# ' | bgzip -c > "$output"
