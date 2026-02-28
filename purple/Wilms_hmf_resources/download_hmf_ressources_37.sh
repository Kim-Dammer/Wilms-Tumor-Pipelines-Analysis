#!/bin/bash

TARGET=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
URL="https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.3.0--2.tar.gz"

if [ "$(ls -A "$TARGET" | grep -vE "(.gitkeep|$(basename "$0"))" | wc -l)" -gt 0 ]; then
    echo "Data already present in $TARGET. Aborting."
    exit 0
fi

echo "Fetching HMF resources into $TARGET..."
wget -qO- "$URL" | tar -xzf - --strip-components=1

echo "Done."
