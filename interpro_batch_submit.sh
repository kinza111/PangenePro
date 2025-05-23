#!/bin/bash

email="$1"
base_dir="$2"
waitevery=30

if [ -z "$email" ] || [ -z "$base_dir" ]; then
  echo "Usage: $0 <email> <seqs_directory>"
  exit 1
fi

count=0
for fa_file in "${base_dir}"/seqRefSet_*/**/*.fa "${base_dir}"/seqRefSet_*/*.fa; do
  [ -f "$fa_file" ] || continue

  xml_file="${fa_file%.fa}.xml"
  if [ -f "$xml_file" ]; then
    echo "Skipping already processed: $xml_file"
    continue
  fi

  echo "Submitting: $fa_file"
  python3 iprscan_client.py \
    --sequence "$fa_file" \
    --email "$email" &

  ((count++))
  if (( count % waitevery == 0 )); then
    wait
  fi
done

wait
echo "All jobs submitted."
exit 0
