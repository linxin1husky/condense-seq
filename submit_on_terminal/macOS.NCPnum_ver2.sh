#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

LOG="NCPnum_ver2.log"
TIRATION_FILE="/Volumes/BackupHaLab/titration_files/mCD8T:WT_NCP_sp_titration.csv"

for i in {1..10}; do
  echo "[$(date)] Start to process chr${i}" >> "$LOG"

  find . -name "chr${i}_*Samp1*_bin.gtab.gz" -print0 |
  while IFS= read -r -d '' f; do
    BASENAME="$(basename -- "$f")"
    PREFIX="${BASENAME/_bin.gtab.gz/}"   # output prefix like chr1_Exp2-Samp11_S5
    
    # Only process the non-existing file
    if [ ! -f ${PREFIX}_num.gtab.gz ]; then
      echo "[$(date)] Start to process $f" >> "$LOG"

      nice -n 3 python3 ../../condense-seq/prepro_scripts/NCPnum_ver2.py \
        -f "$f" \
        -t ${TIRATION_FILE}\
        -o "$PREFIX" \
        | tee -a "$LOG" && wait
    else
      echo -e "[$(date)] The output ${PREFIX}_num.gtab.gz has already existed. Skip to next file."
    fi
  done
done