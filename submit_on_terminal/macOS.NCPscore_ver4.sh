#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

LOG="NCPscore_ver4.log"

for i in {17..22} X Y M; do
  echo "[$(date)] Start to process chr${i}" >> "$LOG"

  find . -name "chr${i}_*Samp11*_bin.gtab.gz" -print0 |
  while IFS= read -r -d '' f; do
    BASENAME="$(basename -- "$f")"
    PREFIX="${BASENAME/_bin.gtab.gz/}"   # output prefix like chr1_Exp2-Samp11_S5

    echo "$PREFIX"  # optional: stdout

    # build control pattern: replace Samp11 -> Samp1
    CONTROL_PRE="${PREFIX/Samp11/Samp1}"

    # remove the trailing _S<number> from CONTROL_PRE
    CONTROL_LOC="${CONTROL_PRE%_S*}"     # e.g. chr1_Exp2-Samp1

    CONTROL=( ./"${CONTROL_LOC}"_S*_bin.gtab.gz )

    # if (( ${#CONTROL[@]} == 0 )); then
    #   echo "[$(date)] WARNING: No control file for $PREFIX (pattern: ${CONTROL_LOC}_S*_bin.gtab.gz). Skipping." >> "$LOG"
    #   continue
    # fi

    # If you expect exactly one control file, enforce it:
    if (( ${#CONTROL[@]} != 1 )); then
      echo "[$(date)] ERROR: Expected 1 control, found ${#CONTROL[@]} for $PREFIX: ${CONTROL[*]}" >> "$LOG"
      continue
    else
      echo "[$(date)] Start to process $f (control: ${CONTROL[*]})" >> "$LOG"

      nice -n 3 python3 ../../condense-seq/prepro_scripts/NCPscore_ver4.py \
        -f "$f" \
        -i "${CONTROL[0]}" \
        -o "$PREFIX" \
        | tee -a "$LOG" && wait
    fi
  done
done