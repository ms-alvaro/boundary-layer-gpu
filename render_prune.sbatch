#!/bin/bash
#SBATCH --account=marlowe-m000252
#SBATCH --partition=preempt
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --job-name=mz_renderprune
#SBATCH --output=/scratch/m000252/tbl_mz_ly55/render_prune.log
# Render frames for any snapshot without one, then prune rendered snapshots
# (keep every 10000th step + the 3 newest as restart safety).
cd /scratch/m000252/tbl_mz_ly55
python3 movie_fst.py data --fps 8 --stride 1 --vmax 0.20 || { echo RENDER_FAILED; exit 1; }
echo "RENDER_OK frames=$(ls data/_frames/s*.png 2>/dev/null | wc -l)"
snaps=$(ls data/BL_prod.[0-9]*[0-9] 2>/dev/null | grep -v stats | sort -t. -k2 -n)
n=$(echo "$snaps" | wc -l)
i=0
for f in $snaps; do
  i=$((i+1))
  step=$(basename "$f" | sed 's/.*\.//')
  stepn=$((10#$step))
  [ $((stepn % 10000)) -eq 0 ] && continue          # keep every 10000th
  [ $i -gt $((n-3)) ] && continue                   # keep 3 newest
  [ -f "data/_frames/s${step}.png" ] && rm -f "$f" && echo "pruned $f"
done
echo "PRUNE_DONE used=$(du -sh /scratch/m000252/tbl_mz_ly55 | cut -f1)"
