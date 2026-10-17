cd "$(dirname "$0")"
for qd in 0.65 0.66 0.67 0.68; do
  for qn in 0.65 0.66 0.67 0.68; do python3 thq2.py $qd $qn & done; wait
done
