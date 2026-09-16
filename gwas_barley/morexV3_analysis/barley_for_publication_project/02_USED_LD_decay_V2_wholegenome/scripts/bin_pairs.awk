# bin_pairs.awk
# Per-worker partial aggregation. Reads PLINK .ld stream from stdin (no header
# expected — the launcher strips it) and emits one row per encountered bin:
#   bin <TAB> SUM(r^2) <TAB> SUM(r^2^2) <TAB> COUNT
# Parameters (set with -v):
#   MAXD = maximum physical distance to keep (bp)
#   BSZ  = bin size (bp)
{
    d = ($5 > $2) ? $5 - $2 : $2 - $5;
    if (d == 0 || d > MAXD) next;
    bin = int(d / BSZ) * BSZ;
    sum[bin]   += $7;
    sumsq[bin] += $7 * $7;
    n[bin]++;
}
END {
    OFS = "\t";
    for (b in n) print b, sum[b], sumsq[b], n[b];
}
