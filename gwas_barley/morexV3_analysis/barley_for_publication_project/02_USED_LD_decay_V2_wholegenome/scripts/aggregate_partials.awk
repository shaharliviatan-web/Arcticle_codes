# aggregate_partials.awk
# Reads the concatenated partial outputs from many bin_pairs.awk workers
# (no header) and emits the final per-bin sums:
#   Bin <TAB> Sum_R2 <TAB> SumSq_R2 <TAB> N_Pairs   <- header
#   <bin> <TAB> ...
{
    sum[$1]   += $2;
    sumsq[$1] += $3;
    n[$1]     += $4;
}
END {
    OFS = "\t";
    print "Bin", "Sum_R2", "SumSq_R2", "N_Pairs";
    for (b in n) print b, sum[b], sumsq[b], n[b];
}
