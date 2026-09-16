# 01_assoc — shared GWAS association files

One file per trait, `SNP` + `P`, 7,110,996 data rows each (7,110,997 with header),
extracted from the EMMAX `.ps` output for the locked configuration (BLUP × 3 PCs).

These are **shared GWAS output, not owned by any one locus definition**, which is why
they sit at this level rather than inside a results directory.

The locus script verifies the full row count before running and
aborts if any file is short:

```bash
NEED=7110997
[[ -s $A && "$(wc -l < "$A")" -eq $NEED ]] || { echo "[31] FAIL: $A incomplete"; exit 1; }
```

Regenerate from `results/emmax_ps/morexV3__<trait>__BLUP__pc3.ps` if ever needed.
