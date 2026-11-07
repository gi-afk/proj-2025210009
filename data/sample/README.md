# Sample data (offline smoke test only)

`sample_dataset.csv` contains 45 windows (15 positives, 15 random negatives,
15 distance-matched negatives) drawn from the full dataset:

- genome: *S. cerevisiae* S288C R64 / sacCer3 (NCBI RefSeq GCF_000146045.2)
- origins: OriDB `sc_ori`, status = Confirmed
- window: 500 bp centred on the annotated origin midpoint

It exists so that `python scripts/smoke_test.py` and the unit tests can run
without network access. It is **not** sufficient for the reported results and
must not be used for AUC estimates or biological conclusions.

The full dataset is rebuilt in `data/processed/dataset.csv` by
`scripts/run_all.py`; raw downloads are git-ignored. See `docs/data-card.md`.
