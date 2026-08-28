# Benchmark datasets

The datasets used in the survey are public GEO series and are not redistributed
here. Download each series and arrange it as below; `setDatasetPaths.m` expects
exactly this layout.

## Directory layout

Each dataset lives in its own folder under `input/`, holding up to seven files:

```
input/<DatasetName>_<GSE>/
├── mix.txt                     # mixture matrix M   (genes x samples)
├── sig.txt                     # signature matrix G (genes x cell types)
├── pure.txt                    # pure profiles H    (genes x pure samples)
├── coef.txt                    # ground-truth proportions C (cell types x samples)
├── feature_annotations.txt     # per-gene annotation
├── reference_annotations.txt   # per-cell-type annotation, with Subclass and Class
└── pure_annotations.txt        # per-pure-sample annotation
```

## Expected folder names

| `setDatasetPaths` name | Folder |
| --- | --- |
| `LiverBrainLung` | `LiverBrainLung_GSE19830` |
| `BreastBlood` | `BreastBlood_GSE29832` |
| `CellLines`, `CellLines_minimal` | `CellLines_GSE11058` |
| `RatBrain` | `RatBrain_GSE19380` |
| `MAQC` | `MAQC_GSE5350` |
| `Retina` | `Retina_GSE33076` |
| `PERT_Cultured` | `PERT_Cultured_GSE16589` |
| `PERT_Uncultured` | `PERT_uncultured_GSE40830` |
| `PBMC_LM22` | `PBMC_Alizadeh` |

## File format

`mix.txt`, `sig.txt`, and `pure.txt` are tab-separated matrices with a header row
of sample or cell-type identifiers and a first column of gene identifiers, read
by `code/utils/my_tblread.m`.

Two things to watch for:

- **Do not log-transform.** `loadDataset` warns and applies `2^x` if the maximum
  value is at most 20, on the assumption that the matrix was log-transformed by
  mistake. Supply values on the linear scale.
- **Gene identifiers must agree** between `mix.txt` and `sig.txt`.
  `loadDataset` intersects them and silently keeps only the overlap, so a
  mismatch in identifier namespace shows up as an unexpectedly small feature set.

`reference_annotations.txt` needs `Reference_ID`, `Subclass`, and `Class`
columns. If the file is absent, `loadDataset` falls back to treating every
reference as its own class, which disables the grouped marker selection in
`select_features`.

`coef.txt` holds the known mixing proportions, used only for evaluation.
