# Benchmark datasets

> These files are stored with [Git LFS](https://git-lfs.com). Run
> `git lfs install` before cloning, or `git lfs pull` in an existing clone,
> or you will get pointer stubs instead of data.

Seven of the nine datasets used in the survey ship with this repository. Two do
not, and must be obtained separately. See the Data section of the top-level
[README](../README.md) for the accession numbers of both groups.

`setDatasetPaths.m` maps each dataset name to the files below, so adding a
dataset of your own means creating a folder in this layout and adding one
`case` there.

## Directory layout

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

`CellLines_GSE11058` additionally carries `mix_minimal.txt`, `sig_minimal.txt`,
and `pure_minimal.txt`, the reduced version selected by the
`CellLines_minimal` dataset name.

## Missing folders

| Dataset name | Expected folder |
| --- | --- |
| `MAQC` | `MAQC_GSE5350` |
| `PBMC_LM22` | `PBMC_Alizadeh` |

## File format

`mix.txt`, `sig.txt`, and `pure.txt` are tab-separated matrices with a header
row of sample or cell-type identifiers and a first column of gene identifiers,
read by `code/utils/my_tblread.m`.

Two things to watch for when adding your own data:

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
