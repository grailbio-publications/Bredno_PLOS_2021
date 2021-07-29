# Command line script

This directory contains the script `predict_ctf.R` that provides model predictions
from user-provided input. This directory receives the prediction models for
breast, lung, and colorectal cancer when the markdown documents in the directory
Rmd are rendered. These models are necessary to run the script.

The script first queries the cancer type for a model:

```
Cancer type
[b]reast cancer / [l]ung cancer [c]olorectal cancer
```

The entry must be "b", "l", or "c", all other entries result in an exit with
error.

The script will then query the number of primary tumor lesions (this is mostly
1):

```
Number of primary lesions: 
```

Entries between 1 and 6 are valid. The script will then query the size of each
primary tumor lesion:

```
Max. size of lesion 1 in mm (n for unknown): 
```

For multiple tumor lesions, the entered data is expected ordered by size, with
the largest / index lesion first. The size in mm must be entered for the first
(index) lesion. The entry "n" for unknown size is valid for evtl. additional tumor
lesions.

Depending on the selected cancer type, the script will then query %Ki-67 positive
for breast cancer (valid are 0 - 100), the FDG SUV_max of the tumor for a
lung cancers (valid values are above 1), or the depth of microinvasion for a
colorectal cancer:

```
[0] Other/missing
[1] Intramucosal
[2] Submucosa
[3] Muscularis propria
[4] Subserosa
[5] Penetrates serosa
[6] Invades adjacent
```

The entry 0 here is valid, then the script will query the T-stage of the
colorectal cancer with valid values 1-4 for T1, T2, T3, T4 respectively.

The script will end with the output of the models (example here):

```
The model predicts a circulating tumor fraction of 0.000391
  (95% CI [0.000206; 0.000575])
Estimated number of tumor-derived GEs in circulation: 1334.9
  (95% CI [687.4; 1982.5])
```

Estimated circulating tumor fraction and total number of genome equivalents for tumor-
derived cfDNA fragments in the circulation of the patient are given together with
their respecitve 95% confidence intervals.
