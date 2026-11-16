# chem-artifact-bench

Synthetic internal consistency test, not external validation. LLM baseline/tool-augmented results exist only after an actual live API run.

| Mode | Cases | Accuracy | Artifact precision | Recall | F1 | Artifact acceptance | Errors |
|---|---:|---:|---:|---:|---:|---:|---:|
| deterministic_oracle | 500 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.0000 | 0 |

The positive class for precision/recall/F1 is ARTIFACT_DETECTED.
Artifact acceptance is the fraction of true artifact cases accepted as VERIFIED.
Negative-result acceptance as a hit is reported separately in JSON.
