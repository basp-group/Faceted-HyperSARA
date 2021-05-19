# Faceting experiments

## Preliminary tests

Image size: 512 x 1024 x 20, 2x2 facets for FHS, no homotopy, 5k max. pdfb iterations

| Run | Algo        | Reg. type | Reg. update | $\alpha$ | $\bar{\alpha}$ | $\upsilon$ | $\bar{\upsilon}$               | $\mu$                     | $\bar{\mu}$               | aSNR  | id              |
| --- | ----------- | --------- | ----------- | -------- | -------------- | ---------- | ------------------------------ | ------------------------- | ------------------------- | ----- | --------------- |
| [ ] | SARA (ch 1) |           |             | 1e-1     | -              |            | -                              |                           | -                         |       |                 |
| [ ] | SARA (ch 1) |           |             | 1        | -              |            | -                              |                           | -                         |       |                 |
| [ ] | SARA (ch 1) |           |             | 10       | -              |            | -                              |                           | -                         |       |                 |
| [ ] | HyperSARA   |           |             | 1e-1     | 1e-1           |            |                                |                           |                           |       |                 |
| [ ] | HyperSARA   |           |             | 1e-1     | 1              |            |                                |                           |                           |       |                 |
| [ ] | HyperSARA   |           |             | 1e-1     | 10             |            |                                |                           |                           |       |                 |
| [ ] | HyperSARA   |           |             | 1        | 1e-1           |            |                                |                           |                           |       |                 |
| [ ] | HyperSARA   |           |             | 1        | 1              |            |                                |                           |                           |       |                 |
| [ ] | HyperSARA   |           |             | 1        | 10             |            |                                |                           |                           |       |                 |
| [ ] | HyperSARA   |           |             | 10       | 1e-1           |            |                                |                           |                           |       |                 |
| [ ] | HyperSARA   |           |             | 10       | 1              |            |                                |                           |                           |       |                 |
| [ ] | HyperSARA   |           |             | 10       | 10             |            |                                |                           |                           |       |                 |
| [x] | FHS (2x2)   | inv       | 0           | 1e-1     | 1e-1           | 2.9079e-07 | [5.8051e-07, 1.3771e-06]       | 1.2522e-04 / -            | 1.0555e-02 /  -           | 29.37 | 795008          |
| [x] | FHS (2x2)   | inv       | 0           | 1e-1     | 1              | "          |                                | "                         | "                         | 29.45 | 795009          |
| [x] | FHS (2x2)   | inv       | 0           | 1e-1     | 10             | "          |                                | "                         | "                         | 29.90 | 795010 / 796646 |
| [x] | FHS (2x2)   | inv       | 0           | 1        | 1e-1           | "          |                                | "                         | "                         | 31.33 | 79502 / 796647  |
| [x] | FHS (2x2)   | inv       | 0           | 1        | 1              | "          |                                | "                         | "                         | 31.14 | 794488          |
| [x] | FHS (2x2)   | inv       | 0           | 1        | 10             | "          |                                | "                         | "                         | 31.14 | 795021 / 796650 |
| [x] | FHS (2x2)   | inv       | 0           | 10       | 1e-1           | "          |                                | "                         | "                         | 31.30 | 795011          |
| [x] | FHS (2x2)   | inv       | 0           | 10       | 1              | "          |                                | "                         | "                         | 31.34 | 795134          |
| [x] | FHS (2x2)   | inv       | 0           | 10       | 10             | "          |                                | "                         | "                         | 32.79 | 795013          |
| [x] | FHS (2x2)   | inv       | 1           | 1e-1     | 1e-1           | 2.9079e-07 | [5.8051e-07, 1.3771e-06] / ... | 1.2522e-04 / 2.084401e-04 | 1.0555e-02 / 1.187855e-04 | 29.37 | 795014          |
| [x] | FHS (2x2)   | inv       | 1           | 1e-1     | 1              | "          | " /                            | " /2.101448e-04           | " / 1.204080e-02          | 29.45 | 795015          |
| [x] | FHS (2x2)   | inv       | 1           | 1e-1     | 10             | "          | " /                            | " / 2.112602e-04          | " / 1.211665e-02          | 29.96 | 795016          |
| [x] | FHS (2x2)   | inv       | 1           | 1        | 1e-1           | "          | " /                            | " / 2.138199e-04          | " / 1.196869e-02          | 31.37 | 795022          |
| [x] | FHS (2x2)   | inv       | 1           | 1        | 1              | "          | " /                            | " / 2.140677e-04          | " / 1.204870e-02          | 31.19 | 794478          |
| [x] | FHS (2x2)   | inv       | 1           | 1        | 10             | "          | " /                            | " / 2.147869e-04          | " / 1.211759e-02          | 31.37 | 795136          |
| [x] | FHS (2x2)   | inv       | 1           | 10       | 1e-1           | "          | " /                            | " / 2.186951e-04          | " / 1.202126e-02          | 31.64 | 795017          |
| [x] | FHS (2x2)   | inv       | 1           | 10       | 1              | "          | " /                            | " / 2.187030e-04          | " / 1.206262e-02          | 31.69 | 795018          |
| [x] | FHS (2x2)   | inv       | 1           | 10       | 10             | "          | " /                            | " / 2.182731e-04          | " / 1.210817e-02          | 33.07 | 795019          |

-> SNR slightly better after update
-> rel. variation starts to oscillate for (10, 10)

Test again with old values: 806976, running in ~/Faceted... (!! BEWARE: need to modify script before launching anything new !!)

---

### Order of magnitude regularization parameters

"Per facet" computation for the low-rank term
| Reg. type | Image approx. | Noise transfer | $\upsilon$ | $\bar{\upsilon}$         | $\mu$      | $\bar{\mu}$ | id     |
| --------- | ------------- | -------------- | ---------- | ------------------------ | ---------- | ----------- | ------ |
| log       | none          | none           | 2.9079e-07 | [1.4224e-06, 2.0568e-06] | 3.3464e-01 | 9.7886e+02  | 796125 |
| log       | none          | precond        | 6.5492e-05 | [8.0022e-05, 1.3153e-04] | 5.9271e-03 | 3.0180e+01  | 796126 |
| log       | precond       | none           | 2.9079e-07 | [5.8051e-07, 1.3771e-06] | 1.3371e-01 | 1.0278e+03  | 796130 |
| log       | precond       | precond        | 6.5492e-05 | [9.2640e-05, 1.7364e-04] | 2.2478e-03 | 1.2166e+01  | 796129 |

Single facet
| Reg. type | Image approx. | Noise transfer | $\upsilon$ | $\bar{\upsilon}$ | $\mu$      | $\bar{\mu}$ | id     |
| --------- | ------------- | -------------- | ---------- | ---------------- | ---------- | ----------- | ------ |
| log       | none          | none           |            |                  |            |             | 806961 |
| log       | none          | precond        | 6.5492e-05 | 1.0377e-04       | 5.9271e-03 | 1.0121e+02  | 798317 |
| log       | precond       | none           | 2.9079e-07 | 1.4449e-06       | 1.3371e-01 | 2.4871e+03  | 798320 |
| log       | precond       | precond        | 6.5492e-05 | 1.2300e-04       | 2.2478e-03 | 4.2927e+01  | 798319 |
| inv       | none          | none           | 2.9079e-07 | 1.4449e-06       | 3.3377e-04 | 9.3655e-02  | 806948 |
| inv       | none          | precond        | 6.5492e-05 | 1.2300e-04       | 3.3377e-04 | 9.3655e-02  | 806949 |
| inv       | precond       | none           | 2.9079e-07 | 1.4449e-06       | 1.2522e-04 | 1.4981e-02  | 806952 |
| inv       | precond       | precond        | 6.5492e-05 | 1.2300e-04       | 1.2522e-04 | 1.4981e-02  | 806950 |

---

## Spatial faceting

Image size: 1024 x 2048 x 20.

"Per facet" computation for the low-rank term (Qx=Qy=2, 50% overlap.)
| Reg. type     | Image approx. | Noise transfer | $\upsilon$ | $\bar{\upsilon}$         | $\mu$      | $\bar{\mu}$ | id     |
| ------------- | ------------- | -------------- | ---------- | ------------------------ | ---------- | ----------- | ------ |
| log           | none          | none           | 6.0010e-07 | [2.0873e-06, 4.5724e-06] | 5.4684e-02 | 4.5052e+02  | 796652 |
| log           | none          | precond        | 9.6485e-05 | [1.5802e-04, 1.7333e-04] | 1.2692e-03 | 1.5704e+01  | 796653 |
| log           | precond       | none           | 6.0010e-07 | [2.1773e-06, 3.9061e-06] | 2.2885e-02 | 3.3885e+02  | 796651 |
| log           | precond       | precond        | 9.6485e-05 | [1.0510e-04, 1.8453e-04] | 5.8640e-04 | 9.9955e+00  | 796639 |

Single facet
| Reg. type | Image approx. | Noise transfer | $\upsilon$ | $\bar{\upsilon}$ | $\mu$      | $\bar{\mu}$ | id     |
| --------- | ------------- | -------------- | ---------- | ---------------- | ---------- | ----------- | ------ |
| log       | none          | none           | 6.0010e-07 | 5.3544e-06       | 5.4684e-02 | 1.0058e+03  | 798285 |
| log       | none          | precond        | 9.6485e-05 | 2.0485e-04       | 1.2692e-03 | 4.3194e+01  | 798286 |
| log       | precond       | none           | 6.0010e-07 | 3.9700e-06       | 2.2885e-02 | 9.2375e+02  | 798284 |
| log       | precond       | precond        | 9.6485e-05 | 1.2259e-04       | 5.8640e-04 | 3.9970e+01  | 798282 |
| inv       | none          | none           | 6.0010e-07 | 5.3544e-06       |            |             | 806954 |
| inv       | none          | precond        | 9.6485e-05 | 2.0485e-04       |            |             | 806955 |
| inv       | precond       | none           | 6.0010e-07 | 3.9700e-06       |            |             | 806956 |
| inv       | precond       | precond        | 9.6485e-05 | 1.2259e-04       |            |             | 806958 |

Test with svd of the noise (ot of the dirty image) to compute upsilon_bar (log option) -> 807235 (ongoing, test/...)

Test (image, noise_transfer) = (precond, precond), reg. parameters computed from a single facet (ongoing, test/...)

| Run | Algo      | Reg. type | Reg. update | $\alpha$ | $\bar{\alpha}$ | $\upsilon$ | $\bar{\upsilon}$ | $\mu$ | $\bar{\mu}$ | aSNR | id     |
| --- | --------- | --------- | ----------- | -------- | -------------- | ---------- | ---------------- | ----- | ----------- | ---- | ------ |
| [ ] | FHS (2x2) | log       | 0           | 1e-1     | 1e-1           |            |                  |       |             |      | 807239 |
| [ ] | FHS (2x2) | log       | 0           | 1e-1     | 1              |            |                  |       |             |      | 807240 |
| [ ] | FHS (2x2) | log       | 0           | 1        | 1e-1           |            |                  |       |             |      | 807241 |
| [ ] | FHS (2x2) | log       | 0           | 1        | 1              |            |                  |       |             |      | 807242 |

---

## Spectral faceting

Image size: 256 x 512 x 100.

| Reg. type | Image approx. | Noise transfer | $\upsilon$ | $\bar{\upsilon}$ | $\mu$ | $\bar{\mu}$ | id  |
| --------- | ------------- | -------------- | ---------- | ---------------- | ----- | ----------- | --- |
| log       | none          | none           |            |                  |       |             |     |
| log       | none          | precond        |            |                  |       |             |     |
| log       | precond       | none           |            |                  |       |             |     |
| log       | precond       | precond        |            |                  |       |             |     |
