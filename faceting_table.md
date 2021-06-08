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
| [x] | FHS (2x2)   | log (old) | 1           | 10       | 10             | 1.2574e-03 | [6.1506e-03, 8.8939e-03]       | 7.7386e-05                | 2.2637e-01                | 31.12 | 806976          |

old values, with stopping criterion raised to 5e-5, pdfb stops really early! (1500 iterations in total over the 30 reweights, which is a problem) -> relaunch with the stopping criterion back to 1e-5? (and w/o min. number of iterations?)

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
| log       | none          | none           | 2.9079e-07 | 3.0812e-06       | 3.3464e-01 | 1.9673e+03  | 809077 |
| log       | none          | precond        | 6.5492e-05 | 1.0377e-04       | 5.9271e-03 | 1.0121e+02  | 798317 |
| log       | precond       | none           | 2.9079e-07 | 1.4449e-06       | 1.3371e-01 | 2.4871e+03  | 798320 |
| log       | precond       | precond        | 6.5492e-05 | 1.2300e-04       | 2.2478e-03 | 4.2927e+01  | 798319 |
| inv       | none          | none           | 2.9079e-07 | 1.4449e-06       | 3.3377e-04 | 9.3655e-02  | 806948 |
| inv       | none          | precond        | 6.5492e-05 | 1.2300e-04       | 3.3377e-04 | 9.3655e-02  | 806949 |
| inv       | precond       | none           | 2.9079e-07 | 1.4449e-06       | 1.2522e-04 | 1.4981e-02  | 806952 |
| inv       | precond       | precond        | 6.5492e-05 | 1.2300e-04       | 1.2522e-04 | 1.4981e-02  | 806950 |

Final runs (spatial faceting, inv, (none, none), SVD of the dirty image to compute $\bar{\upsilon}$)
| Run | Algo      | Reg. type        | Reg. update | $\alpha$ | $\bar{\alpha}$ | $\upsilon$ | $\bar{\upsilon}$ | $\mu$ | $\bar{\mu}$ | aSNR | id  |
| --- | --------- | ---------------- | ----------- | -------- | -------------- | ---------- | ---------------- | ----- | ----------- | ---- | --- |
| [R] | SARA      | inv (none, none) | 0           | 1        | -              |            |                  |       |             |      |     |
| [R] | HS        | inv (none, none) | 0           | 1        | 1              |            |                  |       |             |      |     |
| [R] | FHS (4x4) | inv (none, none) | 0           | 1        | 1              |            |                  |       |             |      |     |

---

## Spatial faceting

Image size: 1024 x 2048 x 20.

"Per facet" computation for the low-rank term (Qx=Qy=2, 50% overlap.)
| Reg. type | Image approx. | Noise transfer | $\upsilon$ | $\bar{\upsilon}$         | $\mu$      | $\bar{\mu}$ | id     |
| --------- | ------------- | -------------- | ---------- | ------------------------ | ---------- | ----------- | ------ |
| log       | none          | none           | 6.0010e-07 | [2.0873e-06, 4.5724e-06] | 5.4684e-02 | 4.5052e+02  | 796652 |
| log       | none          | precond        | 9.6485e-05 | [1.5802e-04, 1.7333e-04] | 1.2692e-03 | 1.5704e+01  | 796653 |
| log       | precond       | none           | 6.0010e-07 | [2.1773e-06, 3.9061e-06] | 2.2885e-02 | 3.3885e+02  | 796651 |
| log       | precond       | precond        | 9.6485e-05 | [1.0510e-04, 1.8453e-04] | 5.8640e-04 | 9.9955e+00  | 796639 |

Single facet
| Reg. type          | Image approx. | Noise transfer | $\upsilon$ | $\bar{\upsilon}$ | $\mu$      | $\bar{\mu}$ | id     |
| ------------------ | ------------- | -------------- | ---------- | ---------------- | ---------- | ----------- | ------ |
| log                | none          | none           | 6.0010e-07 | 5.3544e-06       | 5.4684e-02 | 1.0058e+03  | 798285 |
| log                | none          | precond        | 9.6485e-05 | 2.0485e-04       | 1.2692e-03 | 4.3194e+01  | 798286 |
| log                | precond       | none           | 6.0010e-07 | 3.9700e-06       | 2.2885e-02 | 9.2375e+02  | 798284 |
| log                | precond       | precond        | 9.6485e-05 | 1.2259e-04       | 5.8640e-04 | 3.9970e+01  | 798282 |
| inv                | none          | none           | 6.0010e-07 | 5.3544e-06       | 8.3710e-05 | 3.7578e-02  | 808666 | 4.054766e+01 -> run all solvers in this configuration |
| inv                | none          | precond        | 9.6485e-05 | 2.0485e-04       | 8.3710e-05 | 3.7578e-02  | 808667 | nothing here                                          |
| inv                | precond       | none           | 6.0010e-07 | 3.9700e-06       | 3.9072e-05 | 7.4795e-03  | 808669 | 3.615950e+01                                          |
| inv                | precond       | precond        | 9.6485e-05 | 1.2259e-04,      | 3.9072e-05 | 7.4795e-03  | 808668 | 3.615950e+01                                          |
| log (svd of noise) | precond       | precond        | 9.6485e-05 | 2.0693e-02       | 5.8640e-04 | 1.3523e+00  | 807235 |

**Remark**: the last option might be more appropriate in term of values (+ a priori no need to play with $\alpha$, $\bar{\alpha}$)

Test (image, noise_transfer) = (precond, precond), reg. parameters computed from a single facet (ongoing, test/...)
Results after 700 pdfb iterations (not even 1st reweight yet)

| Run | Algo      | Reg. type | Reg. update | $\alpha$ | $\bar{\alpha}$ | $\upsilon$ | $\bar{\upsilon}$ | $\mu$      | $\bar{\mu}$ | aSNR         | id     |
| --- | --------- | --------- | ----------- | -------- | -------------- | ---------- | ---------------- | ---------- | ----------- | ------------ | ------ |
| [R] | FHS (2x2) | log       | 0           | 1e-1     | 1e-1           | 9.6485e-05 | 1.2259e-04       | 5.8640e-04 | 3.9970e+01  | 3.524659e+01 | 808670 |
| [R] | FHS (2x2) | log       | 0           | 1e-1     | 1              | "          | "                | "          | "           | 3.127044e+01 | 808671 |
| [R] | FHS (2x2) | log       | 0           | 1        | 1e-1           | "          | "                | "          | "           | 3.535832e+01 | 808672 |
| [R] | FHS (2x2) | log       | 0           | 1        | 1              | 9.6485e-05 | 1.2259e-04       | 5.8640e-04 | 3.9970e+01  | 3.132294e+01 | 808673 |

Final runs (spatial faceting, inv, (none, none), SVD of the dirty image to compute $\bar{\upsilon}$)
| Run | Algo      | Reg. type        | Reg. update | $\alpha$ | $\bar{\alpha}$ | $\upsilon$ | $\bar{\upsilon}$ | $\mu$      | $\bar{\mu}$ | aSNR  | id              |
| --- | --------- | ---------------- | ----------- | -------- | -------------- | ---------- | ---------------- | ---------- | ----------- | ----- | --------------- |
| [x] | SARA      | inv (none, none) | 0           | 1        | -              | ...        | ...              | ...        | ...         | 35.28 | 810790 - 810809 |
| [x] | HS        | inv (none, none) | 0           | 1        | 1              | 6.0010e-07 | 5.3544e-06       | 8.3710e-05 | 3.7578e-02  | 34.90 | 810313          |
| [x] | FHS (4x4) | inv (none, none) | 0           | 1        | 1              | 6.0010e-07 | 5.3544e-06       | 8.3710e-05 | 3.7578e-02  | 34.90 | 810314          |

-> reweighting sort of kill the reconstruction quality in this case (40 dB before first reweight, then decreases significantly...)
-> same behaviour actually observed with SARA
-> remove $\upsilon$ from the numerator when defining the weights? actually, $\upsilon$ and $\bar{\upsilon}$ are too small.

- After latest fix (add heuristic, normaliza by operator norm per channel, ...)
| Run | Algo      | Rw. type  | Reg. type | Noise transfer | Xdirty  | reg. option (SVD of noise, or Xdirty) | $\upsilon$ | $\bar{\upsilon}$         | $\mu$      | $\bar{\mu}$ | aSNR  | id     |
| --- | --------- | --------- | --------- | -------------- | ------- | ------------------------------------- | ---------- | ------------------------ | ---------- | ----------- | ----- | ------ |
| [x] | HS        | dirty     | inv       | psf            | none    | none                                  | 7.8432e-03 | 1.9261e+00               | 7.7626e-05 | 3.4569e-02  | 42.29 | 847555 |
| [x] | HS        | dirty     | log       | precond        | precond | none                                  | 1.0025e-04 | 2.0726e-02               | 5.6404e-04 | 4.7098e-01  | 42.97 | 847556 |
| [x] | HS        | heuristic | inv       | none           | none    | none                                  | 1.0965e-04 | 1.0652e-01               | 7.7626e-05 | 3.4569e-02  | 40.92 | 847557 |
| [R] | HS        | heuristic | log       | none           | none    | none                                  |            |                          |            |             |       | 877816 |
| [R] | FHS (4x4) | dirty     | inv       | none           | none    | none                                  |            |                          |            |             |       |        | -> solve  |
| [R] | FHS (4x4) | dirty     | inv       | none           | none    | dirty                                 |            |                          |            |             |       |        | -> solve  |
| [x] | FHS (4x4) | dirty     | inv       | psf            | none    | none                                  | 7.8432e-03 | 1.9261e+00               | 7.7626e-05 | 3.4569e-02  | 42.25 | 848745 |
| [K] | FHS (4x4) | dirty     | log       | precond        | precond | none                                  | 1.0025e-04 | 2.0726e-02               | 5.6404e-04 | 4.7098e-01  | 43.51 | 848746 | iter 5200 |
| [x] | FHS (4x4) | heuristic | inv       | none           | none    | none                                  | 1.0965e-04 | [8.8810e-03, 3.5512e-02] | 7.7626e-05 | 3.4569e-02  | 41.62 | 848748 |
| [K] | FHS (4x4) | heuristic | log       | none           | none    | none                                  | 1.0965e-04 | [8.8810e-03, 3.5512e-02] | 1.1355e-03 | 2.0797e-01  | 43.02 | 848749 | iter 5200 |

-> for inv: could move stopping criterion to 1e-5 for pdfb (would allow to run a bit longer)

Using the heuristic for both $\mu = \upsilon}$, $\bar{\mu} = \bar{\upsilon}$ (stoppign crit: 5e-5) (in ~/test/Faceted...)
| Run | Algo      | Rw. type  | Reg. type | Noise transfer | Xdirty | reg. option | $\upsilon$               | $\bar{\upsilon}$         | $\mu$                    | $\bar{\mu}$              | aSNR  | id              |
| --- | --------- | --------- | --------- | -------------- | ------ | ----------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ----- | --------------- |
| [x] | SARA      | heuristic | heuristic | none           | none   | none        | [1.7789e-05, 3.3730e-05] | -                        | [1.7789e-05, 3.3730e-05] | -                        | 35.87 | 878221 - 878239 |
| [x] | HS        | heuristic | heuristic | none           | none   | none        | 1.0868e-04               | 1.0652e-01               | 1.0868e-04               | 1.0652e-01               | 41.87 | 878185          |
| [x] | FHS (4x4) | heuristic | heuristic | none           | none   | none        | 1.0868e-04               | [8.8810e-03, 3.5512e-02] | 1.0868e-04               | [8.8810e-03, 3.5512e-02] | 41.72 | 878184          |

Using the heuristic for both $\mu = \upsilon}$, $\bar{\mu} = \bar{\upsilon}$ (stoppign crit: 1e-5) (in ~/test/Faceted...)
| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$               | $\bar{\upsilon}$         | $\mu$                    | $\bar{\mu}$              | aSNR         | id                |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------ | ----------------- |
| [x] | SARA (a=1)     | heuristic  | heuristic  | none         | none   | none        | [1.7789e-05, 3.3730e-05] | -                        | [1.7789e-05, 3.3730e-05] | -                        | 36.78        | 893579 - 893607   |
| [x] | SARA (a=2)     | heuristic  | heuristic  | none         | none   | none        | [1.7789e-05, 3.3730e-05] | -                        | [1.7789e-05, 3.3730e-05] | -                        | 39.98        |                   |
| [R] | SARA (a=3)     | heuristic  | heuristic  | none         | none   | none        | [1.7789e-05, 3.3730e-05] | -                        | [1.7789e-05, 3.3730e-05] | -                        | 41.33        | 1110145 - 1110164 |
| [x] | HS (a,ab=1)    | heuristic  | heuristic  | none         | none   | none        | 1.0868e-04               | 1.0652e-01               | 1.0868e-04               | 1.0652e-01               | 42.80        |                   |
| [x] | FHS (a,ab=1)   | heuristic  | heuristic  | none         | none   | none        | 1.0868e-04               | [8.8810e-03, 3.5512e-02] | 1.0868e-04               | [8.8810e-03, 3.5512e-02] | 44.07        | 893342            |
| [x] | HS  (a,ab=3)   | heuristic  | heuristic  | none         | none   | none        | 1.0868e-04               | 1.0652e-01               | 1.0868e-04               | 1.0652e-01               | 43.11        | 1306789           |
| [K] | FHS (a,ab=3)   | heuristic  | heuristic  | none         | none   | none        | 1.0868e-04               | [8.8810e-03, 3.5512e-02] | 1.0868e-04               | [8.8810e-03, 3.5512e-02] | 43.12        | 1120950           |
| [x] | HS  (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 1.5997e-04               | 1.0652e-01               | 1.5997e-04               | [1.0652e-01, 1.0652e-01] | 42.94        | 1308520           |
| [K] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 1.5997e-04               | [8.8810e-03, 3.5512e-02] | 1.5997e-04               | [8.8810e-03, 3.5512e-02] | 44.01 (rw 5) | 1308860           |
| [x] | FHS (a=1,ab=1) | heuristic2 | heuristic2 | none         | none   | none        | 1.5997e-04               | [8.8810e-03, 3.5512e-02] | 1.5997e-04               | [8.8810e-03, 3.5512e-02] | 44.17        | 1312320           |

Average "target" SNR: 47.93 dB (alpha = 1), 41.91 dB (alpha = 2), 38.40 dB (alpha = 3)

Setting 500 iterations max per pdfb (in ~/test/Faceted...)
| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$               | $\bar{\upsilon}$         | $\mu$                    | $\bar{\mu}$              | aSNR  | id                 |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ----- | ------------------ |
| [x] | SARA (a=3)     | heuristic  | heuristic  | none         | none   | none        | [1.7789e-05, 3.3730e-05] | -                        | [1.7789e-05, 3.3730e-05] | -                        | 41.33 | 1313495 - 1313514  |
| [x] | HS (a=1,ab=3)  | heuristic2 | heuristic2 | none         | none   | none        | 1.5997e-04               | 1.0652e-01               | 1.5997e-04               | 1.0652e-01               | 42.94 | 1313494            |
| [R] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 1.5997e-04               | [8.8810e-03, 3.5512e-02] | 1.5997e-04               | [8.8810e-03, 3.5512e-02] | 44.15 | 1313824 (in rw 19) |
-> restart from rw 15 (1319645)

New data set with a differet noise realization (different folder, seed = 54321), 2k pdfb iterations (in ~/Faceted...)
data generation: job 1312342
| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$ | $\bar{\upsilon}$         | $\mu$      | $\bar{\mu}$              | aSNR  | id                                                  |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ---------- | ------------------------ | ---------- | ------------------------ | ----- | --------------------------------------------------- |
| [E] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 1.5997e-04 | [8.8810e-03, 3.5512e-02] | 1.5997e-04 | [8.8810e-03, 3.5512e-02] | 43.90 | 1313822 (in rw 6, already more than 10k iterations) | HDF5 file error (??) | -> no need to restart or retrieve (same behaviour as before) |

### New data set A+C configuration (test_configAC, not much difference to be expected in terms of quality) (500 iter per pdfb, 30 reweights)

| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$               | $\bar{\upsilon}$         | $\mu$                    | $\bar{\mu}$              | aSNR  | id                 |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ----- | ------------------ |
| [K] | SARA (a=3)     | heuristic  | heuristic  | none         | none   | none        | [2.2002e-05, 4.0750e-05] | -                        | [2.2002e-05, 4.0750e-05] | -                        | 41.01 | 1318405  - 1318424 | -> channel 11 missing  1319652     |
| [K] | HS (a=1,ab=3)  | heuristic2 | heuristic2 | none         | none   | none        | 1.9527e-04               | 1.3003e-01               | 1.9527e-04               | 1.3003e-01               | 42.14 | 1317496, 1319676   | error in reweighting 1, hdf5 error |
| [K] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 1.9527e-04               | [1.0841e-02, 4.3347e-02] | 1.9527e-04               | [1.0841e-02, 4.3347e-02] | 44.30 | 1318404            | still running     (rw 5)           |

### A configuration (test_configA) (500 iter per pdfb, 30 reweights)

| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$               | $\bar{\upsilon}$         | $\mu$                    | $\bar{\mu}$              | aSNR  | id                |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ----- | ----------------- |
| [x] | SARA (a=3)     | heuristic  | heuristic  | none         | none   | none        | [2.8344e-05, 4.0439e-05] | -                        | [2.8344e-05, 4.0439e-05] | -                        | 37.25 | 1317262 - 1317281 | -> retrieved files                |
| [K] | HS (a=1,ab=3)  | heuristic2 | heuristic2 | none         | none   | none        | 2.2254e-04               | 1.4819e-01               | 2.2254e-04               | 1.4819e-01               | 24.67 | 1314934           | -> error in rw 2, restart 1319650 | -> investigate |
| [K] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 2.2254e-04               | [1.2355e-02, 4.9401e-02] | 2.2254e-04               | [1.2355e-02, 4.9401e-02] | 41.00 | 1315277           | -> error in rw 9, restart 1319649 |

2k iterations, 5 reweights
| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$ | $\bar{\upsilon}$         | $\mu$      | $\bar{\mu}$              | aSNR  | id                |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ---------- | ------------------------ | ---------- | ------------------------ | ----- | ----------------- |
| [x] | SARA (a=3)     | heuristic  | heuristic  | none         | none   | none        |            | -                        |            | -                        |       | 1320764 - 1320783 |
| [R] | HS (a=1,ab=3)  | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          |       | 1323195           |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 2.2254e-04 | [1.2355e-02, 4.9401e-02] | 2.2254e-04 | [1.2355e-02, 4.9401e-02] | 40.75 | 1320763, 1322375  | 4x4, restart rw2 |
| [R] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          |       | 1321020, 1322373  | 2x2, restart rw0 | stagnates around 33 dB...? (issue here?)
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 40.75 | 1321021, 1322374  | 3x3, restart rw1 |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 38.39 | 1323194           | 0                |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 41.01 | 1322378           | 0.1              |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 41.47 | 1322379           | 0.25             |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 41.22 | 1322380           | 0.4              |




### Control simulations (test/..., config A)

2k iterations, 5 reweights
| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$ | $\bar{\upsilon}$ | $\mu$ | $\bar{\mu}$ | aSNR | id      |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ---------- | ---------------- | ----- | ----------- | ---- | ------- |
| [R] | HS (a=1,ab=3)  | heuristic2 | heuristic2 | none         | none   | none        |            |                  |       |             |      | 1323437 |
| [R] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                  |       |             |      | 1323438 |
