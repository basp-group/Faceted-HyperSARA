# Faceting experiments

## Preliminary tests (ABCD uv configurations)

Image size: 512 x 1024 x 20, 2x2 facets for FHS, no homotopy, 5k max. pdfb iterations

| Run | Algo        | Reg. type | Reg. update | $\alpha$ | $\bar{\alpha}$ | $\upsilon$ | $\bar{\upsilon}$               | $\mu$                     | $\bar{\mu}$               | aSNR  | id              |
| --- | ----------- | --------- | ----------- | -------- | -------------- | ---------- | ------------------------------ | ------------------------- | ------------------------- | ----- | --------------- |
| [ ] | SARA (ch 1) |           |             | 1e-1     | -              |            | -                              |                           | -                         |       |                 |
| [ ] | HyperSARA   |           |             | 1e-1     | 1e-1           |            |                                |                           |                           |       |                 |
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

### Spatial faceting

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
| inv                | none          | none           | 6.0010e-07 | 5.3544e-06       | 8.3710e-05 | 3.7578e-02  | 808666 |
| inv                | none          | precond        | 9.6485e-05 | 2.0485e-04       | 8.3710e-05 | 3.7578e-02  | 808667 |
| inv                | precond       | none           | 6.0010e-07 | 3.9700e-06       | 3.9072e-05 | 7.4795e-03  | 808669 |
| inv                | precond       | precond        | 9.6485e-05 | 1.2259e-04,      | 3.9072e-05 | 7.4795e-03  | 808668 |
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

## Current simulations: A configuration (test_configA)

### 500 iter per pdfb, 30 reweights

| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$               | $\bar{\upsilon}$         | $\mu$                    | $\bar{\mu}$              | aSNR  | id                |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ | ----- | ----------------- |
| [x] | SARA (a=3)     | heuristic  | heuristic  | none         | none   | none        | [2.8344e-05, 4.0439e-05] | -                        | [2.8344e-05, 4.0439e-05] | -                        | 37.25 | 1317262 - 1317281 | -> retrieved files                |
| [K] | HS (a=1,ab=3)  | heuristic2 | heuristic2 | none         | none   | none        | 2.2254e-04               | 1.4819e-01               | 2.2254e-04               | 1.4819e-01               | 24.67 | 1314934           | -> error in rw 2, restart 1319650 |
| [K] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 2.2254e-04               | [1.2355e-02, 4.9401e-02] | 2.2254e-04               | [1.2355e-02, 4.9401e-02] | 41.00 | 1315277           | -> error in rw 9, restart 1319649 |

### 2k iterations, 5 reweights

| Run | Algo           | Rw.        | Reg.       | Noise trans. | Xdirty | reg. option | $\upsilon$ | $\bar{\upsilon}$         | $\mu$      | $\bar{\mu}$              | aSNR  | id                |
| --- | -------------- | ---------- | ---------- | ------------ | ------ | ----------- | ---------- | ------------------------ | ---------- | ------------------------ | ----- | ----------------- |
| [x] | SARA (a=3)     | heuristic  | heuristic  | none         | none   | none        |            | -                        |            | -                        |       | 1320764 - 1320783 |
| [R] | HS (a=1,ab=3)  | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          |       | 1323195           |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        | 2.2254e-04 | [1.2355e-02, 4.9401e-02] | 2.2254e-04 | [1.2355e-02, 4.9401e-02] | 40.75 | 1320763, 1322375  | 4x4, restart rw2 |
| [R] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          |       | 1321020, 1322373  | 2x2, restart rw0 | stagnates around 33 dB...? (issue here?) |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 40.75 | 1321021, 1322374  | 3x3, restart rw1 |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 38.39 | 1323194           | 0                |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 41.01 | 1322378           | 0.1              |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 41.47 | 1322379           | 0.25             |
| [x] | FHS (a=1,ab=3) | heuristic2 | heuristic2 | none         | none   | none        |            |                          |            |                          | 41.22 | 1322380           | 0.4              |

| Run | Algo            | $\mu$ | rw0   | rw 1  | rw 2  | rw 3  | rw 4  | rw 5             | asnr_log (x10, /10)        | arun             | acpu              | total run | total_cpu | id                                                |
| --- | --------------- | ----- | ----- | ----- | ----- | ----- | ----- | ---------------- | -------------------------- | ---------------- | ----------------- | --------- | --------- | ------------------------------------------------- |
| [x] | SARA            |       | -     | -     | -     | 35.10 | -     | 36.28 (6.77e-01) | 16.90                      | 3.34 (1.53e-01)  | 7.25 (4.28e-01)   |           |           | 1320764-1320783, 1330682-1330701, 1330682-1330701 | ok |
| [E] | HS              |       | 22.37 | 30.62 | 32.16 |       |       |                  |                            |                  |                   |           |           | 1323195, 1332672,                                 |
| [K] | HS (restart)    |       |       |       |       |       |       |                  |                            |                  |                   |           |           | 1333864                                           |
| [K] | FHS (2x2, 0.25) |       | 24.37 | 36.11 | 37.01 | 36.65 | 36.29 |                  |                            |                  |                   |           |           | 1327578, 1330711, 1332233                         |
| [x] | FHS (3x3 0.25)  |       | 24.71 | 37.27 | 39.60 | 39.99 | 40.05 | 40.05            |                            |                  |                   |           |           | ... 1332232                                       |
| [x] | FHS (2x2, 0.5)  |       | 23.78 | 33.07 | 33.87 | 32.42 | 31.30 | 31.30 (3.34e+00) | 15.69 (2.03e+01, 1.09e+01) | 31.03 (3.84e+00) | 212.84 (2.56e+01) |           |           | 1321020, 1322373                                  |
| [x] | FHS (3x3, 0.5)  |       | 24.72 | 36.27 | 39.83 | 40.50 | 40.69 | 40.75 (2.87e+00) | 21.57 (2.78e+01, 1.64e+01) | 17.20 (7.75e-01) | 237.49 (9.46e+00) |           |           | 1321021, 1322374                                  |

### 5k pdfb iterations

| Run | Algo            | $\mu$ | rw0   | rw 1  | rw 2  | rw 3  | rw 4  | rw 5             | asnr_log (x10, /10)        | arun             | acpu              |
| --- | --------------- | ----- | ----- | ----- | ----- | ----- | ----- | ---------------- | -------------------------- | ---------------- | ----------------- |
| [x] | FHS (4x4, 0)    |       | 24.31 | 36.80 | 38.17 | 38.37 | 38.44 | 38.39 (1.59e+00) | 21.41 (2.76e+01, 1.55e+01) | 11.47 (5.66e-01) | 281.29 (1.29e+01) |
| [x] | FHS (4x4, 0.1)  |       | 24.37 | 38.34 | 40.55 | 40.93 | 41.01 | 41.02 (2.73e+00) | 21.64 (2.81e+01, 1.57e+01) | 11.52 (4.26e-01) | 282.69 (1.01e+01) |
| [x] | FHS (4x4, 0.25) |       | 24.16 | 38.34 | 40.99 | 41.40 | 41.47 | 41.48 (2.55e+00) | 22.34 (2.89e+01, 1.65e+01) | 11.95 (4.63e-01) | 291.61 (1.17e+01) |
| [x] | FHS (4x4, 0.4)  |       | 24.03 | 38.81 | 40.86 | 41.18 | 41.23 | 41.22 (2.28e+00) | 22.31 (2.86e+01, 1.70e+01) | 12.84 (5.37e-01) | 306.42 (1.38e+01) |
| [x] | FHS (4x4, 0.5)  |       | 25.28 | 38.26 | 40.31 | 40.69 | 40.75 | 40.76 (2.75e+00) | 21.37 (2.76e+01, 1.62e+01) | 13.60 (5.28e-01) | 317.53 (1.58e+01) |

| --- | --------------- | ------- | ----- | ----- | ----- |
| [R] | HS              | 1338171 | 22.13 |       |       | 2295, 27.70 (in rw 1, rel_val = 9.953951e-05), mu = 2.2254e-04, mu_bar = [1.4819e-01, 1.4819e-01] |
| [R] | FHS (4x4, 0)    | 1336079 | 24.35 | 36.64 | 38.11 | -> mu = 2.2254e-04, mu_bar = [3.7047e-02, 3.7047e-02]                                             |
| [R] | FHS (4x4, 0.1)  | 1336080 | 24.15 | 38.03 | 40.38 | ->                                                                                                |
| [R] | FHS (4x4, 0.25) | 1336081 | 24.07 | 38.22 | 40.79 | -> mu = 1.7461e-04, mu_bar = [2.8860e-02, 4.1161e-02] to be checked                               |
| [R] | FHS (4x4, 0.4)  | 1336082 | 24.19 | 38.70 | 40.68 | -> mu = 2.2254e-04, mu_bar = [2.0632e-02, 4.5269e-02]                                             |
| [R] | FHS (4x4, 0.5)  | 1336083 | 25.45 | 38.06 | 40.09 | -> mu = 2.2254e-04, mu_bar = [1.2355e-02, 4.9401e-02]                                             |
| [R] | FHS (2x2, 0.25) | 1336084 | 24.37 | 35.38 | 35.43 | -> mu = 2.2254e-04, mu_bar = [5.7678e-02, 8.2316e-02]                                             |
| [R] | FHS (3x3, 0.25) | 1336085 | 24.71 | 36.87 | 39.36 | -> mu = 2.2254e-04, mu_bar = [3.8470e-02, 5.4859e-02]                                             |

### Current runs: 2k iterations, 5 rw, heuristic2

| Run | Algo                                      | rw0   | rw 1  | rw 2  | rw3   | rw4   | rw5   | mu                       | mu_bar                   | mu_chi, sigma_chi      | id              |
| --- | ----------------------------------------- | ----- | ----- | ----- | ----- | ----- | ----- | ------------------------ | ------------------------ | ---------------------- | --------------- |
| [x] | SARA (1xsigma)                            |       |       |       | 28.21 |       | 30.25 | [2.8344e-05, 4.0439e-05] | -                        | -                      | 1341729-1341748 |
| [x] | SARA (3xsigma)                            |       |       |       | 35.10 |       | 36.28 |                          | -                        | -                      | 1345149-1345168 |
| [K] | HS (1xsigma, 1xsigma)                     | 23.25 | 28.64 |       |       |       |       | 1.7461e-04               | 1.4819e-01               | 1.5065e-04, 2.3964e-05 | 1341721         |
| [R] | HS (v4, 1xsigma, 1/2 sigma)               | 23.60 | 28.48 | 34.02 | 34.90 | 35.53 |       | 1.7461e-04               | 7.4095e-02               | 1.5065e-04, 2.3964e-05 | 1351657         | test_configA |
| [R] | HS (v4, 1xsigma, 3xsigma bar, heuristic3) | 23.51 | 28.09 | 33.41 | 34.45 |       |       | 1.7461e-04               | 3x3.3136e-02             | 1.5065e-04, 2.3964e-05 | 1351656         | test_configA |
| [x] | FHS (4x4, 0)                              | 24.90 | 35.12 | 40.14 | 40.33 | 40.44 | 40.50 | 1.7461e-04               | [3.7047e-02, 3.7047e-02] | 1.5065e-04, 2.3964e-05 | 1341722         |
| [x] | FHS (4x4, 0.1)                            | 24.70 | 34.80 | 40.27 | 41.36 | 41.60 | 41.65 | 1.7461e-04               | [3.4363e-02, 3.8415e-02] | 1.5065e-04, 2.3964e-05 | 1341723         |
| [x] | FHS (4x4, 0.25)                           | 24.45 | 33.70 | 40.08 | 41.33 | 41.57 | 41.58 | 1.7461e-04               | [2.8860e-02, 4.1161e-02] | 1.5065e-04, 2.3964e-05 | 1341724         |
| [x] | FHS (4x4, 0.4)                            | 24.49 | 33.36 | 39.84 | 41.07 | 41.13 | 41.16 | 1.7461e-04               | [2.0632e-02, 4.5269e-02] | 1.5065e-04, 2.3964e-05 | 1341725         |
| [x] | FHS (4x4, 0.5)                            | 25.12 | 33.59 | 39.48 | 39.75 | 39.88 | 39.93 | 1.7461e-04               | [1.2355e-02, 4.9401e-02] | 1.5065e-04, 2.3964e-05 | 1341726         |
| [x] | FHS (2x2, 0.25)                           | 24.71 | 32.97 | 39.02 | 39.97 | 40.07 | 40.10 | 1.7461e-04               | [5.7678e-02, 8.2316e-02] | 1.5065e-04, 2.3964e-05 | 1341727         |
| [x] | FHS (3x3, 0.25)                           | 25.08 | 34.29 | 39.97 | 40.45 | 40.51 | 40.54 | 1.7461e-04               | [3.8470e-02, 5.4859e-02] | 1.5065e-04, 2.3964e-05 | 1341728         |

heuristic 3 (test/, 3xsigma everywhere) (probably 3x mu_chi+3xsigma_chi, 1xsigma)
| Run | Algo            | rw0   | rw 1  | rw 2  | rw3 | rw4 | rw5 | mu         | mu_bar                   | mu_chi, sigma_chi      | id      |
| --- | --------------- | ----- | ----- | ----- | --- | --- | --- | ---------- | ------------------------ | ---------------------- | ------- |
| [R] | HS (v4)         | 22.99 | 33.54 | 38.28 |     |     |     | 2.2254e-04 | 3x3.3136e-02             | 1.5065e-04, 2.3964e-05 | 1352436 | test_configA |
| [R] | FHS (4x4, 0.25) |       |       |       |     |     |     | 1.7461e-04 | [2.8860e-02, 4.1161e-02] | 1.5065e-04, 2.3964e-05 | 1362906 | test         |
| [R] | FHS (2x2, 0.25) |       |       |       |     |     |     | 1.7461e-04 | [5.7678e-02, 8.2316e-02] | 1.5065e-04, 2.3964e-05 | 1362907 | test         |
| [R] | FHS (3x3, 0.25) |       |       |       |     |     |     | 1.7461e-04 | [3.8470e-02, 5.4859e-02] | 1.5065e-04, 2.3964e-05 | 1362908 | test         |

heuristic 4 (test/, 3xsigma everywhere, NL instead of min(N,L)) (nuclear norm seems barely active) (idem)
| Run | Algo            | rw0   | rw 1  | rw 2  | rw3   | rw4   | rw5   | mu         | mu_bar                   | mu_chi, sigma_chi      | id      |
| --- | --------------- | ----- | ----- | ----- | ----- | ----- | ----- | ---------- | ------------------------ | ---------------------- | ------- |
| [K] | HS              | 22.94 |       |       |       |       |       | 2.2254e-04 |                          |                        | 1362470 |
| [x] | FHS (4x4, 0.25) | 22.94 | 32.27 | 36.91 | 38.66 | 39.33 | 39.62 | 2.2254e-04 | [6.8377e-05, 8.5353e-05] | 1.5065e-04, 2.3964e-05 | 1355609 |
| [x] | FHS (2x2, 0.25) | 22.94 | 32.27 | 36.91 | 38.66 | 39.33 | 39.62 | 2.2254e-04 | [7.9658e-05, 8.5315e-05] | 1.5065e-04, 2.3964e-05 | 1355610 |
| [x] | FHS (3x3, 0.25) | 22.94 | 32.27 | 36.91 | 38.66 | 39.33 | 39.62 | 2.2254e-04 | [6.8389e-05, 8.5349e-05] | 1.5065e-04, 2.3964e-05 | 1355611 |

Q=2, ovl=2.50e-01
cw, asnr = 39.62 (2.08e+00), asnr_log = 19.69 (1.43e+00), cpu_cores = 25, iteration_number = 30
$, asnr = 2731.00 (arun (s) = 24.97 (2.10e+00), total_runtime (h) = 18.94
acpu (s) = 178.24 (8.87e+00) , total_cpu_time (h) = 135.21
snr_log, 10 upsilon: 2.59e+01 (1.69e+00) , upsilon/10: 1.47e+01 (1.16e+00)
Q=3, ovl=2.50e-01
cw, asnr = 39.62 (2.08e+00), asnr_log = 19.69 (1.43e+00), cpu_cores = 25, iteration_number = 30
$, asnr = 2730.00 (arun (s) = 14.85 (8.11e-01), total_runtime (h) = 11.26
acpu (s) = 213.07 (7.33e+00) , total_cpu_time (h) = 161.58
snr_log, 10 upsilon: 2.59e+01 (1.69e+00) , upsilon/10: 1.47e+01 (1.16e+00)
Q=4, ovl=2.50e-01
cw, asnr = 39.62 (2.08e+00), asnr_log = 19.69 (1.43e+00), cpu_cores = 25, iteration_number = 30
$, asnr = 2730.00 (arun (s) = 11.89 (5.27e-01), total_runtime (h) = 9.01
acpu (s) = 287.04 (1.43e+01) , total_cpu_time (h) = 217.67
snr_log, 10 upsilon: 2.59e+01 (1.69e+00) , upsilon/10: 1.47e+01 (1.16e+00)

SARA, 1xsigma
rw 5
sara, asnr = 30.25 (3.61e-01), asnr_log = 13.59 (2.85e-01), cpu_cores = 260, iteration_number = 5190
arun (s) = 3.36 (1.77e-01), total_runtime (h) = 5.40
acpu (s) = 7.19 (4.89e-01) , total_cpu_time (h) = 207.45
snr_log, 10 upsilon: 1.84e+01 (6.20e-01) , upsilon/10: 9.05e+00 (2.69e-01)

rw 3
sara, asnr = 28.21 (4.20e-01), asnr_log = 12.13 (3.81e-01), cpu_cores = 260, iteration_number = 4419
arun (s) = 3.36 (1.79e-01), total_runtime (h) = 4.69
acpu (s) = 7.20 (4.95e-01) , total_cpu_time (h) = 176.74
snr_log, 10 upsilon: 1.67e+01 (8.71e-01) , upsilon/10: 7.75e+00 (2.53e-01)

SARA 3xsigma
rw 5
sara, asnr = 36.28 (6.77e-01), asnr_log = 16.90 (2.35e-01), cpu_cores = 260, iteration_number = 2590
arun (s) = 3.31 (1.96e-01), total_runtime (h) = 2.69
acpu (s) = 7.17 (5.00e-01) , total_cpu_time (h) = 103.09
snr_log, 10 upsilon: 2.26e+01 (4.20e-01) , upsilon/10: 1.21e+01 (1.62e-01)

rw 3
sara, asnr = 35.10 (5.10e-01), asnr_log = 16.20 (2.16e-01), cpu_cores = 260, iteration_number = 2376
arun (s) = 3.32 (1.72e-01), total_runtime (h) = 2.46
acpu (s) = 7.18 (5.04e-01) , total_cpu_time (h) = 94.82
snr_log, 10 upsilon: 2.17e+01 (4.54e-01) , upsilon/10: 1.13e+01 (1.51e-01)

in test_configA, heuristic3, 3xsigma everywhere
| Run | Algo            | rw0 | rw 1 | rw 2 | rw3 | rw4 | rw5 | mu  | mu_bar | mu_chi, sigma_chi | id  |
| --- | --------------- | --- | ---- | ---- | --- | --- | --- | --- | ------ | ----------------- | --- |
| [R] | HS (v4)         |     |      |      |     |     |     |     |        |                   |     |
| [R] | FHS (4x4, 0.25) |     |      |      |     |     |     |     |        |                   |     |
| [R] | FHS (2x2, 0.25) |     |      |      |     |     |     |     |        |                   |     |
| [R] | FHS (3x3, 0.25) |     |      |      |     |     |     |     |        |                   |     |

## Spectral faceting (config C only, 2k iterations per pdfb)

heuristic2
| Run | Algo                 | rw 3  | rw 5             | asnr_log (x10, /10)        | arun  | acpu  | total run (rw3, rw5) | total_cpu (rw3, rw5) | id              |
| --- | -------------------- | ----- | ---------------- | -------------------------- | ----- | ----- | -------------------- | -------------------- | --------------- |
| [x] | SARA (Qc=1, 3 sigma) | 19.69 | 19.73 (3.18e+00) | 23.38 (2.43e+01, 1.96e+01) | 0.44  | 0.76  | 0.30, 0.32           | 43.73, 46.93         | 1333870-1333969 | ok    |
| [R] | HS (Qc=1, 3 sigma)   |       |                  |                            |       |       |                      |                      | 1332234         | rerun |
| [x] | HS (Qc=2, 3 sigma)   | 23.35 | 23.59            | 27.89 (2.82e+01, 2.56e+01) | 14.36 | 36.83 | 29.85, 39.98         | 151.83, 199.25       | 1332673-1332674 | rerun |
| [x] | HS (Qc=4, 3 sigma)   | 23.40 | 23.63            | 27.67 (2.83e+01, 2.51e+01) | 7.27  | 19.09 | 15.64, 20.19         | 155.09, 202.09       | 1332676-1332678 | rerun |
| [R] | HS (Qc=5, 3 sigma)   |       |                  |                            |       |       |                      |                      | 1334321-1334325 |
| [x] | HS (Qc=10, 3 sigma)  | 22.55 | 22.63            | 26.95 (2.79e+01, 2.30e+01) | 3.03  | 7.46  | 4.24, 5.94           | 72.19, 79.44         | 1332679-1332688 | ok    |

heuristic 3, 3xsigma, 3xsigma (10 cores for the data + 2, v4)
| Run | Algo                 | rw 3 | rw 5 | asnr_log (x10, /10) | arun | acpu | total run (rw3, rw5) | total_cpu (rw3, rw5) | id              |
| --- | -------------------- | ---- | ---- | ------------------- | ---- | ---- | -------------------- | -------------------- | --------------- |
| [ ] | SARA (Qc=1, 3 sigma) |      |      |                     |      |      |                      |                      |                 |
| [R] | HS (Qc=1, 3 sigma)   |      |      |                     |      |      |                      |                      | 1363961         |
| [R] | HS (Qc=2, 3 sigma)   |      |      |                     |      |      |                      |                      | 1363962-1363963 |
| [R] | HS (Qc=5, 3 sigma)   |      |      |                     |      |      |                      |                      | 1363964-1363968 |
| [R] | HS (Qc=10, 3 sigma)  |      |      |                     |      |      |                      |                      | 1363969-1363978 |
