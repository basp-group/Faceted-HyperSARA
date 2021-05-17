# Faceting experiments

## Preliminary tests

Image size: 512 x 1024 x 20, 2x2 facets for FHS, no homotopy

| Run | Algo        | Reg. type | Reg. update | $\alpha$ | $\bar{\alpha}$ | $\upsilon$ | $\bar{\upsilon}$               | $\mu$                     | $\bar{\mu}$               | aSNR  | id     |
| --- | ----------- | --------- | ----------- | -------- | -------------- | ---------- | ------------------------------ | ------------------------- | ------------------------- | ----- | ------ |
| [ ] | SARA (ch 1) |           |             | 1e-1     | -              |            | -                              |                           | -                         |       |        |
| [ ] | SARA (ch 1) |           |             | 1        | -              |            | -                              |                           | -                         |       |        |
| [ ] | SARA (ch 1) |           |             | 10       | -              |            | -                              |                           | -                         |       |        |
| [ ] | HyperSARA   |           |             | 1e-1     | 1e-1           |            |                                |                           |                           |       |        |
| [ ] | HyperSARA   |           |             | 1e-1     | 1              |            |                                |                           |                           |       |        |
| [ ] | HyperSARA   |           |             | 1e-1     | 10             |            |                                |                           |                           |       |        |
| [ ] | HyperSARA   |           |             | 1        | 1e-1           |            |                                |                           |                           |       |        |
| [ ] | HyperSARA   |           |             | 1        | 1              |            |                                |                           |                           |       |        |
| [ ] | HyperSARA   |           |             | 1        | 10             |            |                                |                           |                           |       |        |
| [ ] | HyperSARA   |           |             | 10       | 1e-1           |            |                                |                           |                           |       |        |
| [ ] | HyperSARA   |           |             | 10       | 1              |            |                                |                           |                           |       |        |
| [ ] | HyperSARA   |           |             | 10       | 10             |            |                                |                           |                           |       |        |
| [X] | FHS (2x2)   | inv       | 0           | 1e-1     | 1e-1           | 2.9079e-07 | [5.8051e-07, 1.3771e-06]       | 1.2522e-04 / -            | 1.0555e-02 /  -           | 29.37 | 795008 | -> SNR slightly better after update |
| [X] | FHS (2x2)   | inv       | 0           | 1e-1     | 1              | "          |                                | "                         | "                         | 29.45 | 795009 |
| [R] | FHS (2x2)   | inv       | 0           | 1e-1     | 10             | "          |                                | "                         | "                         | 29.97 | 795010 |
| [R] | FHS (2x2)   | inv       | 0           | 1        | 1e-1           | "          |                                | "                         | "                         | 31.38 | 795020 |
| [X] | FHS (2x2)   | inv       | 0           | 1        | 1              | "          |                                | "                         | "                         | 31.14 | 794488 |
| [R] | FHS (2x2)   | inv       | 0           | 1        | 10             | "          |                                | "                         | "                         | 28.50 | 795021 |
| [X] | FHS (2x2)   | inv       | 0           | 10       | 1e-1           | "          |                                | "                         | "                         | 31.30 | 795011 |
| [R] | FHS (2x2)   | inv       | 0           | 10       | 1              | "          |                                | "                         | "                         | 28.61 | 795012 | restarted                           |
| [R] | FHS (2x2)   | inv       | 0           | 10       | 10             | "          |                                | "                         | "                         | 33.07 | 795013 |
| [x] | FHS (2x2)   | inv       | 1           | 1e-1     | 1e-1           | 2.9079e-07 | [5.8051e-07, 1.3771e-06] / ... | 1.2522e-04 / 2.084401e-04 | 1.0555e-02 / 1.187855e-04 | 29.37 | 795014 |
| [X] | FHS (2x2)   | inv       | 1           | 1e-1     | 1              | "          | " /                            | " /2.101448e-04           | " / 1.204080e-02          | 29.45 | 795015 |
| [R] | FHS (2x2)   | inv       | 1           | 1e-1     | 10             | "          | " /                            | " / 2.112602e-04          | " / 1.211665e-02          | 30.07 | 795016 |
| [R] | FHS (2x2)   | inv       | 1           | 1        | 1e-1           | "          | " /                            | " / 2.138199e-04          | " / 1.196869e-02          | 31.41 | 795022 | restarted                           |
| [x] | FHS (2x2)   | inv       | 1           | 1        | 1              | "          | " /                            | " / 2.140677e-04          | " / 1.204870e-02          | 31.19 | 794478 |
| [R] | FHS (2x2)   | inv       | 1           | 1        | 10             | "          | " /                            | " / 2.147869e-04          | " / 1.211759e-02          | 28.50 | 795136 | restarted                           |
| [x] | FHS (2x2)   | inv       | 1           | 10       | 1e-1           | "          | " /                            | " / 2.186951e-04          | " / 1.202126e-02          | 31.64 | 795017 |
| [R] | FHS (2x2)   | inv       | 1           | 10       | 1              | "          | " /                            | " / 2.187030e-04          | " / 1.206262e-02          | 28.61 | 795018 | restarted                           |
| [R] | FHS (2x2)   | inv       | 1           | 10       | 10             | "          | " /                            | " / 2.182731e-04          | " / 1.210817e-02          | 33.34 | 795019 | rel. variation starts to oscillate  |

-> take (1,1) (much more stable in the end), no update really necessary

---

### Order of magnitude regularization parameters

| Reg. type | Image approx. | Noise transfer | $\upsilon$ | $\bar{\upsilon}$ | $\mu$ | $\bar{\mu}$ | id  |
| --------- | ------------- | -------------- | ---------- | ---------------- | ----- | ----------- | --- |
| log       | none          | none           |            |                  |       |             |     |
| log       | none          | precond        |            |                  |       |             |     |
| log       | precond       | none           |            |                  |       |             |     |
| log       | precond       | precond        |            |                  |       |             |     |

---

## Spatial faceting

Image size: 1024 x 2048 x 20.

---

## Spectral faceting

Image size: 256 x 512 x 100.
