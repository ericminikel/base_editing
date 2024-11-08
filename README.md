This repository holds raw data and source code for a manuscript in progress:

An & Davis 2024. **In vivo base editing extends lifespan of a humanized mouse model of prion disease.** In preparation.

Data from survival studies and off-target analyses run by the Vallabh/Minikel lab are housed here and the relevant figures and supplementary tables can be reproduced by running `. src/script_in_one.sh`. This runs in 10 seconds on a 2021 MacBook Air and will regenerate all the files in `output/` and `display_items/`

The additional off-target analyses in revision require running 3 scripts: `src/01_rhampSeq_dedup.sh`, `src/02_crispresso_analysis.sh`, and `src/03_edit_eff.py`. The output is summarized in `output/PRNP_OT_summary.xlsx`

Software versions used in this analysis:

-   R 4.2.0

    -   tidyverse 2.0.0
    -   dplyr 1.1.4
    -   janitor 2.2.0
    -   openxlsx 4.2.5.2
    -   readxl 1.4.3
    -   survival 3.6.4
    -   BSgenome 1.72.0
    -   GenomicRanges 1.56.1
    -   Biostrings 2.72.1

-   Python 3.9.13

-   [CRISPResso2](https://github.com/pinellolab/CRISPResso2)

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons License" style="border-width:0"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
