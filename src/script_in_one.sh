#!/bin/bash

echo "C survival study"
Rscript src/c_figures.R

echo "Analyzing off targets C"
Rscript src/plot_off_targets_C.R

echo "Analyzing off targets Nextseq"
Rscript src/plot_off_targets_nextseq_35_100dpi.R

echo "Analyzing off targets HEK BE39max treated"
Rscript src/plot_off_targets_HEK_BE3.9max.R

echo "Analyzing off targets HEK TadCBEd treated"
Rscript src/plot_off_targets_HEK_TadCBEd.R

echo "Analyzing off targets VEP result"
Rscript src/vep_result.R


echo "All plots and summary are generated."
