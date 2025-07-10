#!/usr/bin/env bash
# ------------------------------------------------------------
# Figure 5 – Deep-learning model benchmarking (IQ-paper)
# ------------------------------------------------------------
# This script lists ONE canonical training command per model /
# organism that was **actually used** to generate the numbers in
# Fig. 5C–D.  Edit the two path variables below before running.
# ------------------------------------------------------------

# TODO – add an S3 link or local location that contains all
# normalised bigWig tracks and BED peak sets.  Example:
#   s3://my-bucket/iq-paper/crested_inputs/
DATA_PREFIX="benchmarks/data"        # <—— fill in
OUT_DIR="output"       # results will be created here

# Common CLI fragments -------------------------------------------------------
GENOME_MM10_FASTA="mm10.fa"      # download or symlink locally
GENOME_HG38_FASTA="hg38.fa"
CHROMSIZES_MM10="mm10.chrom.sizes"
CHROMSIZES_HG38="hg38.chrom.sizes"
PROJECT="IQ-comparisons"
LOGGER="wandb"
EPOCHS=60     # CNNs / simple models

##############################################################################
# 1) MOUSE GASTRULATION (6 trajectories)
##############################################################################

# Dilated CNN (ChromBPNet) – train on **cell-state** ATAC values
python benchmarks/crested/train_crested.py \
    --bigwigs-folder  "$DATA_PREFIX/bw-cell-type-no-tss" \
    --regions-file    "$DATA_PREFIX/gastru_peaks_no_tss.bed" \
    --output-dir      "$OUT_DIR/gastru-chrombpnet" \
    --project-name    "$PROJECT" \
    --split-strategy  chr \
    --val-chroms      chr8 chr10 \
    --test-chroms     chr9 chr18 \
    --architecture    chrombpnet \
    --logger          $LOGGER \
    --epochs          $EPOCHS

# DeepTopic CNN – direct **differential** training
python benchmarks/crested/train_crested.py \
    --bigwigs-folder  "$DATA_PREFIX/bw" \
    --regions-file    "$DATA_PREFIX/gastru_peaks.bed" \
    --output-dir      "$OUT_DIR/gastru-deeptopic_cnn-diff" \
    --project-name    "$PROJECT" \
    --split-strategy  chr \
    --val-chroms      chr8 chr10 \
    --test-chroms     chr9 chr18 \
    --architecture    deeptopic_cnn \
    --logger          $LOGGER \
    --epochs          $EPOCHS

# Simple ConvNet – baseline CNN (differential training)
python benchmarks/crested/train_crested.py \
    --bigwigs-folder  "$DATA_PREFIX/bw" \
    --regions-file    "$DATA_PREFIX/gastru_peaks.bed" \
    --output-dir      "$OUT_DIR/gastru-simple_convnet-diff" \
    --project-name    "$PROJECT" \
    --split-strategy  chr \
    --val-chroms      chr8 chr10 \
    --test-chroms     chr9 chr18 \
    --architecture    simple_convnet \
    --logger          $LOGGER \
    --epochs          $EPOCHS

# Transformer (Borzoi) – two-phase finetuning
python benchmarks/crested/train_borzoi_finetune.py \
    --bigwigs-folder  "$DATA_PREFIX/bw-cell-type-no-tss" \
    --regions-file    "$DATA_PREFIX/gastru_peaks_no_tss.bed" \
    --output-dir      "$OUT_DIR/gastru-borzoi-finetune" \
    --project-name    "$PROJECT" \
    --split-strategy  chr \
    --val-chroms      chr8 chr10 \
    --test-chroms     chr9 chr18 \
    --architecture    borzoi \
    --batch-size      256 \
    --finetune \
    --borzoi-model    Borzoi_mouse_rep0 \
    --seq-len         2048 \
    --initial-lr      1e-5 \
    --finetune-lr     5e-5 \
    --top-k-percent   0.03 \
    --two-phase \
    --finetune-epochs 5 \
    --gini-threshold  1.0

##############################################################################
# 2) HUMAN HEMATOPOIESIS (3 trajectories)
##############################################################################

# Dilated CNN (ChromBPNet) – cell-state training
python benchmarks/crested/train_crested.py \
    --bigwigs-folder  "$DATA_PREFIX/bw-bone-marrow" \
    --regions-file    "$DATA_PREFIX/bm_peaks.bed" \
    --output-dir      "$OUT_DIR/bm-chrombpnet" \
    --project-name    "$PROJECT" \
    --split-strategy  chr \
    --val-chroms      chr8 chr10 \
    --test-chroms     chr9 chr18 \
    --chromsizes-file $CHROMSIZES_HG38 \
    --genome-file     $GENOME_HG38_FASTA \
    --architecture    chrombpnet \
    --logger          $LOGGER \
    --epochs          $EPOCHS

# DeepTopic CNN – differential training
python benchmarks/crested/train_crested.py \
    --bigwigs-folder  "$DATA_PREFIX/bw-bone-marrow" \
    --regions-file    "$DATA_PREFIX/bm_peaks.bed" \
    --output-dir      "$OUT_DIR/bm-deeptopic_cnn-diff" \
    --project-name    "$PROJECT" \
    --split-strategy  chr \
    --val-chroms      chr8 chr10 \
    --test-chroms     chr9 chr18 \
    --chromsizes-file $CHROMSIZES_HG38 \
    --genome-file     $GENOME_HG38_FASTA \
    --architecture    deeptopic_cnn \
    --logger          $LOGGER \
    --epochs          $EPOCHS

# Simple ConvNet – differential training
python benchmarks/crested/train_crested.py \
    --bigwigs-folder  "$DATA_PREFIX/bw-bone-marrow" \
    --regions-file    "$DATA_PREFIX/bm_peaks.bed" \
    --output-dir      "$OUT_DIR/bm-simple_convnet-diff" \
    --project-name    "$PROJECT" \
    --split-strategy  chr \
    --val-chroms      chr8 chr10 \
    --test-chroms     chr9 chr18 \
    --chromsizes-file $CHROMSIZES_HG38 \
    --genome-file     $GENOME_HG38_FASTA \
    --architecture    simple_convnet \
    --logger          $LOGGER \
    --epochs          $EPOCHS

# Transformer (Borzoi) – two-phase finetuning
python benchmarks/crested/train_borzoi_finetune.py \
    --bigwigs-folder  "$DATA_PREFIX/bw-bone-marrow" \
    --regions-file    "$DATA_PREFIX/bm_peaks.bed" \
    --output-dir      "$OUT_DIR/bm-borzoi-finetune" \
    --project-name    "$PROJECT" \
    --split-strategy  chr \
    --val-chroms      chr8 chr10 \
    --test-chroms     chr9 chr18 \
    --architecture    borzoi \
    --batch-size      256 \
    --chromsizes-file $CHROMSIZES_HG38 \
    --genome-file     $GENOME_HG38_FASTA \
    --finetune \
    --borzoi-model    Borzoi_human_rep0 \
    --seq-len         2048 \
    --initial-lr      1e-5 \
    --finetune-lr     5e-5 \
    --top-k-percent   0.03 \
    --two-phase \
    --finetune-epochs 5 \
    --gini-threshold  1.0

##############################################################################
# 3) COLLECT R² NUMBERS ------------------------------------------------------
# After every run is finished:
#   python benchmarks/scripts/collect_r2.py --runs ${OUT_DIR}/* > all_models_r2.csv
############################################################################## 