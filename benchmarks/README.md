# IQ-paper – Figure 5 Benchmarking

This sub-directory contains **all** code and helper scripts required to reproduce the deep-learning benchmarks used in Figure 5 of the manuscript.

## 1. What is shown in Figure 5?
* **5A – IQ model variants** R² obtained with three IQ motif-based models on the mouse gastrulation dataset: (i) 90 motifs w/o interactions (blue), (ii) 90 motifs with pair-wise interactions (green) and (iii) an expanded 180-motif IQ with interactions (red).
* **5B – Observed vs predicted scatter** Density scatter plots of differential ATAC-seq accessibility (dAP) for one example trajectory comparing IQ, Borzoi and DeepTopic CNN.
* **5C – Model comparison per trajectory** Boxplots of per-trajectory R² for Borzoi, IQ and DeepTopic CNN. Grey lines connect the same trajectory across models (paired Wilcoxon, FDR-adj. p = 0.81 – no significant difference).
* **5D – All models & ensembles** Per-trajectory R² for every individual model that was trained (CNNs, Transformers, IQ) together with simple linear ensembles. Models are sorted by mean performance.
* **5E – Error attribution** Kolmogorov–Smirnov D statistics comparing motif energies inside model-specific errors (> 0.1 RMS) against well-predicted regions (< 0.05 RMS).

The current repository focuses on *training & evaluation* (panels 5C–5D).

## 2. Directory layout
```
benchmarks/
├── crested/                 # CREsted training scripts
│   ├── train_crested.py
│   └── train_borzoi_finetune.py
├── scripts/                 # Lightweight helpers 
│   └── collect_r2.py
├── commands.sh              # One run-command per model / dataset
├── requirements.txt         
└── README.md                # (this file)
```

## 3. Quick start
1. **Hardware requirements**
   - **For IQ models**: 32-core CPU system, 30GB RAM (peak usage)
   - **For deep learning benchmarks**: 2x L40 GPUs (or equivalent), additional CPU/RAM for training
   
2. **Install deps**
   ```bash
   python -m venv venv && source venv/bin/activate
   pip install -r benchmarks/requirements.txt
   # The data download script uses the AWS CLI
   # pip install awscli
   ```
3. **Download reference genomes**
   A helper script is provided to download the `hg38` and `mm10` reference genomes into the project's root directory.
   ```bash
   ./benchmarks/scripts/download_genomes.sh
   ```
4. **Download input data**
   The BigWig files and BED peak lists are hosted in a public AWS S3 bucket at:
   - Bucket: `s3://iceqream/`
   - Public URL: `https://iceqream.s3.eu-west-1.amazonaws.com/`
   
   Download the data into the `benchmarks/data/` directory:
   ```bash
   ./benchmarks/scripts/download_data.sh
   ```
   This will download all required files from the S3 bucket. The script requires the AWS CLI to be installed.

5. **Run benchmarks**
   ```bash
   ./benchmarks/commands.sh   # runs *all* commands, OR
   # copy a single command from commands.sh and run individually
   ```
6. **Collect R² values**
   ```bash
   python benchmarks/scripts/collect_r2.py --runs output/* > all_models_r2.csv
   ```

---
Feel free to open an issue if you encounter any problems. 