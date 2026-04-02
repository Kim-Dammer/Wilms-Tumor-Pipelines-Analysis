# Wilms Tumor Pipelines & Analysis

This repository contains the code and analysis workflows associated with a project investigating the **clonal history** of Wilms tumors.

## Overview

The workflows are designed to process Wilms tumor samples (aligned to **hg19**) on the **DKFZ ODCF cluster**. The pipeline focuses on high-resolution copy number profiling and robust somatic variant calling.

---

## Analysis Pipelines

### 1. Copy Number Profiling
We use the **HMF (Hartwig Medical Foundation)** suite to estimate purity, ploidy, and copy number alterations:
* **Amber:** B-allele frequency (BAF) generation.
* **Cobalt:** Read depth ratios.
* **Purple:** Integration of Amber and Cobalt for final purity and copy number calling.

### 2. Somatic Variant Calling
To ensure high specificity, we employ a consensus approach for SNVs and Indels:
* Runs both **Mutect2** and **Strelka2**.
* Final calls are derived from the **intersection** of both callers.


### Installation
To run the analysis notebook locally, set up a virtual environment:

```bash
# Create and activate environment
python -m venv .venv
source .venv/bin/activate  # On Windows use: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt