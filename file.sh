#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

#########################################
# CONFIGURATION
#########################################

SRR_LIST=("SRR32915532" "SRR32302310" "SRR31734152" "SRR36272134")
THREADS=1            # Low-memory: use 1 thread
WORKDIR="$PWD/sra_cytoscape_lowmem"

REF_DIR="$WORKDIR/reference"
FASTQ_DIR="$WORKDIR/fastq"
FASTQC_DIR="$WORKDIR/fastqc_reports"
ALIGNED_DIR="$WORKDIR/aligned"
COUNTS_DIR="$WORKDIR/counts"

mkdir -p "$REF_DIR" "$FASTQ_DIR" "$FASTQC_DIR" "$ALIGNED_DIR" "$COUNTS_DIR"

# Ensembl URLs for GRCh38 primary assembly and GTF
FASTA_URL="ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"

REF_GENOME="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="$REF_DIR/Homo_sapiens.GRCh38.109.gtf"

#########################################
# STEP 0: Download reference genome if missing
#########################################

if [ ! -f "$REF_GENOME" ]; then
    echo "Downloading GRCh38 primary assembly FASTA..."
    wget -O "$REF_GENOME.gz" "$FASTA_URL"
    gunzip -f "$REF_GENOME.gz"
fi

if [ ! -f "$GTF_FILE" ]; then
    echo "Downloading GRCh38 GTF annotation..."
    wget -O "$GTF_FILE.gz" "$GTF_URL"
    gunzip -f "$GTF_FILE.gz"
fi

#########################################
# STEP 1: Stream SRA to FASTQ directly (low-memory)
#########################################

echo "Downloading SRA and converting to FASTQ (streaming)..."
for srr in "${SRR_LIST[@]}"; do
    echo "Processing $srr"
    # Use fastq-dump with streaming from SRA toolkit without storing .sra
    fasterq-dump "$srr" \
        -O "$FASTQ_DIR" \
        --split-3 \
        --skip-technical \
        --threads "$THREADS" \
        --mem 500 \
        --progress
done

#########################################
# STEP 2: Quality Control with FastQC
#########################################

echo "Running FastQC..."
for fq in "$FASTQ_DIR"/*.fastq; do
    [ -s "$fq" ] || continue
    fastqc "$fq" -o "$FASTQC_DIR"
done

#########################################
# STEP 3: Build HISAT2 index and align
#########################################

if [ ! -f "$REF_DIR/hg38_index.1.ht2" ]; then
    echo "Building HISAT2 index..."
    hisat2-build -p "$THREADS" "$REF_GENOME" "$REF_DIR/hg38_index"
fi

echo "Aligning reads with HISAT2 (streaming, low-memory)..."
for fq in "$FASTQ_DIR"/*.fastq; do
    [ -s "$fq" ] || continue
    base=$(basename "$fq" .fastq)
    sorted_bam="$ALIGNED_DIR/$base.sorted.bam"

    hisat2 -p "$THREADS" --no-spliced-alignment -x "$REF_DIR/hg38_index" -U "$fq" \
        | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "$sorted_bam"

    samtools index "$sorted_bam"
done

#########################################
# STEP 4: Count reads with FeatureCounts
#########################################

echo "Counting reads per gene..."
for bam in "$ALIGNED_DIR"/*.sorted.bam; do
    featureCounts -T "$THREADS" -a "$GTF_FILE" -o "$COUNTS_DIR/$(basename "$bam" .sorted.bam)_counts.txt" \
        -M --fraction "$bam"
done

#########################################
# STEP 5: Merge counts with Python
#########################################

python3 - <<'EOF'
import pandas as pd
import glob, os

counts_dir = "counts"
files = glob.glob(os.path.join(counts_dir, "*_counts.txt"))
merged = None

for f in files:
    df = pd.read_csv(f, sep='\t', comment='#', index_col=0)
    df = df.iloc[:,5:]  # drop annotation columns
    df.columns = [os.path.basename(f).replace('_counts.txt','')]
    merged = df if merged is None else merged.join(df, how='outer')

merged.fillna(0, inplace=True)
merged.to_csv(os.path.join(counts_dir, "counts_merged.txt"))
EOF

#########################################
# STEP 6: DEGs + KEGG enrichment
#########################################

cat > "$WORKDIR/run_DEGs_KEGG.py" << 'EOF'
import pandas as pd, numpy as np, os
from scipy.stats import ttest_ind
from gseapy import enrichr

counts = pd.read_csv("counts/counts_merged.txt", sep=',', index_col=0)
conditions = ['control','control','treated','treated']
counts.columns = conditions

control_cols = [c for c in counts.columns if 'control' in c]
treated_cols = [c for c in counts.columns if 'treated' in c]

deg_results = []
for gene in counts.index:
    c_vals = counts.loc[gene, control_cols].values
    t_vals = counts.loc[gene, treated_cols].values
    if len(c_vals) > 1 and len(t_vals) > 1:
        t_stat, p_val = ttest_ind(t_vals, c_vals)
        log2fc = np.log2((t_vals.mean()+1)/(c_vals.mean()+1))
    else:
        t_stat, p_val, log2fc = np.nan, np.nan, np.nan
    deg_results.append([gene, log2fc, p_val])

deg_df = pd.DataFrame(deg_results, columns=['Gene','log2FC','pvalue'])
deg_df['padj'] = deg_df['pvalue'] * len(deg_df)
deg_df.to_csv('DEG_results.csv', index=False)

deg_genes = deg_df[(deg_df['padj']<0.05) & (deg_df['log2FC'].abs()>1)]['Gene'].tolist()
if deg_genes:
    os.makedirs("KEGG_results", exist_ok=True)
    enrichr(gene_list=deg_genes, gene_sets=['KEGG_2021_Human'], organism='Human',
             description='KEGG_enrichment', outdir='KEGG_results', cutoff=0.5)
    print("KEGG enrichment completed. Results in KEGG_results/")
else:
    print("No DEGs passed thresholds; skipping KEGG enrichment.")
EOF

#########################################
# STEP 7: Run DEGs + KEGG
#########################################

echo "Running DEGs + KEGG analysis..."
python3 "$WORKDIR/run_DEGs_KEGG.py"

echo "Pipeline complete!"
echo "DEGs in DEG_results.csv"
echo "KEGG results in KEGG_results/"
