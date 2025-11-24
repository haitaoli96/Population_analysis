# Population_analysis

好嘞，我直接按你给的**绝对路径**把一套“C. parvum 单倍体 GATK SNP calling pipeline”改成**可一键批量跑**的版本。你复制到服务器即可跑。

我默认：

* 你的参考是染色体级别 `*.fna`（FASTA）
* fastq 文件在一个目录里，命名可能是 `_1/_2` 或 `_R1/_R2`
* 走 **BWA-MEM2 → MarkDuplicates → HaplotypeCaller(gVCF, ploidy=1) → Joint Genotype → Hard Filter**

---

# 0) 你需要先建好环境（一次）

```bash
conda create -n SNP -c conda-forge -c bioconda \
  gatk4 bwa-mem2 samtools picard fastp bcftools tabix -y
conda activate SNP
```

---

# 1) 写一个总控脚本：`run_gatk_cp.sh`

把下面整段保存为 `run_gatk_cp.sh`（路径随便放，比如你项目根目录）。
**我已经把你的 REF 和 FASTQ 路径填进去了。**

```bash
#!/usr/bin/env bash
set -euo pipefail

########################################
# 0. 固定输入路径（你给的）
########################################
REF="/home/liht/lihaitao/04.Cryptosporidium_refgenome/Cryptosporidium_parvum_IOWA-ATCC/GCA_015245375.1_ASM1524537v1_genomic.fna"
RAW="/home/liht/NCBI_sra_data/C_parvum/fastq_downloads"

########################################
# 1. 输出目录（自动创建）
########################################
WORKDIR=$(pwd)
CLEAN="$WORKDIR/clean_fastq"
BAM="$WORKDIR/bam"
GVCF="$WORKDIR/gvcf"
JOINT="$WORKDIR/joint"
FILTER="$WORKDIR/filter"

LOG_FASTP="$WORKDIR/logs_fastp"
LOG_MAP="$WORKDIR/logs_map"
LOG_GVCF="$WORKDIR/logs_gvcf"
LOG_JOINT="$WORKDIR/logs_joint"
LOG_FILTER="$WORKDIR/logs_filter"

mkdir -p "$CLEAN" "$BAM" "$GVCF" "$JOINT" "$FILTER"
mkdir -p "$LOG_FASTP" "$LOG_MAP" "$LOG_GVCF" "$LOG_JOINT" "$LOG_FILTER"

THREADS=12
PLOIDY=1   # C. parvum 单倍体

echo "=================================================="
echo "REF   = $REF"
echo "RAW   = $RAW"
echo "WORK  = $WORKDIR"
echo "=================================================="

########################################
# Step A. 参考索引（只要缺就建）
########################################
echo "==> [A] Index reference if needed..."

if [[ ! -f "${REF}.bwt.2bit.64" && ! -f "${REF}.0123" ]]; then
  bwa-mem2 index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
  samtools faidx "$REF"
fi

DICT="${REF%.*}.dict"
if [[ ! -f "$DICT" ]]; then
  picard CreateSequenceDictionary R="$REF" O="$DICT"
fi


########################################
# Step B. fastp 清洗（可选但推荐）
########################################
echo "==> [B] fastp trimming..."

shopt -s nullglob
for r1 in "$RAW"/*_1.fastq.gz "$RAW"/*_R1.fastq.gz; do
  base=$(basename "$r1")

  if [[ "$base" == *_1.fastq.gz ]]; then
    sample=${base%%_1.fastq.gz}
    r2="$RAW/${sample}_2.fastq.gz"
  else
    sample=${base%%_R1.fastq.gz}
    r2="$RAW/${sample}_R2.fastq.gz"
  fi

  [[ -f "$r2" ]] || { echo "[WARN] Missing mate for $r1, skip."; continue; }

  out1="$CLEAN/${sample}_R1.clean.fastq.gz"
  out2="$CLEAN/${sample}_R2.clean.fastq.gz"

  if [[ -f "$out1" && -f "$out2" ]]; then
    echo "   fastp done already: $sample"
    continue
  fi

  echo "   -> fastp $sample"
  fastp \
    -i "$r1" -I "$r2" \
    -o "$out1" -O "$out2" \
    -q 20 -u 30 -n 5 -l 50 \
    -w 8 \
    -h "$LOG_FASTP/${sample}.html" \
    -j "$LOG_FASTP/${sample}.json" \
    > "$LOG_FASTP/${sample}.log" 2>&1
done


########################################
# Step C. 比对 + 排序 + MarkDuplicates
########################################
echo "==> [C] Mapping + MarkDuplicates..."

for r1 in "$CLEAN"/*_R1.clean.fastq.gz; do
  sample=$(basename "$r1" _R1.clean.fastq.gz)
  r2="$CLEAN/${sample}_R2.clean.fastq.gz"

  sorted="$BAM/${sample}.sorted.bam"
  mkdup="$BAM/${sample}.mkdup.bam"

  if [[ ! -f "$mkdup" ]]; then
    echo "   -> bwa-mem2 mem $sample"
    bwa-mem2 mem -t $THREADS "$REF" "$r1" "$r2" \
      | samtools sort -@ $THREADS -o "$sorted"

    samtools index "$sorted"

    echo "   -> MarkDuplicates $sample"
    picard MarkDuplicates \
      I="$sorted" \
      O="$mkdup" \
      M="$LOG_MAP/${sample}.dup.metrics.txt" \
      VALIDATION_STRINGENCY=SILENT \
      REMOVE_DUPLICATES=false \
      > "$LOG_MAP/${sample}.markdup.log" 2>&1

    samtools index "$mkdup"
    samtools flagstat "$mkdup" > "$LOG_MAP/${sample}.flagstat.txt"
  else
    echo "   mkdup done already: $sample"
  fi
done


########################################
# Step D. HaplotypeCaller (gVCF, haploid)
########################################
echo "==> [D] HaplotypeCaller gVCF..."

for bam in "$BAM"/*.mkdup.bam; do
  sample=$(basename "$bam" .mkdup.bam)
  out="$GVCF/${sample}.g.vcf.gz"

  if [[ -f "$out" ]]; then
    echo "   gVCF done already: $sample"
    continue
  fi

  echo "   -> HaplotypeCaller $sample"
  gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R "$REF" \
    -I "$bam" \
    -O "$out" \
    -ERC GVCF \
    --sample-ploidy $PLOIDY \
    --native-pair-hmm-threads 8 \
    > "$LOG_GVCF/${sample}.hc.log" 2>&1
done


########################################
# Step E. 联合分型 CombineGVCFs + GenotypeGVCFs
########################################
echo "==> [E] Joint genotyping..."

ls "$GVCF"/*.g.vcf.gz > "$JOINT/gvcf.list"

if [[ ! -f "$JOINT/combined.g.vcf.gz" ]]; then
  echo "   -> CombineGVCFs"
  gatk --java-options "-Xmx16g" CombineGVCFs \
    -R "$REF" \
    $(awk '{print "-V "$1}' "$JOINT/gvcf.list") \
    -O "$JOINT/combined.g.vcf.gz" \
    > "$LOG_JOINT/combine.log" 2>&1
fi

if [[ ! -f "$JOINT/raw.vcf.gz" ]]; then
  echo "   -> GenotypeGVCFs"
  gatk --java-options "-Xmx16g" GenotypeGVCFs \
    -R "$REF" \
    -V "$JOINT/combined.g.vcf.gz" \
    -O "$JOINT/raw.vcf.gz" \
    > "$LOG_JOINT/genotype.log" 2>&1
fi


########################################
# Step F. SNP/INDEL 分离 + Hard Filter
########################################
echo "==> [F] Hard filtering..."

RAWVCF="$JOINT/raw.vcf.gz"

if [[ ! -f "$FILTER/raw.snps.vcf.gz" ]]; then
  gatk SelectVariants \
    -R "$REF" -V "$RAWVCF" \
    --select-type-to-include SNP \
    -O "$FILTER/raw.snps.vcf.gz" \
    > "$LOG_FILTER/select_snp.log" 2>&1
fi

if [[ ! -f "$FILTER/raw.indels.vcf.gz" ]]; then
  gatk SelectVariants \
    -R "$REF" -V "$RAWVCF" \
    --select-type-to-include INDEL \
    -O "$FILTER/raw.indels.vcf.gz" \
    > "$LOG_FILTER/select_indel.log" 2>&1
fi

if [[ ! -f "$FILTER/filtered.snps.vcf.gz" ]]; then
  gatk VariantFiltration \
    -R "$REF" -V "$FILTER/raw.snps.vcf.gz" \
    --filter-name "QD2"       --filter-expression "QD < 2.0" \
    --filter-name "FS60"      --filter-expression "FS > 60.0" \
    --filter-name "MQ40"      --filter-expression "MQ < 40.0" \
    --filter-name "SOR3"      --filter-expression "SOR > 3.0" \
    --filter-name "MQRankSum-12.5" --filter-expression "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum-8" --filter-expression "ReadPosRankSum < -8.0" \
    -O "$FILTER/filtered.snps.vcf.gz" \
    > "$LOG_FILTER/filter_snp.log" 2>&1
fi

if [[ ! -f "$FILTER/filtered.indels.vcf.gz" ]]; then
  gatk VariantFiltration \
    -R "$REF" -V "$FILTER/raw.indels.vcf.gz" \
    --filter-name "QD2"  --filter-expression "QD < 2.0" \
    --filter-name "FS200" --filter-expression "FS > 200.0" \
    --filter-name "SOR10" --filter-expression "SOR > 10.0" \
    -O "$FILTER/filtered.indels.vcf.gz" \
    > "$LOG_FILTER/filter_indel.log" 2>&1
fi

# 只保留 PASS SNP
if [[ ! -f "$FILTER/final.snps.vcf.gz" ]]; then
  bcftools view -f PASS -v snps "$FILTER/filtered.snps.vcf.gz" -Oz -o "$FILTER/final.snps.vcf.gz"
  tabix -p vcf "$FILTER/final.snps.vcf.gz"
fi

echo "=================================================="
echo "Pipeline finished!"
echo "Final SNP VCF: $FILTER/final.snps.vcf.gz"
echo "=================================================="
```

---

# 2) 运行方式

```bash
conda activate gatk_cp
bash run_gatk_cp.sh
```

它会在你当前目录生成：

```
clean_fastq/
bam/
gvcf/
joint/
filter/
logs_*/
```

最终结果：

```
filter/final.snps.vcf.gz
```

---

# 3) 你接下来最可能要做的两件事

### A. 做缺失率/MAF过滤（给群体遗传/机器学习矩阵）

```bash
bcftools view -i 'F_MISSING<0.2 && MAF>0.01' \
  filter/final.snps.vcf.gz -Oz -o filter/final.snps.pruned.vcf.gz
tabix -p vcf filter/final.snps.pruned.vcf.gz
```

### B. 导出 SNP 矩阵

```bash
bcftools query -f '%CHROM\t%POS[\t%GT]\n' \
  filter/final.snps.pruned.vcf.gz > snp_matrix.txt
```

---

# 4) 小提醒（按你的研究方向很关键）

* 你要做**共感染/重组/平行微进化**这些强创新点的话：
  **统计每个样本杂合位点比例**就是第一张“证据图”。
  我可以在你跑完 VCF 后给你一条命令直接输出 per-sample het rate + 热图。

---

你现在直接跑这个脚本就行。
跑完如果出现任何报错，把 `logs_*/*.log` 最后 30 行贴我，我会按你的 HPC 环境快速定位并修正。
