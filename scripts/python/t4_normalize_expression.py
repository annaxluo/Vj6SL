#!/usr/bin/env python3

import argparse, os, sys, glob, re
import pandas as pd


def read_sample_list(sample_list_file):
    """
    Read sample_list.txt . 
    """

    with open(sample_list_file, "r") as f:
        samples = [line.rstrip("\n") for line in f]

    return samples

def parse_samplesheet(samplesheet):
    """
    parse sample information. 
    """
    import csv
 
    sample_info = {}   
    with open(samplesheet, "r") as file: 
    
        reader = csv.DictReader(file)
        sample_keys = reader.fieldnames[0]
        sample_info = {row[sample_keys]: row for row in reader}
    
    return sample_info


def read_mapped_reads(mapped_reads_file):
    """
    Expected format:
        sample_id    mapped_reads
        sample1      12345678
    """

    df = pd.read_csv(mapped_reads_file, sep="\t")

    out = {}
    for _, row in df.iterrows():
        sid = str(row["sample_id"]).strip()
        mapped = float(row["mapped_reads"])

        out[sid] = mapped

    return out

def cpm(count, mapped_reads):
    """
    Counts per million mapped reads.
    """

    if mapped_reads is None:
        return float("nan")

    if pd.isna(mapped_reads):
        return float("nan")

    if mapped_reads <= 0:
        return float("nan")

    return count / mapped_reads * 1e6


def normalize_one_sample(sample_id, condition, bam_file, mapped_reads, ase_file, min_depth):
    """
    Normalize reads for one sample. 
    """
    
    df = pd.read_csv(ase_file, sep="\t", comment="#")

    rows = []

    for _, row in df.iterrows():

        contig = str(row["contig"])
        position = int(row["position"])
        variant_id = str(row["variantID"])
        ref_allele = str(row["refAllele"])
        alt_allele = str(row["altAllele"])

        ref_count = int(row["refCount"])
        alt_count = int(row["altCount"])
        total_count = int(row["totalCount"])
        low_mapq_depth = int(row["lowMAPQDepth"])
        low_baseq_depth = int(row["lowBaseQDepth"])
        raw_depth = int(row["rawDepth"])
        other_bases = int(row["otherBases"])
        improper_pairs = int(row["improperPairs"])

        depth_ = ref_count + alt_count

        if depth_ > 0:
            ref_fraction_ref_alt = ref_count / depth_
            alt_fraction_ref_alt = alt_count / depth_
        else:
            ref_fraction_ref_alt = float("nan")
            alt_fraction_ref_alt = float("nan")

        if total_count > 0:
            ref_fraction_total = ref_count / total_count
            alt_fraction_total = alt_count / total_count
            other_bases_fraction_total = other_bases / total_count
        else:
            ref_fraction_total = float("nan")
            alt_fraction_total = float("nan")
            other_bases_fraction_total = float("nan")

        if raw_depth > 0:
            low_mapq_fraction_raw = low_mapq_depth / raw_depth
            low_baseq_fraction_raw = low_baseq_depth / raw_depth
            improper_pairs_fraction_raw = improper_pairs / raw_depth
        else:
            low_mapq_fraction_raw = float("nan")
            low_baseq_fraction_raw = float("nan")
            improper_pairs_fraction_raw = float("nan")

        if mapped_reads is None or pd.isna(mapped_reads) or mapped_reads <= 0:
            status = "missing_or_zero_mapped_reads"
        elif depth_ < min_depth:
            status = "low_depth"
        else:
            status = "ok"

        rows.append({
            "sample_id": sample_id,
            "condition": condition,
            "bam_file": bam_file,

            "contig": contig,
            "position": position,
            "variantID": variant_id,
            "refAllele": ref_allele,
            "altAllele": alt_allele,

            "mapped_reads": mapped_reads,

            "REF": ref_count,
            "ALT": alt_count,
            "REF_ALT_total": depth_,

            "totalCount": total_count,
            "rawDepth": raw_depth,
            "otherBases": other_bases,
            "lowMAPQDepth": low_mapq_depth,
            "lowBaseQDepth": low_baseq_depth,
            "improperPairs": improper_pairs,

            "REF_CPM": cpm(ref_count, mapped_reads),
            "ALT_CPM": cpm(alt_count, mapped_reads),
            "REF_ALT_total_CPM": cpm(depth_, mapped_reads),
            "totalCount_CPM": cpm(total_count, mapped_reads),
            "rawDepth_CPM": cpm(raw_depth, mapped_reads),
            "otherBases_CPM": cpm(other_bases, mapped_reads),

            "REF_fraction_REF_ALT": ref_fraction_ref_alt,
            "ALT_fraction_REF_ALT": alt_fraction_ref_alt,

            "REF_fraction_totalCount": ref_fraction_total,
            "ALT_fraction_totalCount": alt_fraction_total,
            "otherBases_fraction_totalCount": other_bases_fraction_total,

            "lowMAPQ_fraction_rawDepth": low_mapq_fraction_raw,
            "lowBaseQ_fraction_rawDepth": low_baseq_fraction_raw,
            "improperPairs_fraction_rawDepth": improper_pairs_fraction_raw,

            "min_depth": min_depth,
            "status": status,
            "ase_file": ase_file,
        })

    return rows


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Normalize REF and ALT variant reads to total library size."
        )
    )

    parser.add_argument(
        "--ase-count-dir",
        required=True,
        help="path to ASEReadCounter output tsv files. "
    )

    parser.add_argument(
        "--sample-list",
        required=True,
        help="sample_list.txt with one sample ID per line."
    )

    parser.add_argument(
        "--mapped-reads",
        required=True,
        help="TSV file with columns: sample_id, mapped_reads."
    )

    parser.add_argument(
        "--samplesheet",
        required=True,
        help="samplesheet.csv with sample_id, bam_file, condition."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="output filename for normalized ASE expression TSV."
    )

    parser.add_argument(
        "--min-depth",
        type=int,
        default=10,
        help="minumum REF + ALT count required for status=ok. Default: 10."
    )

    args = parser.parse_args()

    # read sample info
    samples = read_sample_list(args.sample_list)
    sample_info = parse_samplesheet(args.samplesheet)
    mapped_reads_dict = read_mapped_reads(args.mapped_reads)

    all_rows = []

    # process each sample   
    for sample_id in samples:
    
        condition = sample_info[sample_id].get("condition", "NA")
        bam_file = sample_info[sample_id].get("bam_file", "NA")
        mapped_reads = mapped_reads_dict[sample_id]
        ase_file =  os.path.join(args.ase_count_dir, f"{sample_id}.ASEReadCounter.tsv")

        rows = normalize_one_sample(
            sample_id=sample_id,
            condition=condition,
            bam_file=bam_file,
            mapped_reads=mapped_reads,
            ase_file=ase_file,
            min_depth=args.min_depth
        )

        all_rows.extend(rows)

    out = pd.DataFrame(all_rows)

    output_dir = os.path.dirname(args.output)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()