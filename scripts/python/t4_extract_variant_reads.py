#!/usr/bin/env python3

import pysam
import sys

bam_path = sys.argv[1]
chrom = sys.argv[2]
pos1 = int(sys.argv[3])  # 1-based position
ref = sys.argv[4].upper()
alt = sys.argv[5].upper()
prefix = sys.argv[6]

pos0 = pos1 - 1  # 0-based coordinates for pysam

bam = pysam.AlignmentFile(bam_path, "rb")

out_ref = open(prefix + ".REF_reads.tsv", "w")
out_alt = open(prefix + ".ALT_reads.tsv", "w")
out_other = open(prefix + ".OTHER_reads.tsv", "w")

header = [
    "read_name",
    "flag",
    "chrom",
    "start_1based",
    "mapq",
    "cigar",
    "strand",
    "base_at_variant",
    "base_quality",
    "query_pos_0based",
    "read_sequence"
]

for out in [out_ref, out_alt, out_other]:
    out.write("\t".join(header) + "\n")

for read in bam.fetch(chrom, pos0, pos0 + 1):
    if read.is_unmapped:
        continue

    # get_aligned_pairs maps read positions to reference positions
    base = None
    bq = None
    qpos_at_variant = None

    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos == pos0:
            qpos_at_variant = qpos
            if qpos is not None:
                base = read.query_sequence[qpos].upper()
                if read.query_qualities is not None:
                    bq = read.query_qualities[qpos]
                else:
                    bq = "NA"
            else:
                base = "DEL"
                bq = "NA"
            break

    if base is None:
        continue

    strand = "-" if read.is_reverse else "+"

    row = [
        read.query_name,
        str(read.flag),
        read.reference_name,
        str(read.reference_start + 1),
        str(read.mapping_quality),
        read.cigarstring,
        strand,
        base,
        str(bq),
        str(qpos_at_variant),
        read.query_sequence
    ]

    if base == ref:
        out_ref.write("\t".join(row) + "\n")
    elif base == alt:
        out_alt.write("\t".join(row) + "\n")
    else:
        out_other.write("\t".join(row) + "\n")

out_ref.close()
out_alt.close()
out_other.close()
bam.close()
