#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    CompoundLocation,
)


# ---------- 工具函式 ----------

def shift_location(location, offset):
    """
    將 FeatureLocation / CompoundLocation 平移 offset
    保留 join / complement 結構，只改數字
    """
    if isinstance(location, FeatureLocation):
        return FeatureLocation(
            location.start - offset,
            location.end - offset,
            strand=location.strand
        )

    elif isinstance(location, CompoundLocation):
        new_parts = []
        for part in location.parts:
            new_parts.append(
                FeatureLocation(
                    part.start - offset,
                    part.end - offset,
                    strand=part.strand
                )
            )
        return CompoundLocation(
            new_parts,
            operator=location.operator  # 通常是 "join"
        )

    else:
        raise TypeError("Unknown location type")


# ---------- 核心邏輯 ----------

def extract_genes_gbk(
    input_gbk,
    output_gbk,
    start_gene,
    end_gene,
    gene_key_priority=("locus_tag", "gene"),
):
    record = SeqIO.read(input_gbk, "genbank")

    start_feature = None
    end_feature = None

    # 找起始與結束基因
    for feature in record.features:
        if feature.type != "CDS":
            continue

        for key in gene_key_priority:
            if key in feature.qualifiers:
                gene_id = feature.qualifiers[key][0]
                if gene_id == start_gene:
                    start_feature = feature
                if gene_id == end_gene:
                    end_feature = feature

    if not start_feature or not end_feature:
        raise ValueError("找不到 start_gene 或 end_gene，請確認 gene ID")

    # 決定擷取區間（0-based, end-exclusive）
    start_pos = min(start_feature.location.start, end_feature.location.start)
    end_pos = max(start_feature.location.end, end_feature.location.end)

    # 擷取序列
    new_seq = record.seq[start_pos:end_pos]

    new_record = SeqRecord(
        new_seq,
        id=f"{record.id}_{start_gene}_{end_gene}",
        name=record.name,
        description=f"Extracted genes {start_gene} to {end_gene}",
    )

    new_record.annotations = record.annotations.copy()
    new_record.annotations["molecule_type"] = "DNA"

    # 重新加入 features（重算座標，但不破壞 join / complement）
    for feature in record.features:
        if (
            feature.location.start >= start_pos
            and feature.location.end <= end_pos
        ):
            new_location = shift_location(feature.location, start_pos)

            new_feature = SeqFeature(
                location=new_location,
                type=feature.type,
                qualifiers=feature.qualifiers,
            )
            new_record.features.append(new_feature)

    SeqIO.write(new_record, output_gbk, "genbank")


# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(
        description="Extract a gene region from a GenBank file, preserving join/complement."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input scaffold GenBank file",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output GenBank file",
    )
    parser.add_argument(
        "-s", "--start",
        required=True,
        help="Start gene ID (e.g. g4)",
    )
    parser.add_argument(
        "-e", "--end",
        required=True,
        help="End gene ID (e.g. g16)",
    )

    args = parser.parse_args()

    extract_genes_gbk(
        input_gbk=args.input,
        output_gbk=args.output,
        start_gene=args.start,
        end_gene=args.end,
    )

    print("Finish")


if __name__ == "__main__":
    main()
