import json
import pysam
import sys

def get_rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def get_canonical(kmer):
    rc = get_rev_comp(kmer)
    return (kmer, False) if kmer <= rc else (rc, True)

def main():
    k = 45
    json_path = "resources/example_resources/eg3/results/sharda_eg3_k45_debug/raw.json"
    bam_path = "resources/example_resources/eg3/results/example_reads.raw.bam"
    target_qnames = {
        "hap1_hap1_del900_121_587_0:0:0_0:0:0_c",
        "hap2_hap2_del30_183_632_0:0:0_0:0:0_d"
    }

    with open(json_path) as f:
        data = json.load(f)

    # backbone_lookup: kmer -> list of (node_id, ref_pos)
    # read_node_lookup: kmer -> node_id
    backbone_lookup = {}
    read_node_lookup = {}

    for node in data["nodes"]:
        nid = node["id"]
        seq = node["sequence"]
        is_backbone = node["is_backbone"]
        
        if is_backbone:
            # Backbone nodes are typically length k
            # Actually, in sharda raw graph, backbone nodes might be longer if collapsed, 
            # but usually they are per-kmer in initial construction or have ref_pos.
            # Let's check keys.
            ref_pos = node.get("ref_pos", -1)
            kmer, _ = get_canonical(seq)
            backbone_lookup.setdefault(kmer, []).append((nid, ref_pos))
        else:
            kmer, _ = get_canonical(seq)
            read_node_lookup[kmer] = nid

    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam:
        if read.query_name not in target_qnames:
            continue
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue

        seq = read.query_sequence
        ref_start = read.reference_start # 0-based POS-1
        
        path = []
        missing = []
        
        for i in range(len(seq) - k + 1):
            kmer_raw = seq[i:i+k]
            kmer, _ = get_canonical(kmer_raw)
            implied_ref_pos = ref_start + i
            
            chosen_node = None
            if kmer in backbone_lookup:
                # Closest backbone node
                candidates = backbone_lookup[kmer]
                # min abs dist, then min ref_pos
                best_node = None
                min_dist = float('inf')
                for nid, rpos in candidates:
                    dist = abs(rpos - implied_ref_pos)
                    if dist < min_dist:
                        min_dist = dist
                        best_node = nid
                    elif dist == min_dist:
                        if best_node is None or rpos < next(node_rpos for node_id, node_rpos in candidates if node_id == best_node):
                             best_node = nid
                chosen_node = (best_node, "B")
            elif kmer in read_node_lookup:
                chosen_node = (read_node_lookup[kmer], "R")
            else:
                chosen_node = (None, "?")
                missing.append(i)
            
            path.append(chosen_node)

        # Summary
        runs = []
        if path:
            curr_id, curr_type = path[0]
            start_idx = 0
            for i in range(1, len(path)):
                nid, ntype = path[i]
                if nid != curr_id or ntype != curr_type:
                    runs.append(f"{curr_id}{curr_type}[{i-start_idx}]")
                    curr_id, curr_type = nid, ntype
                    start_idx = i
            runs.append(f"{curr_id}{curr_type}[{len(path)-start_idx}]")

        non_backbone_nodes = sorted(list(set(nid for nid, ntype in path if ntype == "R")))

        print(f"READ: {read.query_name}")
        print(f"FLAG: {read.flag} POS: {read.reference_start + 1} CIGAR: {read.cigarstring}")
        print(f"KMERS: {len(path)} MISSING: {len(missing)}")
        print(f"PATH: {' '.join(runs)}")
        print(f"READ_NODES: {non_backbone_nodes}")
        print("-" * 20)

main()
