#!/usr/bin/env python3

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple


@dataclass(frozen=True)
class GfaGraph:
    sequences: Dict[int, str]
    edges: Dict[int, List[int]]
    edge_overlap: Dict[Tuple[int, int], int]
    in_degree: Dict[int, int]
    out_degree: Dict[int, int]

    @property
    def nodes(self) -> List[int]:
        return sorted(self.sequences)

    def sources(self) -> List[int]:
        return [node for node in self.nodes if self.in_degree.get(node, 0) == 0]

    def sinks(self) -> List[int]:
        return [node for node in self.nodes if self.out_degree.get(node, 0) == 0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Trace deletion-supporting unitigs and source-to-sink paths from a unitig GFA."
        )
    )
    parser.add_argument("--unitig-gfa", required=True, type=Path, help="Unitig GFA file")
    parser.add_argument("--reference-fasta", required=True, type=Path, help="Reference haplotype FASTA")
    parser.add_argument("--deleted-fasta", required=True, type=Path, help="Deleted haplotype FASTA")
    parser.add_argument(
        "--max-paths",
        type=int,
        default=1000,
        help="Maximum number of source-to-sink paths to enumerate (default: 1000)",
    )
    parser.add_argument(
        "--junction-k",
        type=int,
        default=55,
        help="Length of deletion-specific anchors to scan for (default: 55)",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> str:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return "".join(line for line in lines if not line.startswith(">"))


def infer_edge_overlap(left: str, right: str) -> int:
    max_len = min(len(left), len(right))
    for overlap in range(max_len, 0, -1):
        if left[-overlap:] == right[:overlap]:
            return overlap
    return 0


def read_gfa(path: Path) -> GfaGraph:
    sequences: Dict[int, str] = {}
    edges: Dict[int, List[int]] = defaultdict(list)
    raw_links: List[Tuple[int, int]] = []
    in_degree: Dict[int, int] = defaultdict(int)
    out_degree: Dict[int, int] = defaultdict(int)

    for line in path.read_text().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        record_type = fields[0]
        if record_type == "S":
            sequences[int(fields[1])] = fields[2]
        elif record_type == "L":
            source = int(fields[1])
            target = int(fields[3])
            edges[source].append(target)
            raw_links.append((source, target))
            out_degree[source] += 1
            in_degree[target] += 1
            in_degree[source] += 0
            out_degree[target] += 0

    for node in sequences:
        edges[node] = list(edges.get(node, []))
        in_degree[node] += 0
        out_degree[node] += 0

    edge_overlap: Dict[Tuple[int, int], int] = {}
    for source, target in raw_links:
        edge_overlap[(source, target)] = infer_edge_overlap(
            sequences[source],
            sequences[target],
        )

    return GfaGraph(
        sequences=sequences,
        edges=dict(edges),
        edge_overlap=edge_overlap,
        in_degree=dict(in_degree),
        out_degree=dict(out_degree),
    )


def path_sequence(graph: GfaGraph, path: Sequence[int]) -> str:
    if not path:
        return ""

    sequence = graph.sequences[path[0]]
    for prev_node, node in zip(path, path[1:]):
        trim = graph.edge_overlap.get((prev_node, node), 0)
        suffix = graph.sequences[node][trim:] if trim else graph.sequences[node]
        sequence += suffix
    return sequence


def enumerate_paths(graph: GfaGraph, max_paths: int) -> List[List[int]]:
    paths: List[List[int]] = []
    for source in graph.sources():
        stack: List[Tuple[int, List[int]]] = [(source, [source])]
        while stack:
            node, path = stack.pop()
            if len(paths) >= max_paths:
                return paths
            next_nodes = graph.edges.get(node, [])
            if not next_nodes:
                paths.append(path)
                continue
            for next_node in reversed(next_nodes):
                if next_node in path:
                    continue
                stack.append((next_node, path + [next_node]))
    return paths


def deletion_specific_windows(reference: str, deleted: str, window: int) -> List[str]:
    if window <= 0:
        raise ValueError("window length must be positive")
    if len(deleted) < window:
        return []

    windows: List[str] = []
    seen = set()
    for start in range(len(deleted) - window + 1):
        chunk = deleted[start : start + window]
        if chunk in reference or chunk in seen:
            continue
        seen.add(chunk)
        windows.append(chunk)
    return windows


def find_nodes_with_chunks(graph: GfaGraph, chunks: Iterable[str]) -> Dict[str, List[int]]:
    hits: Dict[str, List[int]] = {}
    for chunk in chunks:
        matches = [node for node, seq in graph.sequences.items() if chunk in seq]
        if matches:
            hits[chunk] = matches
    return hits


def summarize_path(graph: GfaGraph, path: Sequence[int], reference: str, deleted: str) -> str:
    seq = path_sequence(graph, path)
    in_reference = seq in reference
    in_deleted = seq in deleted
    return (
        f"nodes={','.join(str(node) for node in path)} len={len(seq)} "
        f"in_deleted={in_deleted} in_reference={in_reference}"
    )


def main() -> int:
    args = parse_args()

    graph = read_gfa(args.unitig_gfa)
    reference = read_fasta(args.reference_fasta)
    deleted = read_fasta(args.deleted_fasta)

    paths = enumerate_paths(graph, args.max_paths)
    chunks = deletion_specific_windows(reference, deleted, args.junction_k)
    chunk_hits = find_nodes_with_chunks(graph, chunks)

    deleted_only_nodes = sorted(
        node
        for node, seq in graph.sequences.items()
        if seq in deleted and seq not in reference
    )
    deleted_paths = [path for path in paths if path_sequence(graph, path) in deleted]
    deleted_only_paths = [path for path in deleted_paths if path_sequence(graph, path) not in reference]

    distinct_overlaps = sorted(set(graph.edge_overlap.values()))
    print(f"graph_nodes={len(graph.nodes)} edge_overlaps={distinct_overlaps}")
    print(f"sources={graph.sources()}")
    print(f"sinks={graph.sinks()}")
    print(f"enumerated_paths={len(paths)}")
    print()

    print("deleted_only_unitigs:")
    if deleted_only_nodes:
        for node in deleted_only_nodes:
            print(f"  node {node}: {graph.sequences[node]}")
    else:
        print("  none")
    print()

    print("deletion_specific_anchor_hits:")
    if chunk_hits:
        for chunk, nodes in chunk_hits.items():
            print(f"  nodes {nodes}: {chunk}")
    else:
        print("  none")
    print()

    print("paths_wholly_in_deleted_haplotype:")
    if deleted_paths:
        for path in deleted_paths:
            print(f"  {summarize_path(graph, path, reference, deleted)}")
    else:
        print("  none")
    print()

    print("paths_unique_to_deleted_haplotype:")
    if deleted_only_paths:
        for path in deleted_only_paths:
            print(f"  {summarize_path(graph, path, reference, deleted)}")
    else:
        print("  none")
    print()

    print("adjacency_for_deleted_only_nodes:")
    if deleted_only_nodes:
        for node in deleted_only_nodes:
            preds = [pred for pred, targets in graph.edges.items() if node in targets]
            succs = graph.edges.get(node, [])
            print(f"  node {node}: preds={preds} succs={succs}")
    else:
        print("  none")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())