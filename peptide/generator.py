"""Rule-based peptide candidate generation. Does different things to the seed motifs to get a bunch of candidates, then applies filters and deduplicates."""

from __future__ import annotations
import itertools
import random
from collections import OrderedDict
from typing import Iterable

from .config import (
    ALLOWED_FLANKING_RESIDUES,
    ALLOWED_SPACER_RESIDUES,
    MAX_TOTAL_CANDIDATES,
    MAX_VARIANTS_PER_MOTIF,
    PEPTIDE_MAX_LENGTH,
    PEPTIDE_MIN_LENGTH,
    RANDOM_SEED,
)
from .filters import apply_all_filters
from .metadata import compute_sequence_metadata
from .schema import PeptideCandidate, SeedMotif


def assign_candidate_id(index: int) -> str:
    return f"pep_{index:05d}"


def _unique_in_order(values: Iterable[str]) -> list[str]:
    return list(OrderedDict.fromkeys(values))


def _pad_to_min_length(sequence: str) -> str:
    if len(sequence) >= PEPTIDE_MIN_LENGTH:
        return sequence

    needed = PEPTIDE_MIN_LENGTH - len(sequence)
    left: list[str] = []
    right: list[str] = []
    for index in range(needed):
        residue = ALLOWED_FLANKING_RESIDUES[index % len(ALLOWED_FLANKING_RESIDUES)]
        if index % 2 == 0:
            left.append(residue)
        else:
            right.append(residue)
    return "".join(reversed(left)) + sequence + "".join(right)


def expand_pattern_motif(
    motif: SeedMotif,
    max_variants: int = MAX_VARIANTS_PER_MOTIF,
) -> list[str]:
    if motif.pattern_type == "exact":
        return [motif.normalized_pattern]

    choices: list[tuple[str, ...]] = []
    for symbol in motif.normalized_pattern:
        if symbol == "B":
            choices.append(("K", "R"))
        elif symbol == "X":
            choices.append(ALLOWED_SPACER_RESIDUES)
        else:
            choices.append((symbol,))

    variants: list[str] = []
    for product in itertools.product(*choices):
        variants.append("".join(product))
        if len(variants) >= max_variants:
            break
    return variants


def build_flanked_variants(
    core_sequence: str,
    rng: random.Random,
    max_variants: int = MAX_VARIANTS_PER_MOTIF,
) -> list[tuple[str, str]]:
    base = _pad_to_min_length(core_sequence)
    variants: list[tuple[str, str]] = [(base, "core")]

    flank_templates = [
        ("", "G"),
        ("G", ""),
        ("", "Q"),
        ("Q", ""),
        ("G", "G"),
        ("Q", "Q"),
        ("S", ""),
        ("", "S"),
        ("A", ""),
        ("", "A"),
        ("G", "Q"),
        ("Q", "G"),
    ]

    ordered_templates = list(flank_templates)
    rng.shuffle(ordered_templates)

    for left_flank, right_flank in ordered_templates:
        candidate = f"{left_flank}{base}{right_flank}"
        if len(candidate) > PEPTIDE_MAX_LENGTH:
            continue
        variants.append((candidate, "flanked"))
        if len(variants) >= max_variants:
            break

    unique_variants: list[tuple[str, str]] = []
    seen_sequences: set[str] = set()
    for sequence, method in variants:
        if sequence in seen_sequences:
            continue
        seen_sequences.add(sequence)
        unique_variants.append((sequence, method))
    return unique_variants[:max_variants]


def _build_candidate(
    sequence: str,
    motif: SeedMotif,
    generation_method: str,
) -> PeptideCandidate:
    metadata = compute_sequence_metadata(sequence)
    filter_summary = apply_all_filters(sequence, metadata)
    return PeptideCandidate(
        candidate_id="",
        sequence=sequence,
        length=int(metadata["length"]),
        source_motif_ids=[motif.motif_id],
        generation_method=generation_method,
        parent_motif_pattern=motif.normalized_pattern,
        approx_net_charge=metadata["approx_net_charge"],
        residue_counts=dict(metadata["residue_counts"]),
        residue_fractions=dict(metadata["residue_fractions"]),
        basic_fraction=float(metadata["basic_fraction"]),
        hydrophobic_fraction=float(metadata["hydrophobic_fraction"]),
        polar_fraction=float(metadata["polar_fraction"]),
        acidic_fraction=float(metadata["acidic_fraction"]),
        longest_basic_run=int(metadata["longest_basic_run"]),
        longest_hydrophobic_run=int(metadata["longest_hydrophobic_run"]),
        basic_positions=list(metadata["basic_positions"]),
        charge_spacing=list(metadata["charge_spacing"]),
        motif_hits=[motif.normalized_pattern],
        passed_filters=bool(filter_summary["passed_filters"]),
        filter_flags=list(filter_summary["filter_flags"]),
        warning_flags=list(filter_summary["warning_flags"]),
    )


def deduplicate_by_sequence(candidates: Iterable[PeptideCandidate]) -> list[PeptideCandidate]:
    deduped: OrderedDict[str, PeptideCandidate] = OrderedDict()

    for candidate in candidates:
        existing = deduped.get(candidate.sequence)
        if existing is None:
            deduped[candidate.sequence] = candidate
            continue

        existing.source_motif_ids = _unique_in_order(
            existing.source_motif_ids + candidate.source_motif_ids
        )
        existing.motif_hits = _unique_in_order(existing.motif_hits + candidate.motif_hits)
        existing.filter_flags = _unique_in_order(existing.filter_flags + candidate.filter_flags)
        existing.warning_flags = _unique_in_order(existing.warning_flags + candidate.warning_flags)
        existing.passed_filters = existing.passed_filters and candidate.passed_filters

    final_candidates = list(deduped.values())
    for index, candidate in enumerate(final_candidates, start=1):
        candidate.candidate_id = assign_candidate_id(index)
    return final_candidates


def generate_candidates_from_motifs(
    motifs: list[SeedMotif],
    *,
    max_variants_per_motif: int = MAX_VARIANTS_PER_MOTIF,
    max_total_candidates: int = MAX_TOTAL_CANDIDATES,
    random_seed: int = RANDOM_SEED,
) -> list[PeptideCandidate]:
    rng = random.Random(random_seed)
    generated: list[PeptideCandidate] = []

    for motif in motifs:
        core_limit = max(1, max_variants_per_motif // 3)
        core_sequences = expand_pattern_motif(motif, max_variants=core_limit)
        motif_sequences: set[str] = set()

        for core_sequence in core_sequences:
            remaining_for_motif = max_variants_per_motif - len(motif_sequences)
            if remaining_for_motif <= 0:
                break

            variant_method = "exact" if motif.pattern_type == "exact" else "pattern_expansion"
            flanked_variants = build_flanked_variants(
                core_sequence,
                rng=rng,
                max_variants=remaining_for_motif,
            )

            for sequence, flank_method in flanked_variants:
                if sequence in motif_sequences:
                    continue
                motif_sequences.add(sequence)
                generated.append(
                    _build_candidate(
                        sequence=sequence,
                        motif=motif,
                        generation_method=f"{variant_method}_{flank_method}",
                    )
                )
                if len(generated) >= max_total_candidates:
                    return deduplicate_by_sequence(generated)
                if len(motif_sequences) >= max_variants_per_motif:
                    break

    return deduplicate_by_sequence(generated)
