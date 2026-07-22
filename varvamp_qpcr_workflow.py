#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ayman Faham
"""
Generic workflow for qPCR primer and probe design from nucleotide sequences.

Steps:
1. Check the FASTA input file and required software.
2. Reduce sequence redundancy with CD-HIT-EST.
3. Perform multiple sequence alignment with MAFFT.
4. Analyse alignment conservation.
5. Design qPCR primers and probes with VarVAMP.
6. Copy and visualise the main results.

Interactive execution from the repository root:

    python varvamp_qpcr_workflow_generic.py

The program then asks for the name or path of the input FASTA file.

Non-interactive execution:

    python varvamp_qpcr_workflow_generic.py \
        --input path/to/sequences.fasta \
        --workdir work \
        --results results

External dependencies:
    cd-hit-est
    mafft
    varvamp

Python libraries:
    biopython
    pandas
    matplotlib
    pillow
    pymupdf
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from collections import Counter
from io import BytesIO
from math import log2
from pathlib import Path
from typing import Iterable

import matplotlib

# Allow execution on a server or in GitHub Codespaces without a display.
matplotlib.use("Agg")

import fitz
import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO
from PIL import Image


VALID_BASES = {"A", "C", "G", "T"}
GAP_CHARACTERS = {"-", "."}


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Sequence dereplication, multiple alignment, conservation analysis "
            "and VarVAMP qPCR design for nucleotide sequences."
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help=(
            "Input FASTA file. If omitted, "
            "the program asks for the path at startup."
        ),
    )
    parser.add_argument(
        "--project-name",
        type=str,
        default=None,
        help=(
            "Project name. If omitted, the FASTA filename "
            "is used automatically."
        ),
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help=(
            "Directory for intermediate files. Default: "
            "work/<project_name>."
        ),
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=None,
        help=(
            "Directory for final results. Default: "
            "results/<project_name>."
        ),
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=0.95,
        help="CD-HIT-EST identity threshold (default: 0.95).",
    )
    parser.add_argument(
        "--min-occupancy",
        type=float,
        default=0.95,
        help="Minimum occupancy threshold for conservation (default: 0.95).",
    )
    parser.add_argument(
        "--min-major-frequency",
        type=float,
        default=0.95,
        help="Minimum major-base frequency (default: 0.95).",
    )
    parser.add_argument(
        "--varvamp-threshold",
        type=float,
        default=0.95,
        help="VarVAMP consensus threshold (default: 0.95).",
    )
    parser.add_argument(
        "--primer-ambiguity",
        type=int,
        default=2,
        help="Maximum number of ambiguous bases allowed in primers.",
    )
    parser.add_argument(
        "--probe-ambiguity",
        type=int,
        default=2,
        help="Maximum number of ambiguous bases allowed in the probe.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads used by VarVAMP (default: 8).",
    )
    parser.add_argument(
        "--skip-varvamp",
        action="store_true",
        help="Run the preceding steps without launching VarVAMP.",
    )

    return parser.parse_args()


def ask_for_input_file() -> Path:
    """Interactively request the name or path of the FASTA file."""

    while True:
        try:
            raw_value = input(
                "Name or path of the input FASTA file: "
            ).strip()
        except EOFError as error:
            raise RuntimeError(
                "No interactive input is available. "
                "Use the --input option."
            ) from error

        # Allow quotation marks copied from the terminal.
        raw_value = raw_value.strip("'\"")

        if not raw_value:
            print("Please enter a filename or file path.")
            continue

        candidate = Path(raw_value).expanduser()

        # If only a filename is provided, search the current directory first,
        # then the data/ subdirectory.
        possible_paths = [candidate]

        if not candidate.is_absolute() and candidate.parent == Path("."):
            possible_paths.append(Path("data") / candidate)

        for possible_path in possible_paths:
            resolved_path = possible_path.resolve()

            if resolved_path.is_file():
                return resolved_path

        print(
            "File not found. Checked paths: "
            + ", ".join(str(path.resolve()) for path in possible_paths)
        )


def sanitize_project_name(value: str) -> str:
    """Convert a name into a safe identifier for files and directories."""

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    cleaned = cleaned.strip("._-")

    if not cleaned:
        raise ValueError("The project name is empty or invalid.")

    return cleaned


def validate_nucleotide_fasta(file_path: Path) -> None:
    """Check that the FASTA file contains nucleotide sequences."""

    allowed_characters = set("ACGTURYSWKMBDHVN-.")

    sequence_count = 0
    sequence_characters = 0
    invalid_characters: set[str] = set()

    with file_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith(">"):
                sequence_count += 1
                continue

            upper_line = stripped.upper().replace(" ", "")
            sequence_characters += len(upper_line)
            invalid_characters.update(
                character
                for character in upper_line
                if character not in allowed_characters
            )

    if sequence_count == 0:
        raise ValueError("The file does not contain any FASTA headers.")

    if sequence_characters == 0:
        raise ValueError("The FASTA file does not contain any sequences.")

    if invalid_characters:
        invalid_display = ", ".join(sorted(invalid_characters))
        raise ValueError(
            "The file appears to contain non-nucleotide characters: "
            f"{invalid_display}. This workflow uses CD-HIT-EST and expects "
            "DNA or RNA sequences."
        )


def validate_fraction(value: float, name: str) -> None:
    """Check that a value is between 0 and 1."""

    if not 0 < value <= 1:
        raise ValueError(f"{name} must be between 0 and 1.")


def check_required_tools(tools: Iterable[str]) -> None:
    """Check that required external programs are available."""

    missing = [tool for tool in tools if shutil.which(tool) is None]

    if missing:
        raise RuntimeError(
            "Program(s) not found in PATH: "
            + ", ".join(missing)
            + ". Install them in the active environment."
        )


def run_command(
    command: list[str],
    *,
    stdout_file: Path | None = None,
) -> None:
    """Run an external command and stop the workflow if it fails."""

    print("\nCommand:", " ".join(command))

    try:
        if stdout_file is None:
            subprocess.run(command, check=True)
        else:
            stdout_file.parent.mkdir(parents=True, exist_ok=True)

            with stdout_file.open("w", encoding="utf-8") as handle:
                subprocess.run(
                    command,
                    check=True,
                    stdout=handle,
                )

    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"The command failed with exit code {error.returncode}: "
            + " ".join(command)
        ) from error


def count_fasta_sequences(file_path: Path) -> int:
    """Count FASTA headers without loading the entire file into memory."""

    if not file_path.is_file():
        raise FileNotFoundError(f"File not found: {file_path}")

    with file_path.open("r", encoding="utf-8") as handle:
        return sum(line.startswith(">") for line in handle)


def read_alignment(file_path: Path):
    """Read and validate a FASTA alignment."""

    if not file_path.is_file():
        raise FileNotFoundError(
            f"Alignment file not found: {file_path.resolve()}"
        )

    alignment = AlignIO.read(file_path, "fasta")

    if len(alignment) == 0:
        raise ValueError("The alignment contains no sequences.")

    return alignment


def calculate_real_lengths(alignment) -> list[int]:
    """Calculate sequence lengths after removing gaps."""

    real_lengths: list[int] = []

    for record in alignment:
        sequence = str(record.seq).upper()
        sequence_without_gaps = "".join(
            character
            for character in sequence
            if character not in GAP_CHARACTERS
        )
        real_lengths.append(len(sequence_without_gaps))

    return real_lengths


def calculate_shannon_entropy(base_counts: Counter[str]) -> float:
    """Calculate Shannon entropy for A, C, G and T bases."""

    total_valid_bases = sum(base_counts.values())

    if total_valid_bases == 0:
        return 0.0

    entropy = 0.0

    for count in base_counts.values():
        if count > 0:
            frequency = count / total_valid_bases
            entropy -= frequency * log2(frequency)

    return entropy


def analyse_alignment_positions(
    alignment,
    min_occupancy: float,
    min_major_frequency: float,
) -> tuple[pd.DataFrame, str]:
    """Analyse each column of the multiple sequence alignment."""

    number_of_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    results: list[dict[str, object]] = []
    consensus_bases: list[str] = []

    for position_index in range(alignment_length):
        column = [
            str(record.seq[position_index]).upper()
            for record in alignment
        ]

        gap_count = sum(
            character in GAP_CHARACTERS
            for character in column
        )

        non_gap_characters = [
            character
            for character in column
            if character not in GAP_CHARACTERS
        ]

        occupied_count = len(non_gap_characters)
        occupancy = occupied_count / number_of_sequences

        valid_bases = [
            character
            for character in non_gap_characters
            if character in VALID_BASES
        ]

        ambiguous_count = sum(
            character not in VALID_BASES
            for character in non_gap_characters
        )

        base_counts = Counter(valid_bases)

        if base_counts:
            major_base, major_base_count = base_counts.most_common(1)[0]
        else:
            major_base = "N"
            major_base_count = 0

        major_base_frequency = (
            major_base_count / occupied_count
            if occupied_count > 0
            else 0.0
        )

        sufficiently_represented = occupancy >= min_occupancy
        strictly_conserved = (
            gap_count == 0
            and ambiguous_count == 0
            and major_base_count == number_of_sequences
        )
        highly_conserved = (
            sufficiently_represented
            and major_base_frequency >= min_major_frequency
        )

        consensus_base = major_base if highly_conserved else "N"
        consensus_bases.append(consensus_base)

        results.append(
            {
                "position_alignment": position_index + 1,
                "A_count": base_counts.get("A", 0),
                "C_count": base_counts.get("C", 0),
                "G_count": base_counts.get("G", 0),
                "T_count": base_counts.get("T", 0),
                "valid_base_count": len(valid_bases),
                "gap_count": gap_count,
                "ambiguous_count": ambiguous_count,
                "occupied_count": occupied_count,
                "occupancy": occupancy,
                "major_base": major_base,
                "major_base_count": major_base_count,
                "major_base_frequency": major_base_frequency,
                "shannon_entropy": calculate_shannon_entropy(base_counts),
                "sufficiently_represented": sufficiently_represented,
                "strictly_conserved": strictly_conserved,
                "highly_conserved": highly_conserved,
                "consensus_base": consensus_base,
            }
        )

    return pd.DataFrame(results), "".join(consensus_bases)


def save_consensus_fasta(
    consensus_sequence: str,
    output_file: Path,
    project_name: str,
) -> None:
    """Save the consensus sequence in FASTA format."""

    with output_file.open("w", encoding="utf-8") as handle:
        handle.write(
            f">{project_name}_consensus_occupancy95_major95\n"
        )

        for start in range(0, len(consensus_sequence), 80):
            handle.write(consensus_sequence[start : start + 80] + "\n")


def perform_conservation_analysis(
    alignment_file: Path,
    results_dir: Path,
    project_name: str,
    min_occupancy: float,
    min_major_frequency: float,
) -> None:
    """Generate the conservation table, consensus and summary."""

    alignment = read_alignment(alignment_file)
    number_of_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    real_lengths = calculate_real_lengths(alignment)
    minimum_length = min(real_lengths)
    maximum_length = max(real_lengths)
    average_length = sum(real_lengths) / len(real_lengths)

    results_df, consensus_sequence = analyse_alignment_positions(
        alignment,
        min_occupancy,
        min_major_frequency,
    )

    total_alignment_characters = number_of_sequences * alignment_length
    total_gap_characters = int(results_df["gap_count"].sum())
    total_ambiguous_characters = int(
        results_df["ambiguous_count"].sum()
    )

    gap_percentage = (
        100 * total_gap_characters / total_alignment_characters
    )
    ambiguous_percentage = (
        100 * total_ambiguous_characters / total_alignment_characters
    )

    represented_count = int(
        results_df["sufficiently_represented"].sum()
    )
    strict_count = int(results_df["strictly_conserved"].sum())
    conserved_count = int(results_df["highly_conserved"].sum())

    represented_percentage = 100 * represented_count / alignment_length
    strict_percentage = 100 * strict_count / alignment_length
    conserved_total_percentage = (
        100 * conserved_count / alignment_length
    )

    conserved_represented_percentage = (
        100 * conserved_count / represented_count
        if represented_count > 0
        else 0.0
    )

    results_dir.mkdir(parents=True, exist_ok=True)

    position_file = (
        results_dir / f"{project_name}_conservation_by_position.csv"
    )
    consensus_file = results_dir / f"{project_name}_consensus.fasta"
    summary_file = (
        results_dir / f"{project_name}_conservation_summary.txt"
    )

    results_df.to_csv(position_file, index=False)
    save_consensus_fasta(
        consensus_sequence,
        consensus_file,
        project_name,
    )

    summary_lines = [
        "ALIGNMENT CONSERVATION ANALYSIS",
        "=" * 55,
        "",
        "ANALYSED DATA",
        f"Sequences analysed: {number_of_sequences}",
        f"Alignment length: {alignment_length} positions",
        (
            "Ungapped sequence length: "
            f"{minimum_length}-{maximum_length} nt"
        ),
        f"Mean ungapped sequence length: {average_length:.2f} nt",
        "",
        "ALIGNMENT QUALITY",
        (
            f"Gaps: {total_gap_characters} / "
            f"{total_alignment_characters} ({gap_percentage:.2f} %)"
        ),
        (
            f"Ambiguous characters: {total_ambiguous_characters} / "
            f"{total_alignment_characters} "
            f"({ambiguous_percentage:.4f} %)"
        ),
        "",
        "SUFFICIENTLY REPRESENTED POSITIONS",
        (
            f"Occupancy >= {min_occupancy:.0%} : "
            f"{represented_count} / {alignment_length} "
            f"({represented_percentage:.2f} %)"
        ),
        "",
        "CONSERVATION",
        (
            f"Strictly conserved positions: "
            f"{strict_count} / {alignment_length} "
            f"({strict_percentage:.2f} %)"
        ),
        (
            f"Positions with occupancy >= {min_occupancy:.0%} "
            f"and major-base frequency >= {min_major_frequency:.0%} : "
            f"{conserved_count} / {alignment_length} "
            f"({conserved_total_percentage:.2f} %)"
        ),
        (
            "Conservation among sufficiently "
            f"represented positions: {conserved_count} / {represented_count} "
            f"({conserved_represented_percentage:.2f} %)"
        ),
    ]

    summary_text = "\n".join(summary_lines)
    summary_file.write_text(summary_text, encoding="utf-8")

    print("\n" + summary_text)
    print("\nConservation files generated:")
    print(f"- {position_file.resolve()}")
    print(f"- {consensus_file.resolve()}")
    print(f"- {summary_file.resolve()}")


def copy_varvamp_results(
    varvamp_dir: Path,
    results_dir: Path,
) -> None:
    """Copy the main VarVAMP files into the results directory."""

    expected_files = [
        "qpcr_primers.tsv",
        "qpcr_design.tsv",
        "primers.bed",
        "amplicons.bed",
        "amplicon_plot.pdf",
        "per_base_mismatches.pdf",
        "ambiguous_consensus.fasta",
        "varvamp_log.txt",
    ]

    for filename in expected_files:
        source = varvamp_dir / filename

        if source.is_file():
            shutil.copy2(source, results_dir / filename)
        else:
            print(f"Warning: missing VarVAMP result: {source}")


def plot_bed_segments(
    bed_file: Path,
    output_file: Path,
    title: str,
    y_label: str,
    line_width: float,
) -> None:
    """Create a simple plot from a BED file."""

    if not bed_file.is_file():
        print(f"Warning: BED file not found: {bed_file}")
        return

    data = pd.read_csv(bed_file, sep="\t", header=None)

    if data.empty:
        print(f"Warning: BED file is empty: {bed_file}")
        return

    plt.figure(figsize=(12, 4))

    for index, row in data.iterrows():
        plt.plot(
            [row.iloc[1], row.iloc[2]],
            [index, index],
            linewidth=line_width,
        )

    plt.xlabel("Alignment position")
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Plot saved: {output_file.resolve()}")


def save_pdf_first_page_as_png(
    pdf_file: Path,
    output_file: Path,
    *,
    zoom: int = 4,
) -> None:
    """Convert the first page of a PDF file to PNG."""

    if not pdf_file.is_file():
        print(f"Warning: PDF file not found: {pdf_file}")
        return

    with fitz.open(pdf_file) as pdf:
        if len(pdf) == 0:
            raise ValueError(f"The PDF file is empty: {pdf_file}")

        page = pdf[0]
        pixmap = page.get_pixmap(
            matrix=fitz.Matrix(zoom, zoom),
            alpha=False,
        )

    image = Image.open(BytesIO(pixmap.tobytes("png"))).copy()
    image.save(output_file)

    print(f"First page converted: {output_file.resolve()}")


def print_result_tables(results_dir: Path) -> None:
    """Print qPCR result tables in the terminal."""

    for filename, label in [
        ("qpcr_primers.tsv", "PRIMERS AND PROBE"),
        ("qpcr_design.tsv", "qPCR DESIGNS"),
    ]:
        file_path = results_dir / filename

        if not file_path.is_file():
            continue

        table = pd.read_csv(file_path, sep="\t")

        print(f"\n{label}")
        print("=" * len(label))
        print(table.to_string(index=False))


def main() -> int:
    """Run the complete workflow."""

    args = parse_arguments()

    try:
        validate_fraction(args.identity, "--identity")
        validate_fraction(args.min_occupancy, "--min-occupancy")
        validate_fraction(
            args.min_major_frequency,
            "--min-major-frequency",
        )
        validate_fraction(
            args.varvamp_threshold,
            "--varvamp-threshold",
        )

        if args.input is None:
            input_file = ask_for_input_file()
        else:
            input_file = args.input.expanduser().resolve()

            if not input_file.is_file():
                raise FileNotFoundError(
                    "The input FASTA file was not found: "
                    f"{input_file}"
                )

        validate_nucleotide_fasta(input_file)

        project_name = sanitize_project_name(
            args.project_name if args.project_name else input_file.stem
        )

        workdir = (
            args.workdir.expanduser().resolve()
            if args.workdir is not None
            else (Path("work") / project_name).resolve()
        )
        results_dir = (
            args.results.expanduser().resolve()
            if args.results is not None
            else (Path("results") / project_name).resolve()
        )

        print(f"Selected FASTA file: {input_file}")
        print(f"Project name: {project_name}")
        print(f"Working directory: {workdir}")
        print(f"Results directory: {results_dir}")

        required_tools = ["cd-hit-est", "mafft"]

        if not args.skip_varvamp:
            required_tools.append("varvamp")

        check_required_tools(required_tools)

        workdir.mkdir(parents=True, exist_ok=True)
        results_dir.mkdir(parents=True, exist_ok=True)

        non_redundant_file = workdir / f"{project_name}_nr.fasta"
        alignment_file = workdir / f"{project_name}_alignment.fasta"
        varvamp_dir = workdir / f"{project_name}_qpcr"

        print("\n1. Sequence dereplication with CD-HIT-EST")
        run_command(
            [
                "cd-hit-est",
                "-i",
                str(input_file),
                "-o",
                str(non_redundant_file),
                "-c",
                str(args.identity),
                "-d",
                "0",
            ]
        )

        initial_count = count_fasta_sequences(input_file)
        non_redundant_count = count_fasta_sequences(non_redundant_file)
        reduction = (
            100 * (initial_count - non_redundant_count) / initial_count
            if initial_count > 0
            else 0.0
        )

        print(f"Initial sequences: {initial_count}")
        print(f"Non-redundant sequences: {non_redundant_count}")
        print(f"Reduction: {reduction:.2f} %")

        print("\n2. Multiple sequence alignment with MAFFT")
        run_command(
            [
                "mafft",
                "--auto",
                str(non_redundant_file),
            ],
            stdout_file=alignment_file,
        )

        print("\n3. Conservation analysis")
        perform_conservation_analysis(
            alignment_file,
            results_dir,
            project_name,
            args.min_occupancy,
            args.min_major_frequency,
        )

        # Keep the final alignment among the reproducible results.
        shutil.copy2(
            alignment_file,
            results_dir / f"{project_name}_alignment.fasta",
        )

        if not args.skip_varvamp:
            print("\n4. qPCR design with VarVAMP")

            if varvamp_dir.exists():
                shutil.rmtree(varvamp_dir)

            run_command(
                [
                    "varvamp",
                    "qpcr",
                    "-t",
                    str(args.varvamp_threshold),
                    "-a",
                    str(args.primer_ambiguity),
                    "-pa",
                    str(args.probe_ambiguity),
                    "-th",
                    str(args.threads),
                    str(alignment_file),
                    str(varvamp_dir),
                ]
            )

            copy_varvamp_results(varvamp_dir, results_dir)
            print_result_tables(results_dir)

            plot_bed_segments(
                results_dir / "amplicons.bed",
                results_dir / "amplicons_overview.png",
                "Amplicon positions",
                "Amplicon",
                6,
            )
            plot_bed_segments(
                results_dir / "primers.bed",
                results_dir / "primers_overview.png",
                "Primer and probe positions",
                "Oligonucleotide",
                4,
            )

            save_pdf_first_page_as_png(
                results_dir / "amplicon_plot.pdf",
                results_dir / "amplicon_plot.png",
            )
            save_pdf_first_page_as_png(
                results_dir / "per_base_mismatches.pdf",
                results_dir / "per_base_mismatches.png",
            )

        print("\nWorkflow completed successfully.")
        print(f"Results: {results_dir}")

        return 0

    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        pd.errors.ParserError,
    ) as error:
        print(f"\nERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
