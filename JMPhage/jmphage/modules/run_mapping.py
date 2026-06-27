#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
from pathlib import Path
import logging
from datetime import datetime
import traceback
import shutil

ENV_NAME = "JM_mapping"
PHAGETERM_ENV = "jmp_Phageterm"

PHAGETERM_SCRIPT_REL = "mapping_profile/ptv-py3_release_1_light/PhageTerm.py"


def setup_logger(output_dir, module_name="mapping"):

    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{module_name}_{timestamp}.log"

    logger = logging.getLogger(module_name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(file_formatter)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger, log_file


def run_command(cmd, description="", logger=None):

    if logger:
        logger.info(f"\n{'='*60}\n{description}\n{'='*60}")

    if isinstance(cmd, str):
        import shlex
        cmd = shlex.split(cmd)

    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        if logger:
            for line in process.stdout:
                line_clean = line.rstrip()
                if line_clean:
                    logger.debug(line_clean)

        process.wait()

        if process.returncode != 0:
            if logger:
                logger.error(f"[FAILED] Command exited with return code {process.returncode}")
            raise subprocess.CalledProcessError(process.returncode, cmd)

        return process.returncode

    except subprocess.CalledProcessError as e:
        if logger:
            logger.error(f"Command execution failed.")
        raise


def run_command_in_env(cmd, env_name, description="", logger=None, cwd=None):

    if logger:
        logger.info(f"\n{'='*60}\n{description}\n{'='*60}")

    if isinstance(cmd, str):
        cmd = cmd.split()

    full_cmd = ["conda", "run", "--no-capture-output", "-n", env_name] + cmd

    if logger:
        logger.debug(f"Executing: {' '.join(full_cmd)}")

    try:
        process = subprocess.Popen(
            full_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=cwd,
            bufsize=1
        )

        for line in process.stdout:
            line_clean = line.rstrip()
            if line_clean and logger:
                logger.debug(f"[PhageTerm] {line_clean}")

        process.wait()

        if process.returncode != 0:
            if logger:
                logger.warning(f"[PhageTerm] exited with code {process.returncode}")
            return process.returncode

        return 0

    except FileNotFoundError:
        if logger:
            logger.warning("[PhageTerm] conda command not found, skipping PhageTerm.")
        return 1
    except Exception as e:
        if logger:
            logger.warning(f"[PhageTerm] Unexpected error: {e}")
        return 1


def run_phageterm(input_fasta, read1, read2, sample_name, mapping_dir, database, threads, logger):

    logger.info(f"\n{'#'*60}\n# PhageTerm: Packaging mechanism & terminus detection\n{'#'*60}\n")

    phageterm_dir = mapping_dir / "phageterm"

    marker_fasta = mapping_dir / f"{sample_name}_phageterm.fasta"
    report_pdf = phageterm_dir / f"{sample_name}_PhageTerm_report.pdf"
    stats_csv = phageterm_dir / f"{sample_name}_statistics.csv"

    if marker_fasta.exists():
        logger.info(f"[SKIP] PhageTerm output already exists: {marker_fasta.name}")
        return marker_fasta

    if report_pdf.exists() or stats_csv.exists():
        logger.info("[SKIP] PhageTerm was previously run but produced no reoriented sequence.")
        return None

    phageterm_dir.mkdir(parents=True, exist_ok=True)

    phageterm_script = Path(database) / PHAGETERM_SCRIPT_REL
    if not phageterm_script.exists():
        logger.warning(f"[SKIP] PhageTerm script not found: {phageterm_script}")
        return None

    input_fasta_abs = str(Path(input_fasta).resolve())
    read1_abs = str(Path(read1).resolve())
    read2_abs = str(Path(read2).resolve())

    cmd = [
        "python", str(phageterm_script),
        "-f", read1_abs,
        "-p", read2_abs,
        "-r", input_fasta_abs,
        "--report_title", sample_name,
        "-c", str(threads)
    ]

    rc = run_command_in_env(
        cmd, PHAGETERM_ENV,
        description=f"Running PhageTerm for {sample_name}...",
        logger=logger,
        cwd=str(phageterm_dir)
    )

    if rc != 0:
        logger.warning(f"[WARN] PhageTerm exited with code {rc}. "
                        "This is non-fatal; continuing with original FASTA.")
        return None

    seq_fasta = phageterm_dir / f"{sample_name}_sequence.fasta"

    if not seq_fasta.exists():

        candidates = list(phageterm_dir.glob("*_sequence.fasta"))
        if candidates:
            seq_fasta = candidates[0]
            logger.info(f"[PhageTerm] Found reoriented sequence: {seq_fasta.name}")
        else:
            logger.info("[PhageTerm] No reoriented sequence produced "
                        "(phage may lack identifiable packaging signal). "
                        "Continuing with original FASTA.")
            return None

    valid = False
    with open(seq_fasta, 'r') as f:
        for line in f:
            if not line.startswith('>') and line.strip():
                valid = True
                break

    if not valid:
        logger.warning("[PhageTerm] _sequence.fasta exists but contains no sequence data. "
                        "This phage may use headful packaging (no fixed terminus). "
                        "Ignoring and continuing with original FASTA.")
        return None

    shutil.copy(seq_fasta, marker_fasta)
    logger.info(f"[PhageTerm] Reoriented sequence saved: {marker_fasta}")
    logger.info(f"[PhageTerm] Downstream modules will use this sequence instead of the original.")

    for pattern in ["*_statistics.csv"]:
        for f in phageterm_dir.glob(pattern):
            f.unlink(missing_ok=True)
            logger.info(f"[PhageTerm] Cleaned intermediate file: {f.name}")
            
    return marker_fasta




def run_mapping(input_fasta, read1, read2, threads, database, output_dir):
    output_dir = Path(output_dir).resolve()
    mapping_dir = output_dir / "1.read_mapping"
    sample_name = Path(input_fasta).stem

    if mapping_dir.exists():
        plot_files = list(mapping_dir.glob("*.pdf"))
        if plot_files:
            print(f"Detected existing result plot: {plot_files[0].name}")
            print(f"[SKIP] Mapping module already completed, skipping execution.")

            _ensure_phageterm(input_fasta, read1, read2, sample_name,
                              mapping_dir, database, threads)
            return
    mapping_dir.mkdir(parents=True, exist_ok=True)
    logger, log_file = setup_logger(output_dir, "mapping")

    try:
        index_prefix = mapping_dir / "index"
        run_command(['bowtie2-build', '-f', input_fasta, str(index_prefix), '--threads', str(threads)],
                    "Step 1: Building Bowtie2 index...", logger)

        sam_file = mapping_dir / f"tmp.{sample_name}.sam"
        run_command(['bowtie2', '-1', read1, '-2', read2, '-p', str(threads), '-x', str(index_prefix), '-S', str(sam_file)],
                    "Step 2: Running Bowtie2 alignment...", logger)

        bam_file = mapping_dir / f"tmp.{sample_name}.bam"
        run_command(['samtools', 'view', '-@', str(threads), '-b', '-S', str(sam_file), '-o', str(bam_file)],
                    "Step 3: Converting SAM to BAM...", logger)

        sorted_bam = mapping_dir / f"{sample_name}.sorted.bam"
        run_command(['samtools', 'sort', '-@', str(threads), '-l', '9', '-O', 'BAM', str(bam_file), '-o', str(sorted_bam)],
                    "Step 4: Sorting BAM file...", logger)

        depth_file = mapping_dir / "base_depth.tsv"
        logger.info("Step 5: Calculating depth...")
        with open(depth_file, 'w') as f:
            subprocess.run(["samtools", "depth", str(sorted_bam)], stdout=f, check=True)

        logger.info("Step 6: Cleaning temporary files...")
        for pattern in ['tmp.*', 'index.*']:
            for f in mapping_dir.glob(pattern):
                f.unlink()

        plot_script = Path(database) / "mapping_profile" / "plot_depth.py"
        if plot_script.exists():
            run_command(['python', str(plot_script), str(depth_file), str(mapping_dir) + '/'],
                        "Step 7: Plotting depth profile...", logger)

        # ─── Step 8: PhageTerm ───
        run_phageterm(input_fasta, read1, read2, sample_name,
                      mapping_dir, database, threads, logger)

    except Exception as e:
        logger.error(f"Mapping pipeline error: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)


def _ensure_phageterm(input_fasta, read1, read2, sample_name,
                      mapping_dir, database, threads):

    marker_fasta = mapping_dir / f"{sample_name}_phageterm.fasta"
    if marker_fasta.exists():
        return  

    phageterm_dir = mapping_dir / "phageterm"
    report_pdf = phageterm_dir / f"{sample_name}_PhageTerm_report.pdf"
    stats_csv = phageterm_dir / f"{sample_name}_statistics.csv"
    if report_pdf.exists() or stats_csv.exists():
        return  

    output_dir = mapping_dir.parent
    logger, _ = setup_logger(output_dir, "mapping_phageterm")
    logger.info("[PhageTerm] Mapping was skipped but PhageTerm has not been run yet. Running now...")
    run_phageterm(input_fasta, read1, read2, sample_name,
                  mapping_dir, database, threads, logger)


def main():
    parser = argparse.ArgumentParser(description='JMPhage Mapping Module')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-1', '--read1', required=True, help='Forward reads (R1)')
    parser.add_argument('-2', '--read2', required=True, help='Reverse reads (R2)')
    parser.add_argument('-d', '--database', required=True, help='Database path')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads (default: 4)')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    args = parser.parse_args()

    run_mapping(args.input, args.read1, args.read2, args.threads, args.database, args.output)

if __name__ == "__main__":
    main()