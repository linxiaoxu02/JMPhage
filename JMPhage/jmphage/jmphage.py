#!/usr/bin/env python3

import sys
import os
import argparse
import hashlib
import shutil
import urllib.request
import urllib.error
from pathlib import Path
from datetime import datetime
import subprocess
import json

def get_module_dir():
    current_file = Path(__file__).resolve()
    module_path = current_file.parent / "modules"
    if module_path.exists() and module_path.is_dir():
        return module_path
    else:
        print(f"FATAL ERROR: Could not find the modules directory!")
        print(f"Expected path: {module_path}")
        print(f"Please check if the installation is complete or if the 'modules' folder was deleted.")
        sys.exit(1)

MODULE_DIR = get_module_dir()
VERSION = "1.0"

FASTA_SUFFIX = ".fasta"
R1_SUFFIX = ".R1.fastq"
R2_SUFFIX = ".R2.fastq"

# ============= jm_db configuration =============
ZENODO_RECORD_ID = "20838733"
DB_FILENAME      = "jm_db.tar.gz"
DB_MD5           = "20f8db131a251b150c5598ac192f6c15"
DB_SIZE_GB       = 4.9
DB_UNPACK_GB     = 16
ZENODO_DOI       = f"10.5281/zenodo.{ZENODO_RECORD_ID}"
ZENODO_URL       = f"https://zenodo.org/records/{ZENODO_RECORD_ID}/files/{DB_FILENAME}"
# ===============================================


def print_banner():
    banner = f"""

                                                                     ████╗
                                                                   ████████╗
   ██╗ ███╗   ███╗ ██████╗  ██╗  ██╗  █████╗   ██████╗  ███████╗   ████████║
   ██║ ████╗ ████║ ██╔══██╗ ██║  ██║ ██╔══██╗ ██╔════╝  ██╔════╝    ██████║
   ██║ ██╔████╔██║ ██████╔╝ ███████║ ███████║ ██║  ███╗ █████╗        ██║ 
██ ██║ ██║╚██╔╝██║ ██╔═══╝  ██╔══██║ ██╔══██║ ██║   ██║ ██╔══╝        ██║  
█████║ ██║ ╚═╝ ██║ ██║      ██║  ██║ ██║  ██║ ╚██████╔╝ ███████╗    ██╗ ██╗
╚════╝ ╚═╝     ╚═╝ ╚═╝      ╚═╝  ╚═╝ ╚═╝  ╚═╝  ╚═════╝  ╚══════╝  ██╗═╝   ██╗
                                                                ██╗═╝       ██╗
                                                                ╚═╝         ╚═╝
      A modular pipeline for tailed phage (Caudoviricetes) analysis
                              v{VERSION}
"""
    print(banner)



def scan_input(input_path, reads_dir=None, require_reads=False):

    input_path = Path(input_path).resolve()
    phages = []
    
    if input_path.is_file():
        
        if not input_path.name.endswith(FASTA_SUFFIX):
            print(f"[ERROR] Single-file input must end with '{FASTA_SUFFIX}': {input_path}")
            sys.exit(1)
        fasta_files = [input_path]
    elif input_path.is_dir():
        
        fasta_files = sorted(input_path.glob(f"*{FASTA_SUFFIX}"))
        if not fasta_files:
            print(f"[ERROR] No '*{FASTA_SUFFIX}' files found in directory: {input_path}")
            sys.exit(1)
    else:
        print(f"[ERROR] Input path does not exist: {input_path}")
        sys.exit(1)

    reads_search_dirs = []
    if reads_dir:
        reads_search_dirs.append(Path(reads_dir).resolve())
    if input_path.is_dir():
        reads_search_dirs.append(input_path)
    else:
        reads_search_dirs.append(input_path.parent)

    for fa in fasta_files:
        name = fa.name[:-len(FASTA_SUFFIX)]
        r1, r2 = None, None
        for d in reads_search_dirs:
            cand_r1 = d / f"{name}{R1_SUFFIX}"
            cand_r2 = d / f"{name}{R2_SUFFIX}"
            if cand_r1.exists() and cand_r2.exists():
                r1, r2 = cand_r1, cand_r2
                break

        has_reads = (r1 is not None and r2 is not None)
        phages.append({
            'name': name,
            'fasta': fa,
            'r1': r1,
            'r2': r2,
            'has_reads': has_reads
        })

    if require_reads:
        missing = [p['name'] for p in phages if not p['has_reads']]
        if missing:
            print(f"[ERROR] The 'mapping' module requires reads, but reads were not found for:")
            for n in missing:
                print(f"        - {n} (expected: {n}{R1_SUFFIX} / {n}{R2_SUFFIX})")
            print(f"  Searched in: {', '.join(str(d) for d in reads_search_dirs)}")
            sys.exit(1)

    return phages


def print_phage_manifest(phages, mode="all"):
    
    print(f"\n{'='*60}")
    print(f"  Detected {len(phages)} phage(s) (mode: {mode})")
    print(f"{'='*60}")
    for p in phages:
        if mode in ("mapping",):
            reads_tag = "[reads OK]"
        elif mode in ("annotation", "characterization", "no_mapping"):
            reads_tag = "[no reads needed]"
        else:  # all
            reads_tag = "[reads OK]" if p['has_reads'] else "[NO READS -> mapping skipped]"
        print(f"  - {p['name']}  {reads_tag}")
    print(f"{'='*60}\n")



def run_subprocess_module(env_name, script_name, extra_args):
    script_path = str(MODULE_DIR / script_name)
    if not os.path.exists(script_path):
        print(f" Error: Submodule script not found: {script_path}")
        sys.exit(1)

    cmd = [
        "conda", "run", "--no-capture-output", "-n", env_name,
        "python", script_path
    ]
    cmd.extend(extra_args)

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\n Execution Failed: Module [{script_name}] exited in environment [{env_name}] with error code: {e.returncode}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        print(f" Error: 'conda' command not found. Please ensure Conda is installed and initialized.")
        sys.exit(1)


# ===================== Install: conda environments =====================

def install_envs(env_dir):
    """Install all required conda environments from .yml files."""
    env_list = ["JM_mapping", "JM_annotation", "JM_characterization", "jmp_VIRIDIC", "jmp_GLUVAB", "jmp_Phageterm"]
    env_dir_path = Path(env_dir).resolve()

    print("\n=== Starting JMPhage Environment Deployment ===")
    print(f"Reading YML config directory: {env_dir_path}")

    try:
        subprocess.run(["conda", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[ERROR] Conda not detected. Please ensure Miniconda/Anaconda is installed and added to PATH.")
        sys.exit(1)

    try:
        res = subprocess.run(["conda", "info", "--envs", "--json"], capture_output=True, text=True)
        envs_info = json.loads(res.stdout)
        existing_envs = [Path(p).name for p in envs_info['envs']]
    except Exception as e:
        print(f"[WARNING] Failed to retrieve conda environments. Will attempt forced installation. Error: {e}")
        existing_envs = []

    for env_name in env_list:
        yml_file = env_dir_path / f"{env_name}.yml"
        print("-" * 60)
        if env_name in existing_envs:
            print(f"[INFO] Environment [{env_name}] already exists. Skipping installation.")
            print(f"       (To reinstall, manually run: conda env remove -n {env_name} first)")
            continue
        if not yml_file.exists():
            print(f"[ERROR] Configuration file not found: {yml_file}")
            sys.exit(1)
        print(f"-> Creating environment: {env_name} ...")
        try:
            subprocess.run(["conda", "env", "create", "-f", str(yml_file), "-n", env_name], check=True)
            print(f"[OK] Environment [{env_name}] created successfully!")
        except subprocess.CalledProcessError:
            print(f"[ERROR] Failed to create environment [{env_name}]. Please check network connection or {yml_file} configuration.")
            sys.exit(1)

    print("=" * 60)
    print("[OK] All conda environments deployed successfully!")


# ===================== Install: jm_db database =====================

def compute_md5(filepath, chunk_size=8 * 1024 * 1024):
    """Compute MD5 of a file with progress display."""
    md5 = hashlib.md5()
    total = os.path.getsize(filepath)
    done = 0
    with open(filepath, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            md5.update(chunk)
            done += len(chunk)
            pct = done * 100 / total if total else 0
            print(f"\r  Verifying MD5 ... {pct:5.1f}%", end="", flush=True)
    print()
    return md5.hexdigest()


def download_with_wget(url, dest_path):
    """Download using wget (resumable, with progress bar). Falls back to urllib if wget missing."""
    dest_path = Path(dest_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        subprocess.run(["wget", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        has_wget = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        has_wget = False

    if has_wget:
        print(f"      Using wget (resumable download)")
        cmd = ["wget", "-c", "-O", str(dest_path), url]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"\n[ERROR] wget failed (exit code {e.returncode})")
            print(f"        Partial file may exist at {dest_path} -- re-run install to resume.")
            sys.exit(1)
        except KeyboardInterrupt:
            print(f"\n[INFO] Interrupted. Re-run install to resume from {dest_path}.")
            sys.exit(1)
        return

    # ---- Fallback: urllib with Range header ----
    print(f"      wget not found, using Python urllib (resumable download)")
    tmp_path = dest_path.with_suffix(dest_path.suffix + ".part")
    existing_size = tmp_path.stat().st_size if tmp_path.exists() else 0
    headers = {}
    if existing_size > 0:
        headers["Range"] = f"bytes={existing_size}-"
        print(f"      Resuming from {existing_size / 1e9:.2f} GB ...")

    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            total_size = int(resp.headers.get("Content-Length", 0)) + existing_size
            mode = "ab" if existing_size > 0 else "wb"
            done = existing_size
            with open(tmp_path, mode) as f:
                while True:
                    chunk = resp.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
                    done += len(chunk)
                    if total_size:
                        pct = done * 100 / total_size
                        print(f"\r      Downloading ... {done/1e9:5.2f}/{total_size/1e9:.2f} GB ({pct:5.1f}%)",
                              end="", flush=True)
            print()
    except urllib.error.URLError as e:
        print(f"\n[ERROR] Download failed: {e}")
        print(f"        Partial file kept at {tmp_path} -- re-run install to resume.")
        sys.exit(1)
    except KeyboardInterrupt:
        print(f"\n[INFO] Interrupted. Partial file at {tmp_path} -- re-run install to resume.")
        sys.exit(1)

    tmp_path.rename(dest_path)


def install_database(db_dir):
    """Install jm_db: download from Zenodo, verify MD5, extract."""
    db_dir = Path(db_dir).resolve()
    db_dir.mkdir(parents=True, exist_ok=True)
    marker = db_dir / ".jm_db_installed"

    if marker.exists():
        print(f"\n[INFO] Database already installed at: {db_dir}")
        print(f"       To force reinstall, remove the marker file: {marker}")
        return

    print("\n=== Starting jm_db Database Installation ===")
    print(f"Zenodo DOI : {ZENODO_DOI}")
    print(f"Target dir : {db_dir}")
    print(f"Archive    : ~{DB_SIZE_GB} GB | Unpacked: ~{DB_UNPACK_GB} GB")
    free_gb = shutil.disk_usage(db_dir).free / 1e9
    print(f"Free space : {free_gb:.1f} GB")
    if free_gb < DB_UNPACK_GB:
        print(f"[WARNING] Free space may be insufficient (need ~{DB_UNPACK_GB} GB).")

    archive_path = db_dir / DB_FILENAME

    # ---- Step 1: Download ----
    print(f"\n[1/3] Downloading from Zenodo")
    print(f"      URL : {ZENODO_URL}")
    print(f"      Dest: {archive_path}")
    print(f"      (Offline servers: download {DB_FILENAME} manually from")
    print(f"       https://zenodo.org/records/{ZENODO_RECORD_ID} and place it")
    print(f"       into {db_dir} before running this command.)")
    if archive_path.exists():
        print(f"      Found existing archive, skipping download.")
    else:
        download_with_wget(ZENODO_URL, archive_path)

    # ---- Step 2: MD5 verification ----
    print(f"\n[2/3] Verifying MD5 ...")
    print(f"      Expected: {DB_MD5}")
    actual = compute_md5(archive_path)
    print(f"      Actual  : {actual}")
    if actual != DB_MD5:
        print(f"[ERROR] MD5 mismatch! The file may be corrupted.")
        print(f"        Delete {archive_path} and re-run install.")
        sys.exit(1)
    print(f"[OK] MD5 verified.")

    # ---- Step 3: Extract ----
    print(f"\n[3/3] Extracting to {db_dir} ...")
    try:
        subprocess.run(["tar", "-xzf", str(archive_path), "-C", str(db_dir)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Extraction failed: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"[ERROR] 'tar' command not found.")
        sys.exit(1)
    print(f"[OK] Extracted.")

    # ---- Cleanup & marker ----
    print(f"\n[INFO] Removing archive to save space: {archive_path}")
    archive_path.unlink()

    marker.write_text(
        f"jm_db v1.0\n"
        f"DOI: {ZENODO_DOI}\n"
        f"MD5: {DB_MD5}\n"
        f"Installed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    )

    print("\n" + "=" * 60)
    print(f"[OK] jm_db installed at: {db_dir}")
    print(f"     Use this path with -d/--database in subsequent commands.")
    print(f"     Example: jmphage all -i my_phage.fasta -d {db_dir} -o results")
    print("=" * 60)


# ===================== Pipeline modules =====================

def phage_output_dir(base_out, phage_name, n_phages):

    if n_phages == 1:
        return Path(base_out)
    return Path(base_out) / phage_name


def run_mapping_for_phages(phages, database, threads, base_out):

    print(f"\n[Step] Mapping ({sum(1 for p in phages if p['has_reads'])} phage(s) with reads)")
    n = len(phages)
    for p in phages:
        if not p['has_reads']:
            print(f"  [SKIP] {p['name']}: no reads available")
            continue
        phage_out = phage_output_dir(base_out, p['name'], n)
        phage_out.mkdir(parents=True, exist_ok=True)
        print(f"\n  >>> Mapping: {p['name']}")
        run_subprocess_module(
            "JM_mapping", "run_mapping.py",
            ["-i", str(p['fasta']),
             "-1", str(p['r1']), "-2", str(p['r2']),
             "-d", database, "-t", str(threads),
             "-o", str(phage_out)]
        )


def update_fasta_from_phageterm(phages, base_out):

    n = len(phages)
    updated = 0
    for p in phages:
        phage_out = phage_output_dir(base_out, p['name'], n)
        pt_fasta = phage_out / "1.read_mapping" / f"{p['name']}_phageterm.fasta"
        if pt_fasta.exists() and pt_fasta.stat().st_size > 100:
            renamed = phage_out / f"{p['name']}.fasta"
            if not renamed.exists():
                shutil.copy(pt_fasta, renamed)
            print(f"  [PhageTerm] {p['name']}: using reoriented sequence -> {renamed.name}")
            p['fasta'] = renamed
            updated += 1
    if updated:
        print(f"  [PhageTerm] {updated} phage(s) switched to PhageTerm-reoriented sequences.")
        skipped = len(phages) - updated
        if skipped > 0:
            print(f"  [PhageTerm] {skipped} phage(s) had no valid reoriented sequence; using original FASTA.")
    else:
        print(f"  [PhageTerm] --use-phageterm is enabled, but no valid reoriented sequences were found.")
        print(f"  [PhageTerm] Possible reasons: no fixed terminus detected (e.g. headful packaging).")
        print(f"  [PhageTerm] All phages will use their original FASTA for downstream analysis.")

def run_annotation_for_phages(phages, database, threads, base_out):

    print(f"\n[Step] Annotation ({len(phages)} phage(s))")
    n = len(phages)
    for p in phages:
        phage_out = phage_output_dir(base_out, p['name'], n)
        phage_out.mkdir(parents=True, exist_ok=True)
        print(f"\n  >>> Annotation: {p['name']}")
        run_subprocess_module(
            "JM_annotation", "run_annotation.py",
            ["-i", str(p['fasta']),
             "-d", database, "-t", str(threads),
             "-o", str(phage_out)]
        )


def run_characterization_for_phages(phages, database, threads, base_out):

    print(f"\n[Step] Characterization ({len(phages)} phage(s))")

    if len(phages) == 1:
        p = phages[0]
        phage_out = phage_output_dir(base_out, p['name'], 1)  # = base_out
        phage_out.mkdir(parents=True, exist_ok=True)
        print(f"  (Single-phage mode: legacy logic, pulling neighbors from DB)")
        run_subprocess_module(
            "JM_characterization", "run_characterization.py",
            ["-i", str(p['fasta']),
             "-d", database, "-t", str(threads),
             "-o", str(phage_out)]
        )
        return


    manifest = {
        'base_output': str(base_out),
        'database': database,
        'threads': threads,
        'phages': [
            {'name': p['name'], 'fasta': str(p['fasta'])} for p in phages
        ]
    }
    manifest_path = base_out / "_joint_manifest.json"
    base_out.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"  (Joint mode: shared vConTACT2 once, then per-phage collinearity + joint ANI/tree)")
    run_subprocess_module(
        "JM_characterization", "run_characterization.py",
        ["--joint", str(manifest_path)]
    )


def write_pipeline_report(phages, base_out, mode, start_time, end_time):

    if len(phages) < 2:
        return
    summary_dir = base_out / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    report = summary_dir / "pipeline_report.tsv"
    with open(report, 'w') as f:
        f.write("phage_name\tfasta\tr1\tr2\thas_reads\tmapping_status\tnotes\n")
        for p in phages:
            if mode in ("all",):
                mapping_status = "done" if p['has_reads'] else "skipped (no reads)"
            elif mode == "mapping":
                mapping_status = "done"
            else:
                mapping_status = "not run (mode=" + mode + ")"
            note = "" if p['has_reads'] else "no paired reads found"
            f.write(f"{p['name']}\t{p['fasta']}\t{p['r1'] or ''}\t{p['r2'] or ''}\t"
                    f"{p['has_reads']}\t{mapping_status}\t{note}\n")
    print(f"\n[Report] Pipeline summary -> {report}")




def create_main_parser():
    naming_rules_epilog = """
======================================================================
CRITICAL NAMING RULES (Strict Mode):
======================================================================
To ensure the pipeline runs smoothly, JMPhage enforces strict naming 
conventions for all input files and internal FASTA headers.

1. File Extensions (Case-Sensitive):
   - Genome Assembly : Must end with '{fasta}' (e.g., sample1{fasta})
   - Forward Reads   : Must end with '{r1}' (e.g., sample1{r1})
   - Reverse Reads   : Must end with '{r2}' (e.g., sample1{r2})

2. Name Consistency:
   - The prefix of the reads MUST exactly match the FASTA file name.
   - The primary ID inside the FASTA header (the text right after '>') 
     MUST be identical to the file name.

Correct Example:
   File Name:   myphage01.fasta
   FASTA Header:>myphage01
   Reads Pairs: myphage01.R1.fastq / myphage01.R2.fastq

Incorrect Example (Pipeline Will Fail):
   File Name:   myphage_v1.fasta
   FASTA Header:>Sequence_1          <-- [ERROR: Mismatch with file name]
   Reads Pairs: myphage.R1.fastq     <-- [ERROR: Mismatch with file name]
======================================================================
    """.format(fasta=FASTA_SUFFIX, r1=R1_SUFFIX, r2=R2_SUFFIX)

    parser = argparse.ArgumentParser(
        description='JMPhage - Modular Phage Analysis Pipeline (supports single phage or batch directory input)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        epilog=naming_rules_epilog
    )
    parser.add_argument('-v', '--version', action='version', version=f'JMPhage v{VERSION}')
    parser.add_argument('-h', '--help', action='store_true', help='Show help message and exit')
    subparsers = parser.add_subparsers(dest='module', help='Select a module to run')

    # install (envs and/or database)
    parser_install = subparsers.add_parser(
        'install',
        help='Install Conda environments and/or jm_db reference database',
        description='Install Conda environments (-e) and/or jm_db reference database (-d). '
                    'At least one of -e and -d must be provided.'
    )
    parser_install.add_argument('-e', '--env_dir', default=None,
                                help='Directory containing JMPhage Conda environment YAML files (jm_envs)')
    parser_install.add_argument('-d', '--db_dir', default=None,
                                help='Directory to install jm_db')

    common_args = argparse.ArgumentParser(add_help=False)
    common_args.add_argument('-i', '--input', required=True,
                             help=f'Input FASTA file OR directory containing *{FASTA_SUFFIX} files (max 5 phages per run)')
    common_args.add_argument('-d', '--database', required=True, help='Path to JMP database')
    common_args.add_argument('-t', '--threads', type=int, default=4, help='Number of threads (default: 4)')
    common_args.add_argument('-o', '--output', required=True, help='Output directory')

    # mapping
    parser_mapping = subparsers.add_parser('mapping', parents=[common_args], help='Read mapping and depth calculation')
    parser_mapping.add_argument('-1', '--read1', default=None,
                                help=f'Forward reads (R1). Only used when -i is a single FASTA file. '
                                     f'For directory input, reads are auto-paired by name (*{R1_SUFFIX}).')
    parser_mapping.add_argument('-2', '--read2', default=None,
                                help=f'Reverse reads (R2). Only used when -i is a single FASTA file.')
    parser_mapping.add_argument('-r', '--reads-dir', default=None,
                                help='Optional directory containing reads (default: same dir as FASTA)')

    # annotation / characterization / no_mapping
    subparsers.add_parser('annotation', parents=[common_args], help='Comprehensive genome annotation')
    subparsers.add_parser('characterization', parents=[common_args], help='Genomic characterization')
    subparsers.add_parser('no_mapping', parents=[common_args], help='Annotation + characterization (no mapping)')

    # all
    parser_all = subparsers.add_parser('all', parents=[common_args], help='Run the complete analysis pipeline')
    parser_all.add_argument('-1', '--read1', default=None,
                            help=f'Forward reads (R1). Only used when -i is a single FASTA file.')
    parser_all.add_argument('-2', '--read2', default=None,
                            help=f'Reverse reads (R2). Only used when -i is a single FASTA file.')
    parser_all.add_argument('-r', '--reads-dir', default=None,
                            help='Optional directory containing reads (default: same dir as FASTA / input dir)')
    parser_all.add_argument('--use-phageterm', action='store_true', default=False,
                            help='Use PhageTerm reoriented sequences for downstream analysis. '
                                 'Default: only generate PhageTerm report without replacing FASTA.')
    return parser


def build_phage_list(args, require_reads):

    input_path = Path(args.input).resolve()
    reads_dir = getattr(args, 'reads_dir', None)

    if input_path.is_file():
        explicit_r1 = getattr(args, 'read1', None)
        explicit_r2 = getattr(args, 'read2', None)
        if explicit_r1 and explicit_r2:
            name = input_path.name
            if name.endswith(FASTA_SUFFIX):
                name = name[:-len(FASTA_SUFFIX)]
            else:
                while '.' in name:
                    name = name.rsplit('.', 1)[0]
            r1, r2 = Path(explicit_r1).resolve(), Path(explicit_r2).resolve()
            if require_reads:
                for r in (r1, r2):
                    if not r.exists():
                        print(f"[ERROR] Reads file not found: {r}")
                        sys.exit(1)
            return [{'name': name, 'fasta': input_path,
                     'r1': r1 if r1.exists() else None,
                     'r2': r2 if r2.exists() else None,
                     'has_reads': r1.exists() and r2.exists()}]

    return scan_input(args.input, reads_dir=reads_dir, require_reads=require_reads)



def main():
    parser = create_main_parser()

    if len(sys.argv) == 1:
        print_banner()
        parser.print_help()
        sys.exit(0)

    args, unknown = parser.parse_known_args()

    if args.help or not args.module:
        print_banner()
        parser.print_help()
        sys.exit(0)

    print_banner()

    if args.module == 'install':
        if not args.env_dir and not args.db_dir:
            print("[ERROR] Specify at least one of: -e (envs) and/or -d (database)")
            print("        Examples:")
            print("          jmphage install -e ./yml                 # envs only")
            print("          jmphage install -d /data/jm_db           # database only")
            print("          jmphage install -e ./yml -d /data/jm_db  # both")
            sys.exit(1)
        if args.env_dir:
            install_envs(args.env_dir)
        if args.db_dir:
            install_database(args.db_dir)
        print("\n[OK] Installation finished. You can now run 'jmphage <module> ...'")
        return

    base_out = Path(args.output).resolve()
    base_out.mkdir(parents=True, exist_ok=True)

    if args.module == 'mapping':
        phages = build_phage_list(args, require_reads=True)
        print_phage_manifest(phages, mode='mapping')
        start_time = datetime.now()
        run_mapping_for_phages(phages, args.database, args.threads, base_out)
        end_time = datetime.now()
        write_pipeline_report(phages, base_out, 'mapping', start_time, end_time)

    elif args.module == 'annotation':
        phages = build_phage_list(args, require_reads=False)
        print_phage_manifest(phages, mode='annotation')
        run_annotation_for_phages(phages, args.database, args.threads, base_out)

    elif args.module == 'characterization':
        phages = build_phage_list(args, require_reads=False)
        print_phage_manifest(phages, mode='characterization')
        run_characterization_for_phages(phages, args.database, args.threads, base_out)

    elif args.module == 'no_mapping':
        phages = build_phage_list(args, require_reads=False)
        print_phage_manifest(phages, mode='no_mapping')
        start_time = datetime.now()
        print(f"\nStarting pipeline (no_mapping): {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

        print("\n[Stage 1/2] Annotation")
        run_annotation_for_phages(phages, args.database, args.threads, base_out)

        print("\n[Stage 2/2] Characterization")
        run_characterization_for_phages(phages, args.database, args.threads, base_out)

        end_time = datetime.now()
        write_pipeline_report(phages, base_out, 'no_mapping', start_time, end_time)
        duration = (end_time - start_time).total_seconds()
        print("\n" + "═" * 60)
        print(f"JMPhage Pipeline (no_mapping) Completed!")
        print(f"Total Time: {duration/60:.2f} minutes")
        print("═" * 60)

    elif args.module == 'all':
        phages = build_phage_list(args, require_reads=False)
        print_phage_manifest(phages, mode='all')
        start_time = datetime.now()
        print(f"\nStarting full pipeline: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

        print("\n[Stage 1/3] Mapping (per-phage, skipped if no reads)")
        run_mapping_for_phages(phages, args.database, args.threads, base_out)

        if args.use_phageterm:
            update_fasta_from_phageterm(phages, base_out)
        else:
            print("  [PhageTerm] Reports generated. Use --use-phageterm to replace FASTA for downstream analysis.")

        print("\n[Stage 2/3] Annotation (per-phage)")
        run_annotation_for_phages(phages, args.database, args.threads, base_out)

        print("\n[Stage 3/3] Characterization")
        run_characterization_for_phages(phages, args.database, args.threads, base_out)

        end_time = datetime.now()
        write_pipeline_report(phages, base_out, 'all', start_time, end_time)
        duration = (end_time - start_time).total_seconds()
        print("\n" + "═" * 60)
        print(f"JMPhage Pipeline Analysis Completed!")
        print(f"Phages analyzed: {len(phages)}")
        print(f"Total Time: {duration/60:.2f} minutes")
        print(f"Output: {base_out}")
        print("═" * 60)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n User interrupted the program. Exiting...")
        sys.exit(1)
