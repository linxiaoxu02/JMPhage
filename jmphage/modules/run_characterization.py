#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
from pathlib import Path
import json
import logging
from datetime import datetime
import traceback
import shutil
import re

ENV_NAME = "JM_characterization"


def setup_logger(output_dir, module_name="characterization"):
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{module_name}_{timestamp}.log"

    logger = logging.getLogger(module_name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(
        '[%(asctime)s] [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_formatter)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger, log_file


def check_file_exists(filepath, description="File", logger=None):
    if not os.path.exists(filepath):
        error_msg = f"[ERROR] {description} not found: {filepath}"
        if logger:
            logger.error(error_msg)
        else:
            print(error_msg)
        sys.exit(1)
    if logger:
        logger.info(f"[CHECK] {description} passed: {filepath}")
    return True


def get_sample_name(input_fasta):
    filename = Path(input_fasta).name
    sample_name = filename
    while '.' in sample_name:
        sample_name = sample_name.rsplit('.', 1)[0]
    return sample_name


def run_command(cmd, description="", logger=None, env_name=ENV_NAME, cwd=None, quiet=False, **kwargs):
    if logger:
        logger.debug(f"\n{'='*20} PROCESS START: {description} {'='*20}")
        if not quiet:
            logger.info(f"\n{'='*60}\n{description}\n{'='*60}")

    env_vars = os.environ.copy()
    for var in ["PERL5LIB", "PYTHONPATH", "PYTHONHOME", "PERL_LOCAL_LIB_ROOT", "LD_LIBRARY_PATH"]:
        env_vars.pop(var, None)

    try:
        search_cmd = ["conda", "info", "--envs", "--json"]
        res = subprocess.run(search_cmd, capture_output=True, text=True)
        if res.returncode == 0:
            envs_info = json.loads(res.stdout)
            target_prefix = next((p for p in envs_info['envs'] if Path(p).name == env_name), None)
            if target_prefix:
                target_bin = os.path.join(target_prefix, "bin")
                env_vars["PATH"] = target_bin + os.pathsep + env_vars.get("PATH", "")
    except Exception as e:
        if logger:
            logger.debug(f"Path priority adjustment failed: {e}")

    if isinstance(cmd, str):
        cmd = cmd.split()

    should_check = kwargs.get('check', True)
    full_cmd = ["conda", "run", "-n", env_name] + cmd if cmd[0] != 'conda' else cmd

    if logger:
        logger.debug(f"Executing command: {' '.join(full_cmd)}")

    try:
        process = subprocess.Popen(
            full_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env_vars,
            cwd=cwd,
            bufsize=1
        )

        for line in process.stdout:
            line_clean = line.rstrip()
            if line_clean and logger:
                logger.debug(f"[{description}] {line_clean}")

        process.wait()

        if process.returncode != 0:
            if logger:
                logger.error(f"[FAILED] {description} (Exit Code: {process.returncode})")
            if should_check:
                raise subprocess.CalledProcessError(process.returncode, full_cmd)

        if logger:
            logger.debug(f"{'='*20} PROCESS END: {description} (Code: {process.returncode}) {'='*20}")

        return process.returncode

    except Exception as e:
        if logger:
            logger.error(f"Unexpected error during '{description}': {str(e)}")
            logger.debug(traceback.format_exc())
        if should_check:
            raise
        return 1


def fix_vcontact2_compatibility(logger=None):
    if logger: logger.info(f"[Auto-Fix] Checking vConTACT2 compatibility in {ENV_NAME}...")

    vcontact2_path = None
    try:
        cmd = ["conda", "run", "-n", ENV_NAME, "python", "-c", "import vcontact2; print(vcontact2.__path__[0])"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            vcontact2_path = Path(result.stdout.strip())
            if logger: logger.info(f"[OK] Source located: {vcontact2_path}")
    except Exception as e:
        if logger: logger.debug(f"Detection failed: {e}")

    if not vcontact2_path:
        return

    fixes = [
        (vcontact2_path / "pcprofiles.py", "zip(*profiles.values)", "list(zip(*profiles.values))"),
        (vcontact2_path / "matrices.py", "zip(*xy)", "list(zip(*xy))"),
        (vcontact2_path / "modules.py", "N[:, m]", "N[:, int(m)]"),
        (vcontact2_path / "modules.py", "zip(*xy)", "list(zip(*xy))"),
        (vcontact2_path / "modules.py", 'suffixes=["_module", "_cluster"]', 'suffixes=["_mod", "_cl"]'),
        (vcontact2_path / "exports" / "summaries.py", "np.warnings.filterwarnings", "import warnings; warnings.filterwarnings")
    ]

    for f_path, old_str, new_str in fixes:
        if not f_path.exists(): continue
        try:
            content = f_path.read_text()
            if old_str in content and new_str not in content:
                f_path.write_text(content.replace(old_str, new_str))
                if logger: logger.info(f"[FIXED] {f_path.name}: {old_str} -> {new_str}")
        except Exception as e:
            if logger: logger.error(f"[ERROR] Failed to fix {f_path.name}: {e}")


def write_gene2genome_csv(faa_file, contig_id, output_csv):
    """Generate vConTACT2-compatible gene-to-genome mapping CSV.

    Replaces `vcontact2_gene2genome --source-type Prodigal-FAA`, which derives
    the contig_id by stripping the trailing `_N` from each protein_id. That
    works for prokka output (`P911113_00001` -> `P911113`) but breaks for
    pharokka/phanotate output (`P911113_CDS_0001` -> `P911113_CDS`), so we
    write the mapping ourselves with the correct contig_id.
    """
    with open(faa_file) as fin, open(output_csv, 'w') as fout:
        fout.write("protein_id,contig_id,keywords\n")
        for line in fin:
            if line.startswith('>'):
                pid = line.strip()[1:].split()[0]
                fout.write(f"{pid},{contig_id},\n")


########single input########################

def build_shared_network(input_fasta, sample_name, output_dir, database, threads, logger):
    logger.info(f"\n{'#'*60}\n# Step 1: Building Shared Network (Custom DB Mode)\n{'#'*60}\n")

    ntw_dir = output_dir / "5.shared_network"
    vcontact2_output = ntw_dir / "vcontact2_output"
    c1_ntw = vcontact2_output / "c1.ntw"

    ntw_dir.mkdir(parents=True, exist_ok=True)

    if c1_ntw.exists():
        logger.info("[INFO] Existing vConTACT2 result found. Skipping...")
        return ntw_dir

    pharokka_output = ntw_dir / "pharokka_output"
    input_path = Path(input_fasta).resolve()
    pharokka_db = Path(database) / "pharokka_db"

    if not pharokka_db.exists():
        logger.error(f"[ERROR] Pharokka database not found: {pharokka_db}")
        return None

    cmd = [
        'pharokka.py',
        '-i', str(input_path),
        '-o', str(pharokka_output),
        '-p', sample_name,
        '-l', sample_name,
        '-t', str(threads),
        '-f',
        '-d', str(pharokka_db)
    ]
    run_command(cmd, "Step 1.1: Running Pharokka...", logger)

    faa_file = pharokka_output / "phanotate.faa"
    if not faa_file.exists():
        logger.error("[ERROR] Pharokka failed to generate phanotate.faa")
        return None

    logger.info("Step 1.2: Cleaning FASTA headers...")
    cleaned_lines = []
    with open(faa_file, 'r') as f:
        for line in f:
            cleaned_lines.append(line.split(' ')[0] + '\n' if line.startswith('>') else line)
    with open(faa_file, 'w') as f:
        f.writelines(cleaned_lines)

    tmp_map_csv = ntw_dir / "out_map_tmp.csv"
    logger.info("Step 1.3: Mapping genes to genomes...")
    write_gene2genome_csv(faa_file, sample_name, tmp_map_csv)
    logger.info(f"[OK] Gene-to-genome mapping generated: {tmp_map_csv.name}")

    coline_db = Path(database) / "coline_profile"
    ref_map, ref_faa = coline_db / "out_map.csv", coline_db / "ref.faa"

    if not ref_map.exists() or not ref_faa.exists():
        logger.error("[ERROR] Missing database files (out_map.csv/ref.faa)")
        return None

    final_map_csv, all_faa = ntw_dir / "out_map.csv", ntw_dir / "all.faa"

    logger.info("Step 1.4: Merging with local database...")
    with open(final_map_csv, 'wb') as out_f:
        if tmp_map_csv.exists():
            with open(tmp_map_csv, 'rb') as in_f: shutil.copyfileobj(in_f, out_f)
        with open(ref_map, 'rb') as in_f: shutil.copyfileobj(in_f, out_f)
    if tmp_map_csv.exists(): tmp_map_csv.unlink()

    with open(all_faa, 'wb') as out_f:
        with open(faa_file, 'rb') as in_f: shutil.copyfileobj(in_f, out_f)
        with open(ref_faa, 'rb') as in_f: shutil.copyfileobj(in_f, out_f)

    fix_vcontact2_compatibility(logger)
    c1_jar = Path(database) / 'cluster_one-1.0.jar'

    cmd = [
        'vcontact2', '--raw-proteins', str(all_faa), '--rel-mode', 'Diamond',
        '--proteins-fp', str(final_map_csv), '--db', 'None',
        '--pcs-mode', 'MCL', '--vcs-mode', 'ClusterONE',
        '--c1-bin', str(c1_jar), '--output-dir', str(vcontact2_output),
        '-t', str(threads)
    ]
    run_command(cmd, "Step 1.5: Running vConTACT2 network analysis...", logger)

    return ntw_dir if c1_ntw.exists() else None


def run_collinearity_analysis(sample_name, ntw_dir, output_dir, database, threads, logger):
    logger.info(f"\n{'#'*60}\n# Step 2: Collinearity Analysis (Local DB Mode)\n{'#'*60}\n")

    base_dir = output_dir / "6.collinearity_analysis"
    output_html = base_dir / f"{sample_name}.html"

    if output_html.exists():
        logger.info("[INFO] Existing result found. Skipping...")
        return output_html

    pharokka_gbk = ntw_dir / "pharokka_output" / f"{sample_name}.gbk"
    pharokka_faa = ntw_dir / "pharokka_output" / "phanotate.faa"
    c1_ntw = ntw_dir / "vcontact2_output" / "c1.ntw"

    if not pharokka_gbk.exists() or not c1_ntw.exists():
        logger.error("[ERROR] Missing input files (GBK or c1.ntw)")
        return None

    gbk_dir, faa_dir = base_dir / "gbk", base_dir / "faa"
    gbk_dir.mkdir(parents=True, exist_ok=True)
    faa_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Step 2.1: Extracting Top 3 related genomes...")
    try:
        sample_lines = []
        with open(c1_ntw, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3 and parts[0] == sample_name:
                    sample_lines.append(parts)
        sample_lines.sort(key=lambda x: float(x[2]), reverse=True)
        top_genomes = [line[1] for line in sample_lines[:3]]
    except Exception as e:
        logger.error(f"Failed to parse c1.ntw: {e}"); return None

    ref_gb_src = Path(database) / "coline_profile" / "ref_gb"
    for genome in top_genomes:
        src = ref_gb_src / f"{genome}.gb"
        if not src.exists(): src = ref_gb_src / f"{genome}.gbk"
        if src.exists(): shutil.copy(src, gbk_dir / f"{genome}.gb")

    gbk_faa_script = Path(database) / "coline_profile" / "gbk_faa.py"
    if gbk_faa_script.exists():
        run_command(['python', str(gbk_faa_script), str(gbk_dir), str(faa_dir)], "Extracting FAA from GBK...", logger, check=False)

    shutil.copy(pharokka_faa, faa_dir / f"{sample_name}.faa")
    shutil.copy(pharokka_gbk, gbk_dir / f"{sample_name}.gb")

    logger.info("Step 2.2: Annotating proteins with HMMscan...")
    hmm_db = Path(database) / "VOG_HMM" / "phage_all.hmm"
    phage_hmm_type = Path(database) / "VOG_HMM" / "phage_hmm_type.tsv"
    hmm_type_dict = {}
    if phage_hmm_type.exists():
        for line in phage_hmm_type.read_text().splitlines():
            parts = line.split('\t')
            if len(parts) >= 2: hmm_type_dict[parts[0]] = '\t'.join(parts[1:])

    for faa_file in faa_dir.glob("*.faa"):
        g_name = faa_file.stem
        sub_dir = faa_dir / g_name
        sub_dir.mkdir(exist_ok=True)
        domtbl = sub_dir / "results.domtblout"

        run_command(['hmmscan', '-E', '1e-5', '--cpu', str(threads), '--domtblout', str(domtbl), str(hmm_db), str(faa_file)], f"HMMscan: {g_name}", logger, check=False)

        if domtbl.exists():
            best_hits = {}
            for line in domtbl.read_text().splitlines():
                if line.startswith('#'): continue
                cols = line.split()
                if len(cols) >= 4 and cols[3] not in best_hits: best_hits[cols[3]] = cols[0]

            with open(sub_dir / "protein_product", 'w') as f_prod:
                for q in sorted(best_hits.keys()):
                    f_prod.write(f"{q}\t{best_hits[q]}\t{hmm_type_dict.get(best_hits[q], '')}\n")

            gbk_mod = Path(database) / "coline_profile" / "gbk_motified.py"
            if gbk_mod.exists():
                run_command(['python', str(gbk_mod), str(sub_dir / "protein_product"), str(gbk_dir / f"{g_name}.gb")], f"Updating GBK: {g_name}", logger, check=False)

    logger.info("Step 2.3: Finalizing Collinearity Visualization...")
    gene_func_csv = base_dir / "gene_functions.csv"
    merge_script = Path(database) / "coline_profile" / "merge_protein_products.py"
    prod_files = list(faa_dir.glob("*/protein_product"))
    if merge_script.exists():
        run_command(['python', str(merge_script), str(gene_func_csv)] + [str(p) for p in prod_files], "Merging products", logger, check=False)

    fix_script = Path(database) / "coline_profile" / "fix_gene_functions.py"
    if fix_script.exists():
        run_command(['python', str(fix_script), str(gene_func_csv), str(gbk_dir)], "Filling missing genes", logger, check=False)

    gbk_files = list(gbk_dir.glob("*.gb")) + list(gbk_dir.glob("*.gbk"))
    if gbk_files and gene_func_csv.exists():
        colour_map = Path(database) / "coline_profile" / "colour_map.csv"
        cmd = ['clinker', '--gene_functions', str(gene_func_csv), '--colour_map', str(colour_map), '-dso'] + [str(g) for g in gbk_files] + ['-p', str(output_html)]
        run_command(cmd, "Running Clinker...", logger, check=True)

    return output_html if output_html.exists() else None


def run_ani_analysis(sample_name, ntw_dir, output_dir, database, threads, logger, input_fasta):
    logger.info(f"\n{'#'*60}\n# Step 3: ANI Analysis (VIRIDIC)\n{'#'*60}\n")

    VIRIDIC_ENV = "jmp_VIRIDIC"
    input_fasta = Path(input_fasta).resolve()
    ani_dir = (output_dir / "7.ANI_analysis").resolve()
    ani_fasta_dir = ani_dir / "ANI_fasta"
    ani_res_dir = ani_dir / "ANI_result"
    merged_fasta = ani_dir / "ANI.fasta"
    viridic_out_dir = ani_res_dir / "04_VIRIDIC_out"
    heatmap_pdf = viridic_out_dir / "Heatmap.PDF"
    result_csv = viridic_out_dir / "sim_MA_genCol.csv"
    if heatmap_pdf.exists() or result_csv.exists():
        logger.info(f"[INFO] Skip: VIRIDIC core results found in {viridic_out_dir.name}. Skipping calculation.")
        return ani_res_dir

    c1_ntw = ntw_dir / "vcontact2_output" / "c1.ntw"
    if not c1_ntw.exists(): return None

    ani_fasta_dir.mkdir(parents=True, exist_ok=True)
    ani_res_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Step 3.1: Gathering Top 20 sequences...")
    try:
        with open(c1_ntw, 'r') as f:
            lines = [l.split() for l in f if l.startswith(sample_name)]
        lines.sort(key=lambda x: float(x[2]))
        top_genomes = [line[1] for line in lines[-20:]]
    except Exception: return None

    ref_src = Path(database) / "coline_profile" / "ref_fasta"
    count = 0
    for genome in top_genomes:
        src = ref_src / f"{genome}.fasta"
        if src.exists():
            shutil.copy(src, ani_fasta_dir / f"{genome}.fasta")
            count += 1

    if input_fasta.exists():
        shutil.copy(input_fasta, ani_fasta_dir / f"{sample_name}.fasta")
        count += 1

    if count < 2: return None

    with open(merged_fasta, 'w') as out_f:
        for f in ani_fasta_dir.glob("*.fasta"):
            out_f.write(f.read_text().strip() + "\n")

    logger.info("Step 3.2: Running VIRIDIC...")
    vir_src = Path(database) / "VIRIDIC"
    vir_tmp = ani_dir / "VIRIDIC_tool_tmp"
    if vir_tmp.exists(): shutil.rmtree(vir_tmp)
    shutil.copytree(vir_src, vir_tmp)

    script_dir = (vir_tmp / "stand_alone" / "viridic_scripts").resolve()
    cmd = ["Rscript", "00_viridic_master.R", f"projdir={ani_res_dir}", f"in={merged_fasta}", f"ncor={threads}"]
    run_command(cmd, "VIRIDIC R-calculation", logger, env_name=VIRIDIC_ENV, cwd=str(script_dir), check=True)

    shutil.rmtree(vir_tmp, ignore_errors=True)
    shutil.rmtree(ani_fasta_dir, ignore_errors=True)

    return ani_res_dir


def run_phylogenetic_tree(sample_name, ntw_dir, output_dir, database, threads, logger, input_fasta):
    logger.info(f"\n{'#'*60}\n# Step 4: Phylogenetic Tree (Small & Big)\n{'#'*60}\n")

    GLUVAB_ENV = "jmp_GLUVAB"
    phy_dir = (output_dir / "8.phylogenetic_tree").resolve()
    ref_dir = phy_dir / "ref"
    big_dir = phy_dir / "big_tree"
    small_dir = phy_dir / "small_tree"

    big_tree_out = big_dir / f"{sample_name}_Tree_With_node_IDs.newick"
    small_tree_out = small_dir / f"{sample_name}_Tree_With_node_IDs.newick"

    if small_tree_out.exists() and big_tree_out.exists():
        logger.info("[Tree] [SKIP] Both small and big trees already exist.")
        return phy_dir

    for d in [phy_dir, ref_dir, big_dir, small_dir]:
        d.mkdir(parents=True, exist_ok=True)

    glu_scr = (Path(database) / "GLUVAB_tree" / "GLUVAB_v0.6.pl").resolve()
    small_fasta = phy_dir / "small_tree_in.fasta"
    tree1 = phy_dir / "tree1.fasta"

    if small_tree_out.exists() and small_fasta.exists():

        logger.info("[Tree] [SKIP] Small tree exists; reusing small_tree_in.fasta for big tree.")
    else:

        if small_tree_out.exists() and not small_fasta.exists():
            logger.info("[Tree] Small tree exists but small_tree_in.fasta was cleaned up; "
                        "rebuilding it for big tree input.")

        logger.info("[Tree] Step 4.1: Gathering reference sequences...")
        with open(ntw_dir / "vcontact2_output" / "c1.ntw", 'r') as f:
            lines = [l.split() for l in f if l.startswith(sample_name)]
        lines.sort(key=lambda x: float(x[2]), reverse=True)
        lines = lines[:50]

        ref_src = Path(database) / "coline_profile" / "ref_fasta"
        for genome in [x[1] for x in lines]:
            src = ref_src / f"{genome}.fasta"
            if src.exists(): shutil.copy(src, ref_dir / f"{genome}.fasta")

        seen_ids = set()
        write_flag = True
        with open(tree1, 'w') as out_f:
            for f_path in [Path(input_fasta)] + list(ref_dir.glob("*.fasta")):
                for line in f_path.read_text().splitlines():
                    if line.startswith('>'):
                        clean_id = re.sub(r'\..*', '', line).strip()
                        if clean_id in seen_ids:
                            write_flag = False
                            logger.debug(f"[Tree] Skipping duplicate in tree1: {clean_id}")
                            continue
                        seen_ids.add(clean_id)
                        write_flag = True
                        out_f.write(clean_id + "\n")
                    elif write_flag:
                        out_f.write(line + "\n")

        rn_scr = Path(database) / "GLUVAB_tree" / "rename.py"
        rn_lst = Path(database) / "GLUVAB_tree" / "rename_list"
        if rn_scr.exists():
            run_command(["python", str(rn_scr), str(tree1), str(rn_lst), str(small_fasta)],
                        "Renaming for Small Tree", logger)
        else:
            shutil.copy(tree1, small_fasta)

        if not small_tree_out.exists():
            logger.info("[Tree] Step 4.3: Running GLUVAB for Small Tree...")
            run_command(["perl", str(glu_scr), "--genomes_file_1", str(small_fasta),
                         "--file_prefix", sample_name, "--threads", str(threads)],
                        "GLUVAB Small Tree Calculation", logger, env_name=GLUVAB_ENV,
                        cwd=str(small_dir), check=False)

            if (small_dir / f"{sample_name}_Dice_Tree.newick").exists():
                shutil.copy(small_dir / f"{sample_name}_Dice_Tree.newick", small_tree_out)
                logger.info("[Tree] Small Tree Newick generated. Starting visualization...")
                r_script = Path(database) / "GLUVAB_tree" / "tree_plot.R"
                tsv_path = Path(database) / "GLUVAB_tree" / "taxonomy.tsv"
                if r_script.exists() and tsv_path.exists():
                    run_command(["Rscript", str(r_script), str(small_tree_out), str(tsv_path),
                                 sample_name, str(phy_dir / f"{sample_name}_tree")],
                                "ggtree visualization (using small tree)", logger,
                                env_name=GLUVAB_ENV, check=False)

    if big_tree_out.exists():
        logger.info("[Tree] [SKIP] Big tree already exists.")
    else:
        if not small_fasta.exists():
            logger.error("[Tree] [ERROR] Cannot build big tree: small_tree_in.fasta missing.")
            return phy_dir

        logger.info("[Tree] Step 4.4: De-duplicating against backbone for Big Tree...")
        bg_ref = Path(database) / "GLUVAB_tree" / "tree_ref.fasta"
        backbone_ids = set()
        if bg_ref.exists():
            with open(bg_ref, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        backbone_ids.add(line.strip()[1:].split()[0])

        big_fasta_in = phy_dir / "big_tree_in.fasta"
        with open(big_fasta_in, 'w') as out_f:
            write_flag = True
            with open(small_fasta, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        norm_id = line.strip()[1:].split()[0]
                        write_flag = norm_id not in backbone_ids
                        if not write_flag:
                            logger.debug(f"[Tree] Removing ID '{norm_id}' from sample set.")
                    if write_flag:
                        out_f.write(line.rstrip() + "\n")
            if bg_ref.exists():
                with open(bg_ref, 'r') as f:
                    for line in f:
                        out_f.write(line.rstrip() + "\n")

        logger.info("[Tree] Step 4.5: Running GLUVAB for Big Tree...")
        run_command(["perl", str(glu_scr), "--genomes_file_1", str(big_fasta_in),
                     "--file_prefix", sample_name, "--threads", str(threads)],
                    "GLUVAB Big Tree Calculation", logger, env_name=GLUVAB_ENV,
                    cwd=str(big_dir), check=False)

        if (big_dir / f"{sample_name}_Dice_Tree.newick").exists():
            shutil.copy(big_dir / f"{sample_name}_Dice_Tree.newick", big_tree_out)
            logger.info(f"[Tree] Big Tree Success → {big_tree_out}")


    if small_tree_out.exists() and big_tree_out.exists():
        shutil.rmtree(ref_dir, ignore_errors=True)
        tree1.unlink(missing_ok=True)

    return phy_dir


def run_characterization(input_fasta, threads, database, output_dir):
    input_fasta, database, output_dir = Path(input_fasta).resolve(), Path(database).resolve(), Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    logger, log_file = setup_logger(output_dir)
    logger.info("="*60 + "\nJMPhage - Genomic characterization (single phage)\n" + "="*60)

    start_time = datetime.now()
    try:
        check_file_exists(input_fasta, "Input FASTA", logger)
        sample_name = get_sample_name(input_fasta)

        ntw = build_shared_network(input_fasta, sample_name, output_dir, database, threads, logger)
        if ntw:
            run_collinearity_analysis(sample_name, ntw, output_dir, database, threads, logger)
            run_ani_analysis(sample_name, ntw, output_dir, database, threads, logger, input_fasta)
            run_phylogenetic_tree(sample_name, ntw, output_dir, database, threads, logger, input_fasta)

        total_time = (datetime.now() - start_time).total_seconds()
        logger.info(f"\n{'='*60}\n[OK] Pipeline completed in {total_time/60:.2f} min\n{'='*60}")

    except Exception as e:
        logger.error(f"Pipeline failed: {e}\n{traceback.format_exc()}")
        sys.exit(1)


#############multiple phages#####################
JOINT_TREE_TOPN = 30
PER_PHAGE_SMALL_TREE_TOPN = 50

def _pharokka_one(phage_fasta, phage_name, pharokka_out_root, threads, database, logger):

    out_dir = pharokka_out_root / phage_name
    faa = out_dir / "phanotate.faa"
    gbk = out_dir / f"{phage_name}.gbk"

    if faa.exists() and gbk.exists() and faa.stat().st_size > 0:
        logger.info(f"  [SKIP] Pharokka for {phage_name} already done.")
        return faa, gbk

    out_dir.parent.mkdir(parents=True, exist_ok=True)
    pharokka_db = Path(database) / "pharokka_db"

    if not pharokka_db.exists():
        logger.error(f"[ERROR] Pharokka database not found: {pharokka_db}")
        return None, None

    cmd = [
        'pharokka.py',
        '-i', str(Path(phage_fasta).resolve()),
        '-o', str(out_dir),
        '-p', phage_name,
        '-l', phage_name,
        '-t', str(threads),
        '-f',
        '-d', str(pharokka_db)
    ]
    run_command(cmd, f"Pharokka: {phage_name}", logger, quiet=True)

    if not faa.exists():
        logger.error(f"[ERROR] Pharokka failed for {phage_name}")
        return None, None

    cleaned = []
    with open(faa, 'r') as f:
        for line in f:
            cleaned.append(line.split(' ')[0] + '\n' if line.startswith('>') else line)
    with open(faa, 'w') as f:
        f.writelines(cleaned)

    return faa, gbk


def build_joint_network(phages, summary_dir, database, threads, logger):

    logger.info(f"\n{'#'*60}\n# Joint Step 1: Shared vConTACT2 Network (N={len(phages)})\n{'#'*60}\n")

    ntw_dir = summary_dir / "shared_network"
    pharokka_root = ntw_dir / "pharokka"
    vcontact2_output = ntw_dir / "vcontact2_output"
    c1_ntw = vcontact2_output / "c1.ntw"

    ntw_dir.mkdir(parents=True, exist_ok=True)
    pharokka_root.mkdir(parents=True, exist_ok=True)

    pharokka_map = {}

    if c1_ntw.exists():
        logger.info("[INFO] Existing joint vConTACT2 result found. Skipping vConTACT2, but verifying Pharokka outputs...")

        for p in phages:
            faa, gbk = _pharokka_one(p['fasta'], p['name'], pharokka_root, threads, database, logger)
            if faa is None:
                return None, None
            pharokka_map[p['name']] = {'faa': faa, 'gbk': gbk}
        return ntw_dir, pharokka_map

    logger.info(f"[1/4] Running Pharokka for {len(phages)} phage(s)...")
    for p in phages:
        faa, gbk = _pharokka_one(p['fasta'], p['name'], pharokka_root, threads, database, logger)
        if faa is None:
            return None, None
        pharokka_map[p['name']] = {'faa': faa, 'gbk': gbk}

    logger.info(f"[2/4] Generating gene-to-genome mappings...")
    tmp_maps = []
    for p in phages:
        name = p['name']
        out_csv = ntw_dir / f"_map_{name}.csv"
        write_gene2genome_csv(pharokka_map[name]['faa'], name, out_csv)
        tmp_maps.append(out_csv)

    logger.info(f"[3/4] Merging with reference database...")
    coline_db = Path(database) / "coline_profile"
    ref_map, ref_faa = coline_db / "out_map.csv", coline_db / "ref.faa"
    if not ref_map.exists() or not ref_faa.exists():
        logger.error("[ERROR] Missing database files (out_map.csv / ref.faa)")
        return None, None

    final_map_csv = ntw_dir / "out_map.csv"
    all_faa = ntw_dir / "all.faa"

    with open(final_map_csv, 'w') as out_f:
        for idx, m in enumerate(tmp_maps):
            with open(m, 'r') as in_f:
                content = in_f.read()
                
            if content and not content.endswith('\n'):
                content += '\n'
                
            lines = content.splitlines(keepends=True)
            
            if idx == 0:
                out_f.writelines(lines)
            else:
                out_f.writelines(lines[1:])
        
        with open(ref_map, 'r') as in_f:
            ref_lines = in_f.readlines()
            if ref_lines and ref_lines[0].lower().startswith(('protein_id', 'gene', '"protein')):
                out_f.writelines(ref_lines[1:])
            else:
                out_f.writelines(ref_lines)

    for m in tmp_maps:
        m.unlink(missing_ok=True)

    with open(all_faa, 'wb') as out_f:
        for p in phages:
            with open(pharokka_map[p['name']]['faa'], 'rb') as in_f:
                shutil.copyfileobj(in_f, out_f)
        with open(ref_faa, 'rb') as in_f:
            shutil.copyfileobj(in_f, out_f)

    logger.info(f"[4/4] Running joint vConTACT2...")
    fix_vcontact2_compatibility(logger)
    c1_jar = Path(database) / 'cluster_one-1.0.jar'

    cmd = [
        'vcontact2', '--raw-proteins', str(all_faa), '--rel-mode', 'Diamond',
        '--proteins-fp', str(final_map_csv), '--db', 'None',
        '--pcs-mode', 'MCL', '--vcs-mode', 'ClusterONE',
        '--c1-bin', str(c1_jar), '--output-dir', str(vcontact2_output),
        '-t', str(threads)
    ]
    run_command(cmd, "Joint vConTACT2", logger)

    if not c1_ntw.exists():
        logger.error("[ERROR] vConTACT2 did not produce c1.ntw")
        return None, None

    return ntw_dir, pharokka_map


def _query_neighbors(c1_ntw, sample_name, top_n):

    pairs = []
    with open(c1_ntw, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            a, b, w = parts[0], parts[1], parts[2]
            try:
                w = float(w)
            except ValueError:
                continue
            if a == sample_name:
                pairs.append((b, w))
            elif b == sample_name:
                pairs.append((a, w))
    pairs.sort(key=lambda x: x[1], reverse=True)
    seen = set()
    result = []
    for n, _ in pairs:
        if n in seen or n == sample_name:
            continue
        seen.add(n)
        result.append(n)
        if len(result) >= top_n:
            break
    return result


def run_per_phage_collinearity_joint(phages, ntw_dir, pharokka_map, base_out, database, threads, logger):

    logger.info(f"\n{'#'*60}\n# Joint Step 2: Per-phage Collinearity (from shared c1.ntw)\n{'#'*60}\n")

    c1_ntw = ntw_dir / "vcontact2_output" / "c1.ntw"
    if not c1_ntw.exists():
        logger.error("[ERROR] c1.ntw not found, cannot do collinearity")
        return

    for p in phages:
        name = p['name']
        phage_out = base_out / name
        coline_dir = phage_out / "6.collinearity_analysis"
        output_html = coline_dir / f"{name}.html"

        if output_html.exists():
            logger.info(f"  [SKIP] Collinearity for {name} already done.")
            continue

        logger.info(f"\n  >>> Collinearity: {name}")

        fake_ntw = phage_out / "5.shared_network"
        fake_pharokka = fake_ntw / "pharokka_output"
        fake_vct = fake_ntw / "vcontact2_output"
        fake_pharokka.mkdir(parents=True, exist_ok=True)
        fake_vct.mkdir(parents=True, exist_ok=True)

        src_faa = pharokka_map[name]['faa']
        src_gbk = pharokka_map[name]['gbk']
        dst_faa = fake_pharokka / "phanotate.faa"
        dst_gbk = fake_pharokka / f"{name}.gbk"
        if not dst_faa.exists():
            shutil.copy(src_faa, dst_faa)
        if not dst_gbk.exists() and src_gbk.exists():
            shutil.copy(src_gbk, dst_gbk)

        dst_ntw = fake_vct / "c1.ntw"
        if not dst_ntw.exists():
            shutil.copy(c1_ntw, dst_ntw)

        try:
            run_collinearity_analysis(name, fake_ntw, phage_out, database, threads, logger)
        except Exception as e:
            logger.error(f"  [WARN] Collinearity failed for {name}: {e}")


def _parse_fasta_ids(fasta_path):
    ids = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                raw = line.strip()[1:].split()[0]
                ids.append(raw)
    return ids


def _normalize_id(raw_id):
    return re.sub(r'\..*', '', raw_id).strip()


def _append_fasta_dedupe(src_fasta, out_handle, seen_ids, logger=None, source_tag=""):

    added = 0
    skipped = 0
    write = False
    with open(src_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                raw = line.strip()[1:].split()[0]
                norm = _normalize_id(raw)
                if norm in seen_ids:
                    write = False
                    skipped += 1
                    if logger:
                        logger.debug(f"  [dedupe:{source_tag}] skip {raw}")
                    continue
                seen_ids.add(norm)
                out_handle.write(f">{norm}\n")
                write = True
                added += 1
            elif write:
                out_handle.write(line)
    return added, skipped


def run_per_phage_ani_joint(phages, ntw_dir, base_out, database, threads, logger):
    logger.info(f"\n{'#'*60}\n# Joint Step 3: Per-phage ANI (from shared c1.ntw)\n{'#'*60}\n")

    c1_ntw = ntw_dir / "vcontact2_output" / "c1.ntw"
    if not c1_ntw.exists():
        logger.error("[ERROR] c1.ntw not found, cannot do ANI")
        return

    for p in phages:
        name = p['name']
        phage_out = base_out / name

        viridic_out_dir = phage_out / "7.ANI_analysis" / "ANI_result" / "04_VIRIDIC_out"
        if (viridic_out_dir / "Heatmap.PDF").exists() or (viridic_out_dir / "sim_MA_genCol.csv").exists():
            logger.info(f"  [SKIP] ANI for {name} already done.")
            continue

        logger.info(f"\n  >>> ANI: {name}")

        fake_ntw = phage_out / "5.shared_network"
        fake_vct = fake_ntw / "vcontact2_output"
        fake_vct.mkdir(parents=True, exist_ok=True)
        dst_ntw = fake_vct / "c1.ntw"
        if not dst_ntw.exists():
            shutil.copy(c1_ntw, dst_ntw)

        try:
            run_ani_analysis(name, fake_ntw, phage_out, database, threads, logger, p['fasta'])
        except Exception as e:
            logger.error(f"  [WARN] ANI failed for {name}: {e}")

def run_joint_ani(phages, ntw_dir, summary_dir, database, threads, logger):

    logger.info(f"\n{'#'*60}\n# Joint Step 3b: Joint ANI Analysis (VIRIDIC)\n{'#'*60}\n")

    VIRIDIC_ENV = "jmp_VIRIDIC"
    c1_ntw = ntw_dir / "vcontact2_output" / "c1.ntw"
    if not c1_ntw.exists():
        logger.error("[ERROR] c1.ntw not found")
        return None

    ani_dir = (summary_dir / "ANI_analysis").resolve()
    ani_fasta_dir = ani_dir / "ANI_fasta"
    ani_res_dir = ani_dir / "ANI_result"
    merged_fasta = ani_dir / "ANI.fasta"
    viridic_out_dir = ani_res_dir / "04_VIRIDIC_out"
    heatmap_pdf = viridic_out_dir / "Heatmap.PDF"
    result_csv = viridic_out_dir / "sim_MA_genCol.csv"

    if heatmap_pdf.exists() or result_csv.exists():
        logger.info(f"[INFO] Skip: Joint VIRIDIC results already exist.")
        return ani_res_dir

    ani_fasta_dir.mkdir(parents=True, exist_ok=True)
    ani_res_dir.mkdir(parents=True, exist_ok=True)

    input_names = {p['name'] for p in phages}
    neighbors_union = set()
    for p in phages:
        for nb in _query_neighbors(c1_ntw, p['name'], 20):
            if nb not in input_names:
                neighbors_union.add(nb)

    logger.info(f"  Inputs: {len(phages)}, Unique neighbors: {len(neighbors_union)}")

    ref_src = Path(database) / "coline_profile" / "ref_fasta"
    count = 0

    for p in phages:
        dst = ani_fasta_dir / f"{p['name']}.fasta"
        if not dst.exists():
            shutil.copy(p['fasta'], dst)
        count += 1

    for nb in sorted(neighbors_union):
        src = ref_src / f"{nb}.fasta"
        if src.exists():
            shutil.copy(src, ani_fasta_dir / f"{nb}.fasta")
            count += 1

    logger.info(f"  Total sequences for joint ANI: {count}")

    if count < 2:
        logger.warning("  [WARN] Less than 2 sequences — skipping joint ANI.")
        return None

    with open(merged_fasta, 'w') as out_f:
        for f in sorted(ani_fasta_dir.glob("*.fasta")):
            out_f.write(f.read_text().strip() + "\n")

    logger.info("  Running joint VIRIDIC...")
    vir_src = Path(database) / "VIRIDIC"
    vir_tmp = ani_dir / "VIRIDIC_tool_tmp"
    if vir_tmp.exists():
        shutil.rmtree(vir_tmp)
    shutil.copytree(vir_src, vir_tmp)

    script_dir = (vir_tmp / "stand_alone" / "viridic_scripts").resolve()
    cmd = ["Rscript", "00_viridic_master.R",
           f"projdir={ani_res_dir}", f"in={merged_fasta}", f"ncor={threads}"]
    run_command(cmd, "Joint VIRIDIC R-calculation", logger,
                env_name=VIRIDIC_ENV, cwd=str(script_dir), check=True)

    shutil.rmtree(vir_tmp, ignore_errors=True)
    shutil.rmtree(ani_fasta_dir, ignore_errors=True)

    if heatmap_pdf.exists() or result_csv.exists():
        logger.info(f"  Joint ANI completed → {ani_res_dir}")
    else:
        logger.warning("  [WARN] VIRIDIC did not produce expected output.")

    return ani_res_dir

def run_joint_phylogenetic_tree(phages, ntw_dir, summary_dir, database, threads, logger):
    logger.info(f"\n{'#'*60}\n# Joint Step 4: Joint Phylogenetic Tree (GLUVAB)\n{'#'*60}\n")

    GLUVAB_ENV = "jmp_GLUVAB"
    c1_ntw = ntw_dir / "vcontact2_output" / "c1.ntw"
    if not c1_ntw.exists():
        logger.error("[ERROR] c1.ntw not found")
        return None

    phy_dir = (summary_dir / "phylogenetic_tree").resolve()
    small_dir = phy_dir / "small_tree"
    big_dir = phy_dir / "big_tree"

    JOINT_PREFIX = "joint"
    small_tree_out = small_dir / f"{JOINT_PREFIX}_Tree_With_node_IDs.newick"
    big_tree_out = big_dir / f"{JOINT_PREFIX}_Tree_With_node_IDs.newick"

    if small_tree_out.exists() and big_tree_out.exists():
        logger.info(f"[INFO] Skip: Joint trees already exist.")
        return phy_dir

    for d in [phy_dir, small_dir, big_dir]:
        d.mkdir(parents=True, exist_ok=True)

    input_names = {p['name'] for p in phages}
    neighbors_union = set()
    for p in phages:
        for nb in _query_neighbors(c1_ntw, p['name'], JOINT_TREE_TOPN):
            if nb in input_names:
                continue
            neighbors_union.add(nb)

    logger.info(f"  Inputs: {len(phages)}, Unique neighbors (vs inputs): {len(neighbors_union)}")

    ref_src = Path(database) / "coline_profile" / "ref_fasta"
    tree1 = phy_dir / "tree1.fasta"
    seen_ids = set()
    n_in, n_nb, n_skip = 0, 0, 0
    missing = []

    with open(tree1, 'w') as out_f:
        for p in phages:
            added, skipped = _append_fasta_dedupe(p['fasta'], out_f, seen_ids, logger, "input")
            n_in += added
            n_skip += skipped
        for nb in sorted(neighbors_union):
            ref_fa = ref_src / f"{nb}.fasta"
            if not ref_fa.exists():
                missing.append(nb)
                continue
            added, skipped = _append_fasta_dedupe(ref_fa, out_f, seen_ids, logger, "neighbor")
            n_nb += added
            n_skip += skipped

    if missing:
        logger.warning(f"  [WARN] {len(missing)} neighbor reference FASTA(s) missing in DB.")
    logger.info(f"  Small tree sequences: {n_in + n_nb} (inputs={n_in}, neighbors={n_nb}, deduped={n_skip})")

    if (n_in + n_nb) < 2:
        logger.warning("  [WARN] Less than 2 sequences after dedupe — skipping tree.")
        return None

    small_fasta = phy_dir / "small_tree_in.fasta"
    rn_scr = Path(database) / "GLUVAB_tree" / "rename.py"
    rn_lst = Path(database) / "GLUVAB_tree" / "rename_list"
    if rn_scr.exists():
        run_command(["python", str(rn_scr), str(tree1), str(rn_lst), str(small_fasta)],
                    "Renaming for Joint Small Tree", logger)
    else:
        shutil.copy(tree1, small_fasta)

    logger.info("  Running GLUVAB for Joint Small Tree...")
    glu_scr = (Path(database) / "GLUVAB_tree" / "GLUVAB_v0.6.pl").resolve()
    run_command(["perl", str(glu_scr), "--genomes_file_1", str(small_fasta),
                 "--file_prefix", JOINT_PREFIX, "--threads", str(threads)],
                "GLUVAB Joint Small Tree", logger,
                env_name=GLUVAB_ENV, cwd=str(small_dir), check=False)

    if (small_dir / f"{JOINT_PREFIX}_Dice_Tree.newick").exists():
        shutil.copy(small_dir / f"{JOINT_PREFIX}_Dice_Tree.newick", small_tree_out)
        logger.info("  Joint Small Tree generated. Visualizing...")
        r_script = Path(database) / "GLUVAB_tree" / "tree_plot.R"
        tsv_path = Path(database) / "GLUVAB_tree" / "taxonomy.tsv"
        if r_script.exists() and tsv_path.exists():
            run_command(["Rscript", str(r_script), str(small_tree_out), str(tsv_path),
                         JOINT_PREFIX, str(phy_dir / f"{JOINT_PREFIX}_tree")],
                        "Joint ggtree visualization", logger,
                        env_name=GLUVAB_ENV, check=False)

    logger.info("  De-duplicating against backbone for Joint Big Tree...")
    bg_ref = Path(database) / "GLUVAB_tree" / "tree_ref.fasta"
    backbone_ids = set()
    if bg_ref.exists():
        with open(bg_ref, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    backbone_ids.add(line.strip()[1:].split()[0])

    big_fasta_in = phy_dir / "big_tree_in.fasta"
    removed = 0
    with open(big_fasta_in, 'w') as out_f:
        write_flag = True
        with open(small_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    norm_id = line.strip()[1:].split()[0]
                    write_flag = norm_id not in backbone_ids
                    if not write_flag:
                        removed += 1
                        logger.debug(f"  [Tree] Removing ID '{norm_id}' (collides with backbone).")
                if write_flag:
                    out_f.write(line.rstrip() + "\n")
        if bg_ref.exists():
            with open(bg_ref, 'r') as f:
                for line in f:
                    out_f.write(line.rstrip() + "\n")
    if removed:
        logger.info(f"  Removed {removed} sequences from sample set (backbone takes priority).")

    logger.info("  Running GLUVAB for Joint Big Tree...")
    run_command(["perl", str(glu_scr), "--genomes_file_1", str(big_fasta_in),
                 "--file_prefix", JOINT_PREFIX, "--threads", str(threads)],
                "GLUVAB Joint Big Tree", logger,
                env_name=GLUVAB_ENV, cwd=str(big_dir), check=False)

    if (big_dir / f"{JOINT_PREFIX}_Dice_Tree.newick").exists():
        shutil.copy(big_dir / f"{JOINT_PREFIX}_Dice_Tree.newick", big_tree_out)
        logger.info(f"  Joint Big Tree generated → {big_tree_out}")

    tree1.unlink(missing_ok=True)

    return phy_dir

def run_per_phage_small_tree_joint(phages, ntw_dir, base_out, database, threads, logger):

    logger.info(f"\n{'#'*60}\n# Joint Step 5: Per-phage Small Tree (from shared c1.ntw)\n{'#'*60}\n")

    GLUVAB_ENV = "jmp_GLUVAB"
    c1_ntw = ntw_dir / "vcontact2_output" / "c1.ntw"
    if not c1_ntw.exists():
        logger.error("[ERROR] c1.ntw not found, cannot build per-phage small trees")
        return

    glu_scr  = (Path(database) / "GLUVAB_tree" / "GLUVAB_v0.6.pl").resolve()
    ref_src  = Path(database) / "coline_profile" / "ref_fasta"
    rn_scr   = Path(database) / "GLUVAB_tree" / "rename.py"
    rn_lst   = Path(database) / "GLUVAB_tree" / "rename_list"
    r_script = Path(database) / "GLUVAB_tree" / "tree_plot.R"
    tsv_path = Path(database) / "GLUVAB_tree" / "taxonomy.tsv"

    input_names = {p['name'] for p in phages}

    for p in phages:
        name = p['name']
        phy_dir   = (base_out / name / "8.phylogenetic_tree").resolve()
        small_dir = phy_dir / "small_tree"
        ref_dir   = phy_dir / "ref"
        small_tree_out = small_dir / f"{name}_Tree_With_node_IDs.newick"

        if small_tree_out.exists():
            logger.info(f"  [SKIP] Small tree for {name} already done.")
            continue

        logger.info(f"\n  >>> Small tree: {name}")
        for d in (phy_dir, small_dir, ref_dir):
            d.mkdir(parents=True, exist_ok=True)

        neighbors = [nb for nb in _query_neighbors(c1_ntw, name,
                                                   PER_PHAGE_SMALL_TREE_TOPN + len(phages))
                     if nb not in input_names][:PER_PHAGE_SMALL_TREE_TOPN]

        tree1 = phy_dir / "tree1.fasta"
        seen_ids = set()
        n_in, n_nb, n_skip = 0, 0, 0
        missing = []
        with open(tree1, 'w') as out_f:
            added, skipped = _append_fasta_dedupe(p['fasta'], out_f, seen_ids, logger, "input")
            n_in += added; n_skip += skipped
            for nb in neighbors:
                ref_fa = ref_src / f"{nb}.fasta"
                if not ref_fa.exists():
                    missing.append(nb)
                    continue
                added, skipped = _append_fasta_dedupe(ref_fa, out_f, seen_ids, logger, "neighbor")
                n_nb += added; n_skip += skipped

        if missing:
            logger.warning(f"  [WARN] {name}: {len(missing)} neighbor reference FASTA(s) missing in DB.")
        logger.info(f"  {name}: small tree sequences = {n_in + n_nb} "
                    f"(input={n_in}, neighbors={n_nb}, deduped={n_skip})")

        if (n_in + n_nb) < 2:
            logger.warning(f"  [WARN] {name}: fewer than 2 sequences — skipping small tree.")
            tree1.unlink(missing_ok=True)
            continue

        small_fasta = phy_dir / "small_tree_in.fasta"
        if rn_scr.exists():
            run_command(["python", str(rn_scr), str(tree1), str(rn_lst), str(small_fasta)],
                        f"Rename small tree: {name}", logger, quiet=True)
        else:
            shutil.copy(tree1, small_fasta)

        run_command(["perl", str(glu_scr), "--genomes_file_1", str(small_fasta),
                     "--file_prefix", name, "--threads", str(threads)],
                    f"GLUVAB Small Tree: {name}", logger,
                    env_name=GLUVAB_ENV, cwd=str(small_dir), check=False)

        dice = small_dir / f"{name}_Dice_Tree.newick"
        if dice.exists():
            shutil.copy(dice, small_tree_out)
            logger.info(f"  {name}: small tree generated. Visualizing...")
            if r_script.exists() and tsv_path.exists():
                run_command(["Rscript", str(r_script), str(small_tree_out), str(tsv_path),
                             name, str(phy_dir / f"{name}_tree")],
                            f"ggtree visualization: {name}", logger,
                            env_name=GLUVAB_ENV, check=False)
        else:
            logger.warning(f"  [WARN] {name}: GLUVAB did not produce a small tree newick.")

        shutil.rmtree(ref_dir, ignore_errors=True)
        tree1.unlink(missing_ok=True)

def run_joint_characterization(manifest_path):

    with open(manifest_path, 'r') as f:
        manifest = json.load(f)

    base_out = Path(manifest['base_output']).resolve()
    database = Path(manifest['database']).resolve()
    threads = int(manifest['threads'])
    phages = manifest['phages']  # [{name, fasta}, ...]

    for p in phages:
        p['fasta'] = Path(p['fasta']).resolve()

    base_out.mkdir(parents=True, exist_ok=True)
    summary_dir = base_out / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    logger, log_file = setup_logger(summary_dir, "joint_characterization")
    logger.info("="*60)
    logger.info(f"JMPhage - Joint Characterization (N={len(phages)})")
    logger.info("="*60)
    logger.info(f"Inputs: {', '.join(p['name'] for p in phages)}")
    logger.info(f"Output: {base_out}")
    logger.info(f"Summary: {summary_dir}")

    start_time = datetime.now()
    try:

        ntw_dir, pharokka_map = build_joint_network(phages, summary_dir, database, threads, logger)
        if ntw_dir is None or pharokka_map is None:
            logger.error("[ABORT] Joint network failed.")
            sys.exit(1)

        run_per_phage_collinearity_joint(phages, ntw_dir, pharokka_map, base_out, database, threads, logger)

        run_per_phage_ani_joint(phages, ntw_dir, base_out, database, threads, logger)

        run_joint_ani(phages, ntw_dir, summary_dir, database, threads, logger)

        run_joint_phylogenetic_tree(phages, ntw_dir, summary_dir, database, threads, logger)

        run_per_phage_small_tree_joint(phages, ntw_dir, base_out, database, threads, logger)
        
        total = (datetime.now() - start_time).total_seconds()
        logger.info(f"\n{'='*60}")
        logger.info(f"[OK] Joint characterization done in {total/60:.2f} min")
        logger.info(f"{'='*60}")

    except Exception as e:
        logger.error(f"Joint characterization failed: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)



def main():
    parser = argparse.ArgumentParser(description='JMPhage Characterization (single + joint)')
    parser.add_argument('--joint', default=None,
                        help='Joint mode: path to a JSON manifest. If set, other args are ignored.')
    parser.add_argument('-i', '--input', default=None)
    parser.add_argument('-d', '--database', default=None)
    parser.add_argument('-t', '--threads', type=int, default=4)
    parser.add_argument('-o', '--output', default=None)
    args = parser.parse_args()

    if args.joint:
        run_joint_characterization(args.joint)
        return

    if not (args.input and args.database and args.output):
        parser.error("Single-phage mode requires -i, -d, and -o.")
    run_characterization(args.input, args.threads, args.database, args.output)


if __name__ == "__main__":
    main()
