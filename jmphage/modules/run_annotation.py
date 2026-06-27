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

ENV_NAME = "JM_annotation"

def setup_logger(output_dir, module_name="annotation"):
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
    file_handler.flush = lambda: file_handler.stream.flush()

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
        logger.info(f"[CHECK] {description} check passed: {filepath}")
    return True


def get_sample_name(input_fasta):
    filename = Path(input_fasta).name
    sample_name = filename.split('.')[0]
    return sample_name


def run_command(cmd, description="", logger=None, cwd=None):
    if logger:
        logger.info(f"\n{'='*60}\n{description}\n{'='*60}")

    env = os.environ.copy()

    current_bin = os.path.dirname(sys.executable)

    old_path = env.get("PATH", "")
    env["PATH"] = current_bin + os.pathsep + old_path

    blacklist = [
        "PERL5LIB", "PERL_LOCAL_LIB_ROOT", "PERL_MM_OPT", "PERL_MB_OPT",
        "PYTHONPATH", "PYTHONHOME", 
        "LD_LIBRARY_PATH" 
    ]
    
    for var in blacklist:
        if var in env:
            if logger:
                logger.debug(f"[Firewall] Masking potential conflict variable: {var}={env[var][:50]}...")
            del env[var]
    if logger:
        paths = env["PATH"].split(os.pathsep)
        logger.debug(f"Path Priority 1: {paths[0]}")
        if len(paths) > 1: logger.debug(f"Path Priority 2: {paths[1]}")
    if isinstance(cmd, str): cmd = cmd.split()
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, 
            text=True,
            cwd=cwd,      
            env=env,      
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
                logger.error(f"--- Execution Failure Diagnosis ---")
                logger.error(f"Failed Command: {' '.join(cmd)}")
                logger.error(f"[FAILED] Command exited with return code {process.returncode}")
            raise subprocess.CalledProcessError(process.returncode, cmd)
        if logger: logger.info(f"[OK] {description} completed")
        return process.returncode
    except subprocess.CalledProcessError as e:
        if logger:
            logger.error(f"Command execution failed.")
        raise


def run_pharokka(input_fasta, sample_name, pre_dir, threads, database, logger):
    combined_faa = pre_dir / "combined_output.faa"
    if combined_faa.exists() and combined_faa.stat().st_size > 0:
        logger.info(f"[SKIP] Pharokka output already exists: {combined_faa.name}")
        return combined_faa

    pharokka_output = pre_dir / "pharokka_output"
    pharokka_db = Path(database) / "pharokka_db"

    if not pharokka_db.exists():
        logger.error(f"[ERROR] Pharokka database not found: {pharokka_db}")
        raise FileNotFoundError(f"Pharokka database not found: {pharokka_db}")

    cmd = [
        'pharokka.py',
        '-i', input_fasta,
        '-o', str(pharokka_output),
        '-p', sample_name,
        '-l', sample_name,
        '-t', str(threads),
        '-f',
        '-d', str(pharokka_db)
    ]
    
    run_command(cmd, "Step 1: Running Pharokka for ORF prediction...", logger)
    
    faa_file = pharokka_output / "phanotate.faa"
    if faa_file.exists():
        logger.info(f"[OK] Pharokka completed")
        shutil.copy(faa_file, combined_faa)
        return combined_faa
    else:
        logger.error("[ERROR] Pharokka failed, phanotate.faa not generated")
        raise FileNotFoundError(f"Pharokka output not found: {faa_file}")


def clean_fasta_headers(faa_file, logger):

    logger.info("\nCleaning FASTA headers for protein sequences...")
    try:
        with open(faa_file, 'r') as f:
            lines = f.readlines()
        
        cleaned_lines = []
        for line in lines:
            if line.startswith('>'):
                cleaned_lines.append(line.split()[0] + '\n')
            else:
                cleaned_lines.append(line)
        
        with open(faa_file, 'w') as f:
            f.writelines(cleaned_lines)
        logger.debug(f"[OK] Header cleaned: {faa_file}")
    except Exception as e:
        logger.warning(f"[WARNING] Failed to clean FASTA header: {e}")


def run_vog_annotation(combined_faa, sample_name, anno_dir, database, logger):
    top_vog = anno_dir / f"{sample_name}_top_vog.tsv"
    if top_vog.exists() and top_vog.stat().st_size > 0:
        logger.info(f"[SKIP] VOG annotation already exists: {top_vog.name}")
        return top_vog

    vog_db = Path(database) / "VOG_HMM" / "vog_db.dmnd"
    if not vog_db.exists():
        logger.warning(f"[WARNING] VOG database not found: {vog_db}")
        return None
    
    vog_output = anno_dir / f"{sample_name}_vog.tsv"
    
    cmd = [
        'diamond', 'blastp',
        '-d', str(vog_db),
        '-q', str(combined_faa),
        '-o', str(vog_output),
        '-e', '1e-5',
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 
                    'evalue', 'bitscore', 'stitle'
    ]
    
    run_command(cmd, "Step 2: Running VOG annotation with Diamond...", logger)
    
    if vog_output.exists():
        top_vog = anno_dir / f"{sample_name}_top_vog.tsv"
        logger.info("Extracting best hits for each protein...")
        try:
            seen = set()
            with open(vog_output, 'r') as infile, open(top_vog, 'w') as outfile:
                for line in infile:
                    fields = line.strip().split('\t')
                    if len(fields) >= 7:
                        qseqid = fields[0]
                        stitle = fields[6]
                        if qseqid not in seen:
                            outfile.write(f"{qseqid}\t{stitle}\n")
                            seen.add(qseqid)
            logger.info(f"[OK] VOG annotation completed")
            return top_vog
        except Exception as e:
            logger.error(f"[ERROR] Failed to process VOG results: {e}")
            return None
    else:
        logger.error("[ERROR] VOG annotation failed")
        return None



def run_hmmscan(faa_file, sample_name, plot_dir, database, threads, logger):

    hmm_db = Path(database) / "VOG_HMM" / "phage_all.hmm"
    if not hmm_db.exists():
        logger.warning(f"[WARNING] HMM database not found: {hmm_db}")
        return None
    
    domtblout = plot_dir / "results.domtblout"
    
    cmd = [
        'hmmscan',
        '-E', '1e-5',
        '--cpu', str(threads),
        '--domtblout', str(domtblout),
        str(hmm_db),
        str(faa_file)
    ]
    
    run_command(cmd, "Step 4: Running HMMscan for domain annotation...", logger)
    
    if not domtblout.exists():
        logger.error("[ERROR] HMMscan failed to generate results")
        return None
    return domtblout


def process_hmm_results(domtblout, plot_dir, database, logger):
    logger.info("\nProcessing HMM results...")
    try:
        protein_hmm_result = plot_dir / "protein_hmm_result"
        best = {}
        with open(domtblout, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.split()
                if len(fields) >= 13:
                    target = fields[0]
                    query = fields[3]
                    if query not in best:
                        best[query] = target
        
        with open(protein_hmm_result, 'w') as f:
            for query in sorted(best.keys()):
                f.write(f"{query}\t{best[query]}\n")
        
        hmm_type_file = plot_dir / "hmm_type"
        phage_hmm_type = Path(database) / "VOG_HMM" / "phage_hmm_type.tsv"
        
        if phage_hmm_type.exists():
            with open(protein_hmm_result, 'r') as infile, open(hmm_type_file, 'w') as outfile:
                for line in infile:
                    hmm_id = line.strip().split('\t')[1]
                    with open(phage_hmm_type, 'r') as type_file:
                        for type_line in type_file:
                            if type_line.startswith(hmm_id):
                                outfile.write(type_line)
                                break
        
        protein_hmm_type = plot_dir / "protein_hmm_type"
        if hmm_type_file.exists():
            hmm_type_dict = {}
            with open(hmm_type_file, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        hmm_type_dict[fields[0]] = '\t'.join(fields[1:])
            
            with open(protein_hmm_result, 'r') as infile, open(protein_hmm_type, 'w') as outfile:
                for line in infile:
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        protein_id = fields[0]
                        hmm_id = fields[1]
                        type_info = hmm_type_dict.get(hmm_id, '')
                        outfile.write(f"{protein_id}\t{hmm_id}\t{type_info}\n")
        
        return protein_hmm_result
    except Exception as e:
        logger.error(f"[ERROR] Failed to process HMM results: {e}")
        return None


def create_genome_features_from_gff(gff_file, protein_hmm_type, plot_dir, logger):
    logger.info("\nGenerating genome feature file (from pharokka GFF)...")
    try:
        protein_start_end = plot_dir / "protein_start_end"
        proteins = []

        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                if line.startswith('>'):
                    break
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 9:
                    continue
                feature_type = fields[2]
                if feature_type != 'CDS':
                    continue
                try:
                    start = int(fields[3])
                    end = int(fields[4])
                except ValueError:
                    continue
                strand = 1 if fields[6] == '+' else -1

                attrs = fields[8]
                protein_id = None
                for kv in attrs.split(';'):
                    kv = kv.strip()
                    if kv.startswith('ID='):
                        protein_id = kv[3:]
                        break
                if protein_id is None:
                    continue

                proteins.append({
                    'id': protein_id,
                    'start': start,
                    'end': end,
                    'strand': strand
                })

        with open(protein_start_end, 'w') as f:
            for p in proteins:
                f.write(f"{p['id']}\t{p['start']}\t{p['end']}\t{p['strand']}\n")

        genome_features = plot_dir / "genome_features_classified.tsv"
        hmm_dict = {}
        if protein_hmm_type.exists():
            with open(protein_hmm_type, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= 3:
                        hmm_dict[fields[0]] = '\t'.join(fields[2:])

        with open(protein_start_end, 'r') as infile, open(genome_features, 'w') as outfile:
            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    protein_id = fields[0]
                    classification = hmm_dict.get(protein_id, 'Unclassified')
                    outfile.write(f"{line.strip()}\t{classification}\n")

        logger.info(f"[OK] Genome feature file generated: {genome_features}")
        return genome_features
    except Exception as e:
        logger.error(f"[ERROR] Failed to generate genome feature file: {e}")
        return None


def extract_hmm_annotations(protein_hmm_result, plot_dir, database, logger):
    logger.info("\nExtracting and integrating HMM functional annotations...")
    try:
        vog_annotations = Path(database) / "VOG_HMM" / "vog.annotations.tsv"
        if not vog_annotations.exists():
            logger.warning(f"[WARNING] VOG annotation file not found: {vog_annotations}")
            return None

        vog_dict = {}
        with open(vog_annotations, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    vog_dict[fields[0]] = fields[-1].replace('REFSEQ', '').strip()

        combined_output = plot_dir / "protein_function_summary.tsv"
        temp_csv = plot_dir / "hmm_annotation.csv" 

        with open(protein_hmm_result, 'r') as infile, \
             open(combined_output, 'w') as f_out, \
             open(temp_csv, 'w') as f_temp:

            f_out.write("Protein_ID\tHMM_ID\tFunctional_Annotation\n")
            
            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    protein_id = fields[0]
                    hmm_id = fields[1]
                    annotation = vog_dict.get(hmm_id, 'Unknown')
                    f_out.write(f"{protein_id}\t{hmm_id}\t{annotation}\n")
                    f_temp.write(f"{hmm_id}\t{annotation}\n")
        
        logger.info(f"[OK] Summary annotation file generated: {combined_output}")
        return combined_output
        
    except Exception as e:
        logger.error(f"[ERROR] Failed to extract HMM annotations: {e}")
        return None


def run_genome_plot(plot_dir, database, logger, sample_fasta):

    plot_script = Path(database).resolve() / "VOG_HMM" / "GC_genome_plot.py"
    
    if not plot_script.exists():
        logger.warning(f"[WARNING] Plotting script not found: {plot_script}")
        return None

    cmd = [
        'python', 
        str(plot_script), 
        str(Path(plot_dir).resolve()), 
        str(Path(sample_fasta).resolve())
    ]
    
    run_command(cmd, "Step 5: Generating genome plot...", logger, cwd=str(plot_dir))
    plot_files = list(Path(plot_dir).glob("*.pdf"))
    if plot_files:
        logger.info(f"[OK] Genome visualization completed, generated files: {[f.name for f in plot_files]}")
        return plot_files
    else:
        logger.warning("[WARNING] Process ended but no PDF files found in the directory")
        return None


def run_annotation(input_fasta, threads, database, output_dir):

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger, log_file = setup_logger(output_dir, "annotation")
    
    logger.info("="*60)
    logger.info("JMPhage - Annotation Module")
    logger.info("="*60)
    logger.info(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Log File: {log_file}")
    logger.info("")
    
    try:
        logger.info("[Checking Input Parameters]")
        check_file_exists(input_fasta, "Input FASTA file", logger)
        if not os.path.isdir(database):
            logger.error(f"[ERROR] Database directory not found: {database}")
            sys.exit(1)
        
        sample_name = get_sample_name(input_fasta)
        logger.info(f"\nSample Name: {sample_name}")
        logger.info(f"Output Dir: {output_dir}")
        logger.info(f"Threads: {threads}")
        
        input_path = Path(input_fasta).resolve()
        expected_file = input_path
        pre_dir = output_dir / "2.ORF_prediction"
        anno_dir = output_dir / "3.function_annotation"
        plot_dir = output_dir / "4.genome_visualization"
        
        pre_dir.mkdir(parents=True, exist_ok=True)
        anno_dir.mkdir(parents=True, exist_ok=True)
        plot_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"\n{'#'*60}")
        logger.info(f"# Starting annotation Pipeline")
        logger.info(f"{'#'*60}\n")
        
        pipeline_start = datetime.now()
        combined_faa = run_pharokka(str(expected_file), sample_name, pre_dir, threads, database, logger)
        clean_fasta_headers(combined_faa, logger)
        top_vog = run_vog_annotation(combined_faa, sample_name, anno_dir, database, logger)
        logger.info(f"\n{'='*60}")
        logger.info("[Genome Visualization Pipeline]")
        logger.info(f"{'='*60}")
        
        pharokka_faa = pre_dir / "pharokka_output" / "phanotate.faa"
        pharokka_gff = pre_dir / "pharokka_output" / f"{sample_name}.gff"
        
        if pharokka_faa.exists(): shutil.copy(pharokka_faa, plot_dir / f"{sample_name}.faa")
        if pharokka_gff.exists(): shutil.copy(pharokka_gff, plot_dir / f"{sample_name}.gff")
        shutil.copy(expected_file, plot_dir / "phage_genome.fasta")
        
        plot_faa = plot_dir / f"{sample_name}.faa"

        final_summary = plot_dir / "protein_function_summary.tsv"
        final_features = plot_dir / "genome_features_classified.tsv"
        plot_files = list(plot_dir.glob("*.pdf"))
        
        if plot_files:
            logger.info(f"[SKIP] Final genome plots (*.pdf) already exist. Skipping the entire visualization pipeline.")
        elif final_summary.exists() and final_features.exists():
            logger.info("[SKIP] HMM annotations already exist. Skipping HMMscan and generating plot directly...")
            run_genome_plot(plot_dir, database, logger, input_fasta)
        else:
            if plot_faa.exists():
                domtblout = run_hmmscan(plot_faa, sample_name, plot_dir, database, threads, logger)
                if domtblout:
                    protein_hmm_result = process_hmm_results(domtblout, plot_dir, database, logger)
                    if protein_hmm_result:
                        gff_file = plot_dir / f"{sample_name}.gff"
                        protein_hmm_type = plot_dir / "protein_hmm_type"
                        
                        if gff_file.exists():
                            create_genome_features_from_gff(gff_file, protein_hmm_type, plot_dir, logger)
                        
                        extract_hmm_annotations(protein_hmm_result, plot_dir, database, logger)
                        run_genome_plot(plot_dir, database, logger, input_fasta)
                        
                        logger.info("\nCleaning up intermediate files...")
                        files_to_delete = [
                            "results.domtblout",       
                            "protein_hmm_result",       
                            "hmm_type",                
                            "protein_hmm_type",        
                            "protein_start_end",       
                            "hmm_annotation.csv",      
                            f"{sample_name}.gff",      
                            f"{sample_name}.faa",      
                            "phage_genome.fasta"       
                        ]
                        
                        deleted_count = 0
                        for fname in files_to_delete:
                            f_path = plot_dir / fname
                            try:
                                if f_path.exists():
                                    f_path.unlink()
                                    deleted_count += 1
                            except Exception:
                                pass
                        
                        logger.info(f"[OK] {deleted_count} intermediate files removed.")
                        logger.info("Retained: genome_features_classified.tsv, protein_function_summary.tsv, *.pdf")
       
        pipeline_end = datetime.now()
        total_time = (pipeline_end - pipeline_start).total_seconds()
        
        logger.info(f"\n{'='*60}")
        logger.info("[Summary of Results]")
        logger.info(f"{'='*60}")
        logger.info("[OK] annotation pipeline completed!")
        logger.info(f"\nTotal Time: {total_time:.2f} s ({total_time/60:.2f} min)")
        logger.info("="*60)
        
    except KeyboardInterrupt:
        logger.error("\n\n[ABORTED] Process interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"\n{'='*60}")
        logger.error("[FAILED] annotation pipeline failed")
        logger.error(f"{'='*60}")
        logger.error(f"Error Info: {str(e)}")
        logger.error(traceback.format_exc())
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='JMPhage Annotation Module - ORF prediction, functional annotation, and genome visualization'
    )
    
    parser.add_argument('-i', '--input', required=True, help='Path to input FASTA file')
    parser.add_argument('-d', '--database', required=True, help='Path to JMP database')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    run_annotation(
        input_fasta=args.input,
        threads=args.threads,
        database=args.database,
        output_dir=args.output
    )


if __name__ == "__main__":
    main()