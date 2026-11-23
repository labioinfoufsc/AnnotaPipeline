#!/usr/bin/python3
"""Legacy CLI for the shell-driven AnnotaPipeline orchestrator."""


####################################################
###                PRIMARY SCRIPT                ###
###   IMPORT PIPELINE SCRIPTS AS FUNCTIONS AND   ###
###   PARSE THROUGH AnnotaPipeline.config FILE   ###
####################################################

# USAGE: python3 AnnotaPipelineMain.py --seq protein.fasta

# --- IMPORT PACKAGES ----------------------------------------------------------

import contextlib
import logging
import os
import pathlib
import re
import shutil
import subprocess
import sys
from shutil import which

import typer
import yaml
from annotation_io import (
    annotate_codingseq,
    fasta_fetcher,
    quantitative_proteomics,
    sequence_cleaner,
)

# --- SCRIPT LOCATION ----------------------------------------------------------

script_pwd = pathlib.Path(sys.argv[0]).absolute()

# AnnotaPipeline location
pipeline_pwd = script_pwd.parents[0]

# initial directory
home_dir_pwd = pathlib.Path.cwd()

# --- PARSER ARGUMENTS ---------------------------------------------------------

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    seq: str = typer.Option(None, "-s", "--seq"),
    protein: str = typer.Option(None, "-p", "--prot"),
    config: str = typer.Option(..., "-c", "--config"),
    gff: str = typer.Option(None, "-gff"),
) -> None:
    """Run the legacy pipeline from a YAML configuration file."""
    if (seq is None) == (protein is None):
        raise typer.BadParameter(
            "exactly one of -s/--seq or -p/--prot is required"
        )

    # Mutable state shared with the nested helper functions.
    aug_parsing = None
    kallisto_method = None
    kallisto_paired_end = None
    use_last_and_first = None

    # ------------------ CREATE LogFile ---------------------------------------

    logger = logging.getLogger("AnnotaPipeline")

    # ----------------------Redirect STDOUT and STDERR to logfile--------------

    class StreamToLogger:
        """
        Fake file-like stream object that redirects writes to a logger instance.
        """

        def __init__(self, logger, log_level=logging.INFO):
            self.logger = logger
            self.log_level = log_level
            self.linebuf = ""

        def write(self, buf):
            temp_linebuf = self.linebuf + buf
            self.linebuf = ""
            for line in temp_linebuf.splitlines(True):
                if line[-1] == "\n":
                    self.logger.log(self.log_level, line.rstrip())
                else:
                    self.linebuf += line

        def flush(self):
            if self.linebuf != "":
                self.logger.log(self.log_level, self.linebuf.rstrip())
            self.linebuf = ""

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
        datefmt="%d/%m/%Y %H:%M:%S",
        filename="AnnotaPipeline_LogFile.log",
        filemode="a",
    )

    stderr_logger = logging.getLogger("AnnotaPipeline")
    sl = StreamToLogger(stderr_logger, logging.ERROR)
    sys.stderr = sl

    # ---------------------- FUNCTIONS ---------------------------------------------
    # Function to close log and quit AnnotaPipeline if some expected file/parameter cant be found
    def log_quit():
        logger.info("Exiting")
        logging.shutdown()
        sys.exit(1)

    # check if softwares are in $PATH
    def is_tool(name):
        if which(name) is None:
            logging.error(f"Program: {name!s} must be avaliable in $PATH")
            log_quit()

    # Function to check if file was generated. It helps if AnnotaPipeline crash
    def check_file(file):
        logger = logging.getLogger("AnnotaPipeline")
        if os.path.isfile(file) == 0:
            logger.error(
                f"File {file!s} does not exist, please check earlier steps"
            )
            log_quit()
        elif os.path.getsize(file) == 0:
            logger.warning(
                f"File {file!s} is empty, please check earlier steps"
            )
            logger.warning("Trying to continue, but failures may occur")
            logger.info(
                "Sometimes things don't work as we expect... and that's ok. "
                "But, seriously, check your inputs and outputs and rerun."
            )

    # ------------------------------------- RUNs -------------------------------------------------
    # Create command line to run augustus (protein prediction)
    def augustus_run(augustus_section, augustus_optional, basename):
        # create AUGUSTUS directory
        logger = logging.getLogger("AUGUSTUS")
        if str(augustus_section.get("augustus-path")).lower() == "conda":
            # Runing with conda
            augustus_bin = "augustus"
            augustus_config = f"{os.environ['CONDA_PREFIX']}/config"
            augustus_script = "getAnnoFasta.pl"
        else:
            # Runing custom
            augustus_bin = (
                pathlib.Path(augustus_section.get("augustus-path"))
                / "bin"
                / "augustus"
            )
            augustus_config = (
                pathlib.Path(augustus_section.get("augustus-path")) / "config"
            )
            augustus_script = (
                pathlib.Path(augustus_section.get("augustus-path"))
                / "scripts"
                / "getAnnoFasta.pl"
            )

        logger.info("AUGUSTUS prediction has started")

        # AUGUSTUS: command line
        aug_config = f"--AUGUSTUS_CONFIG_PATH={augustus_config!s}"
        aug_command = f"{augustus_bin!s} {aug_config!s}"

        # Get mandatory arguments
        for param in augustus_section:
            if param != "augustus-path":
                # Check if it's a flag
                if str(augustus_section.get(param)).lower() == "flag":
                    aug_command += f" -{param!s}"
                else:
                    aug_command += (
                        f" --{param!s}={augustus_section.get(param)!s}"
                    )

        # Get optional arguments
        for param in augustus_optional:
            if str(augustus_optional.get(param)).lower() == "flag":
                aug_command += f" -{param!s}"
            else:
                # Prevent getting empty arguments
                if augustus_optional.get(param) is not None:
                    aug_command += (
                        f" --{param!s}={augustus_optional.get(param)!s}"
                    )

        aug_command += f" {seq_file!s} > AUGUSTUS_{basename!s}.gff"

        logger.debug(str(aug_command))

        subprocess.getoutput(aug_command)

        # Check if expected file exists
        check_file(f"AUGUSTUS_{basename!s}.gff")

        logger.info("AUGUSTUS prediction is finished")

        # Parsing AUGUSTUS files
        logger.info("AUGUSTUS parsing has started")

        aug_file = f"AUGUSTUS_{basename!s}"

        subprocess.run(
            [
                str(augustus_script),
                str(aug_file + ".gff"),
                str("--seqfile=" + str(seq_file)),
            ]
        )

        logger.info("AUGUSTUS parsing is finished")

        # Declare as global, so you can modify out of fuction
        nonlocal aug_parsing
        aug_parsing = aug_file + ".aa"

    # Create command line to run interproscan (functional prediction)
    def interpro_run(
        type,
        basename,
        blast_path,
        augustus_path,
        augustus_file,
        interpro_section,
        interpro_exe,
    ):
        logger = logging.getLogger("INTERPROSCAN")
        logger.info(
            f"Preparing file for INTERPROSCAN {type.capitalize()} Proteins execution"
        )

        if type == "hypothetical":
            hypothetical_id = str(
                blast_path / str(f"{basename}_{type}_products.txt")
            )
            no_hit_id = str(blast_path / str(f"{basename}_no_hit_products.txt"))
            with open(hypothetical_id) as source:
                hypothetical_id_strip = [line.strip() for line in source]
            with open(no_hit_id) as source:
                no_hit_id_strip = [line.strip() for line in source]
            concatenate_list = [*hypothetical_id_strip, *no_hit_id_strip]

            fasta_fetcher(
                str(augustus_path / str(f"Clear_{augustus_file}")),
                concatenate_list,
                f"{type.capitalize()}_Products.fasta",
            )
        else:
            annotated_file = str(
                blast_path / str(f"{basename}_{type}_products.txt")
            )
            with open(annotated_file) as source:
                annotated_id = [line.strip().split()[0] for line in source]
            fasta_fetcher(
                str(augustus_path / str(f"Clear_{augustus_file}")),
                annotated_id,
                f"{type.capitalize()}_Products.fasta",
            )

        # Check if expected file exists
        check_file(f"{type.capitalize()}_Products.fasta")

        logger.info(f"{type.capitalize()} Proteins file preparation complete")

        logger.info(f"Running with {type.capitalize()} Proteins")
        # INTERPROSCAN: commandline

        # General
        interpro_command_line = (
            f"{interpro_exe} -i {type.capitalize()}_Products.fasta "
            f"-o {basename}_interproscan_{type}_output.gff3 "
            f"-f GFF3 -t p -goterms -iprlookup"
        )

        # Optionals
        for variable in interpro_section:
            if str(interpro_section.get(variable)).lower() == "flag":
                interpro_command_line += f" -{variable!s}"
            elif (
                variable != "interproscan-exe"
                and interpro_section.get(variable) is not None
            ):
                interpro_command_line += (
                    f" -{variable!s} {interpro_section.get(variable)!s}"
                )

        logger.debug(str(interpro_command_line))

        subprocess.getoutput(interpro_command_line)

        # INTERPROSCAN parser (info_parser.py) can run without this result, but must be a valid file.
        if os.path.isfile(f"{basename}_interproscan_{type}_output.gff3") == 0:
            # Generate valid file
            open(f"{basename}_interproscan_{type}_output.gff3", "w").close()
            logger.warning(
                "INTERPROSCAN analysis return no results, moving on without this results."
            )
            logger.warning(
                "Check if your sequences have special characters (like *), remove it and rerun"
            )

        logger.info(f"INTERPROSCAN finished for {type.capitalize()} Proteins")

    # Create command line to run kallisto (transcriptomics)
    def kallisto_run(
        python_path, kallisto_exe, paired_end, method, basename, fasta
    ):

        logger.info("KALLISTO index has started")
        kallisto_command_index = (
            f"{kallisto_exe} index -i {basename}_kallisto_index.idx {fasta}"
        )
        logger.debug(f"{kallisto_command_index}")
        subprocess.getoutput(kallisto_command_index)
        check_file(f"{basename}_kallisto_index.idx")

        # Standart command line for kallisto index
        # kallisto index -i transcripts.idx transcripts.fasta
        if paired_end:
            if kallisto.get("l") is not None:
                l_flag = f"-l {kallisto.get('l')}"
                logger.info(f"KALLISTO will run with -l {kallisto.get('l')}")
            else:
                l_flag = ""
            if kallisto.get("s") is not None:
                s_flag = f"-s {kallisto.get('s')}"
                logger.info(f"KALLISTO will run with -s {kallisto.get('s')}")
            else:
                s_flag = ""

            kallisto_command_quant = (
                f"{kallisto_exe} quant -i {basename}_kallisto_index.idx "
                f"{s_flag} {l_flag} "
                f"-o {basename}_kallisto_output -b {kallisto.get('bootstrap')} {' '.join(kallisto.get('rna-seq'))}"
            )
            logger.info("KALLISTO quant has started")
            logger.debug(kallisto_command_quant)
            subprocess.getoutput(kallisto_command_quant)
            # Standart command line for kallisto quant paired end
            # kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq reads_2.fastq
        else:
            kallisto_command_quant = (
                f"{kallisto_exe} quant -i {basename}_kallisto_index.idx "
                f"-l {kallisto.get('s')} -s {kallisto.get('l')} --single "
                f"-o kallisto_output -b {kallisto.get('bootstrap')} {' '.join(kallisto.get('rna-seq'))}"
            )
            logger.info("KALLISTO quant has started")
            logger.debug(kallisto_command_quant)
            subprocess.getoutput(kallisto_command_quant)
            # Standart command line for kallisto quant single end
            # kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 reads_1.fastq
        check_file(f"{basename}_kallisto_output/abundance.tsv")
        # Define method to parse
        if method == "median":
            kallisto_parser_flag = "-tpmmd"
        elif method == "mean":
            kallisto_parser_flag = "-tpmavg"
        else:
            kallisto_parser_flag = f"-tpmval {kallisto.get('threshold')}"
        # Run parser
        # inside 4_Transcript_Quantification_
        kallisto_parser_path = str(pipeline_pwd / "kallisto_parser.py")
        kallisto_parser_command = (
            f"{python_path} {kallisto_parser_path} "
            f"-ktfile {basename}_kallisto_output/abundance.tsv "  # kallisto output is always "abundance.tsv"
            f"-basename {basename} {kallisto_parser_flag}"
        )
        logger.info("KALLISTO parsing has started")
        subprocess.getoutput(kallisto_parser_command)
        check_file(f"{basename}_Transcript_Quantification.tsv")

    # Run fastatogff parser
    def run_fastatogff(python_path):
        subprocess.run(
            [
                python_path,
                str(pipeline_pwd / "fastatogff.py"),
                "-gff",
                str(gff_file),
                "-all",
                "All_Annotated_Products.txt",
                "-b",
                str(basename),
            ]
        )

    # Run gff to fasta parser
    def gfftofasta(python_path):
        subprocess.run(
            [
                python_path,
                str(pipeline_pwd / "gfftofasta_parser.py"),
                "-gff",
                str(gff_file),
                "-annot",
                "All_Annotated_Products.txt",
                "-b",
                str(basename),
                "-faf",
                str(augustus_folder / str("Clear_" + aug_parsing)),
                "-org",
                str('"{}"'.format(str(pipeline_settings.get("organism")))),
            ]
        )

    # ------------------------------------------ -------------------------------------------------
    # ------------------------- Check Parameters -------------------------------------------------
    # Part of check parameters function >> check specific entries for kallisto (transcriptomics)
    def kallisto_check_parameters(kallisto_section):
        # kallisto method indicate wich method to parse (mean, median, value)
        # kallisto paired_end indicate rnaseq type (important to run kallisto)
        nonlocal kallisto_method, kallisto_paired_end
        # This variable will store method to parse kallisto output
        kallisto_method = None
        # if kallisto path were given, check other arguments, if not pass.
        if kallisto_section.get("kallisto-exe") is None:
            logger.info(
                "Arguments for KALLISTO are empty. This step will be skipped."
            )
        else:
            kallisto_check = []
            for argument in ("rna-seq", "threshold"):
                if argument == "rna-seq":
                    if len(kallisto_section.get("rna-seq")) > 2:
                        logger.error(
                            "[KALLISTO]: there are more arguments than required for rnaseq-data (KALLISTO). "
                            "Pass one if your data is from single-end, "
                            "and two files if your data is from paired-end!"
                        )
                        log_quit()
                    ###### This box check how many files were given in rnaseq-data #####
                    elif len(kallisto_section.get("rna-seq")) == 2:
                        kallisto_paired_end = True
                        logger.info("KALLISTO will run with paired end data")
                        kallisto_check.append(argument)
                    elif len(kallisto_section.get("rna-seq")) == 1:
                        kallisto_paired_end = False
                        # Check required arguments for single end data
                        if kallisto_section.get("l") is None:
                            logger.error(
                                "[KALLISTO]: mandatory argument for single-end data 'l' is empty"
                            )
                            log_quit()
                        elif kallisto_section.get("s") is None:
                            logger.error(
                                "[KALLISTO]: mandatory argument for single-end data 's' is empty"
                            )
                            log_quit()
                        logger.info("KALLISTO will run with single-end data")
                        kallisto_check.append(argument)
                    else:
                        logger.error(
                            "[KALLISTO]: check values for rnaseq-data (KALLISTO). "
                            "Pass one if your data is from single-end, "
                            "and two files if your data is from paired-end!"
                        )
                        log_quit()
                    ####################################################################
                # Check if argument is Empty
                elif argument == "threshold":
                    if kallisto_section.get("threshold") is None:
                        logger.error(
                            "[KALLISTO]: mandatory argument 'threshold' is empty"
                        )
                        log_quit()
                    else:
                        if kallisto_section.get("threshold") == "mean":
                            kallisto_method = "mean"
                        elif kallisto_section.get("threshold") == "median":
                            kallisto_method = "median"
                        elif type(kallisto_section.get("threshold")) in (
                            int,
                            float,
                        ):
                            kallisto_method = "value"
                        else:
                            logger.error(
                                "[KALLISTO]: Unrecognized type for argument 'threshold'"
                            )
                            log_quit()
            # Check kallisto arguments
            # If all kallisto arguments are empty, then don't run this guy
            if len(kallisto_check) == 0:
                logger.info(
                    "Arguments for KALLISTO are empty. This step will be skipped."
                )
                kallisto_method = None
            # Rna-seq data is required for kallisto
            if "rna-seq" in kallisto_check:
                # if there is rna-seq data, check if method is correctly given
                if kallisto_section.get("bootstrap") is None:
                    logger.error(
                        "[KALLISTO]: bootstrap is empty, default value is 0. At least pass this value"
                    )
                    log_quit()
                else:
                    logger.info(
                        f"KALLISTO will run with method: {kallisto_method}"
                    )

    # Part of check parameters function >> check specific entries for comet (proteomics)
    def proteomics_check_parameters(proteomic_section):
        if proteomic_section.get("comet-exe") is None:
            logger.info(
                "Arguments for COMET are empty. This step will be skipped."
            )
        else:
            # ------------ check F and L -------------------
            last_check = False
            first_check = False
            if proteomic_section.get("first") is not None:
                first_check = True
            if proteomic_section.get("last") is not None:
                last_check = True
            # ----------------------------------------------
            for argument in ("comet-params", "comet-spectrometry"):
                if proteomic_section.get(argument) == 0:
                    logger.error(
                        f"[COMET] Parameter [{argument}] from section [COMET] is null"
                    )
                    log_quit()
            # ------------ check extension mass files ------
            if len(proteomic_section.get("comet-ext").split()) < 1:
                logger.error(
                    f"[COMET] Parameter [{proteomic_section.get('mass_files_ext')}] from section [PROTEOMICS] is null"
                )
                log_quit()
            elif len(proteomic_section.get("comet-ext").split()) > 1:
                logger.error(
                    f"[COMET] Parameter [{proteomic_section.get('mass_files_ext')}] from section [PROTEOMICS] have more than one argument"
                )
                log_quit()
            # ------------ Specific conditions -----------------------
            # XOR conditional >> check if one value is True and another is False
            if last_check ^ first_check:
                logger.error(
                    "[COMET]: both arguments from comet, first and last, must be given"
                )
                logger.error("[COMET]: Leave both empty or give both")
                log_quit()
            elif last_check == first_check is True:
                logger.warning(
                    "COMET values for 'fist' and 'last' will overwrite those in params file"
                )
                nonlocal use_last_and_first
                use_last_and_first = True
            else:
                use_last_and_first = False
            # -------------------------------------------------
            # ------ Check parameters for percolator ----------
            # Only check percolator if comet parames are ok
            if proteomic_section.get("percolator-exe") is None:
                logger.error(
                    "[PERCOLATOR]: path to software is empty, but COMET was provided"
                )
                logger.warning(
                    "PERCOLATOR leave COMET fields empty or pass parameters for PERCOLATOR"
                )
                log_quit()
            if proteomic_section.get("percolator-qvalue") is None:
                logger.error(
                    "[PERCOLATOR]: qvalue cutoff is empty, check this parameter"
                )
                log_quit()
            # check value for parser
            if not (
                0 <= float(proteomic_section.get("percolator-qvalue")) <= 1
            ):
                logger.error(
                    "[PERCOLATOR] qvalue cutoff invalid. Must be float between [0-1]"
                )
                log_quit()

    # Function to check if all parameters in config file are correct
    def check_parameters(sections):
        for section, list_section in sections.items():
            if str(section) == "kallisto":
                kallisto_check_parameters(list_section)
            elif str(section) == "proteomics":
                proteomics_check_parameters(list_section)
            elif str(section) == "databases":
                params_empty = [
                    param
                    for param in list_section
                    if list_section.get(param) is None
                ]
                # When databse is not custom, parameters customsep and customcolumn aren't used
                if list_section.get("secondary-format") != "custom":
                    with contextlib.suppress(ValueError):
                        params_empty.remove("customsep")
                    with contextlib.suppress(ValueError):
                        params_empty.remove("customcolumn")
                if params_empty:
                    logger.error(
                        f"[DATABASE]: there is some secondary database parameters missing. Params empty: {' '.join(params_empty)}. Please review config file!"
                    )
                    log_quit()
                else:
                    # Check if format of secondary database is correct
                    if list_section["secondary-format"] not in (
                        "eupathdb",
                        "nrdb",
                        "trembldb",
                        "custom",
                    ):
                        logger.error(
                            "[DATABASE]: Unrecognized type for 'secondary-format'. Please review config file!"
                        )
                        log_quit()
            elif str(section) == "augustus":
                # Arguments for augustus don't need to be checked if protein file were given
                if protein is not None:
                    pass
                else:
                    for key in list_section:  # get variable for each box
                        # MS is optional parameter for INTERPROSCAN
                        if str(key) == "ms":
                            pass
                        elif list_section.get(key) is None:
                            # Crash pipeline if some required variable is empty
                            logger.error(
                                f"Variable [{key}] from section [{section!s}] is null"
                            )
                            log_quit()
            # augustus optionals are optionals
            elif str(section) in "augustus-optional":
                pass
            else:
                # check if variable in list each section
                for key in list_section:  # get variable for each box
                    # MS is optional parameter for INTERPROSCAN
                    if str(key) == "ms":
                        pass
                    elif list_section.get(key) is None:
                        # Crash pipeline if some required variable is empty
                        logger.error(
                            f"Variable [{key}] from section [{section!s}] is null"
                        )
                        log_quit()

    # -------------------------------------------------------------------------------------------

    # Function to change essential parameters to run percolator and parse output files
    def modify_comet_params(comet_params):
        with open(f"{comet_params}") as comet_params_read:
            new_comet_params_read = comet_params_read.read()
            new_comet_params_read = re.sub(
                r"decoy_search\s*=\s*[0-9]",
                "decoy_search = 1",
                new_comet_params_read,
            )
            new_comet_params_read = re.sub(
                r"output_pepxmlfile\s*=\s*[0-9]",
                "output_pepxmlfile = 0",
                new_comet_params_read,
            )
            new_comet_params_read = re.sub(
                r"output_percolatorfile\s*=\s*[0-9]",
                "output_percolatorfile = 1",
                new_comet_params_read,
            )
            # \s escape
            # this is not a bug, its a feature
            whitespace = (
                " " * 18
            )  # as of python3.8.5 this is the only way to write \s using re package
            try:
                new_comet_params_read = re.sub(
                    r"decoy_prefix\s*=\s*.*#",
                    rf"decoy_prefix = DECOY_{whitespace}#",
                    new_comet_params_read,
                )
            except Exception as warn:
                logger.warning(
                    f"Cannot replace decoy_prefix in {comet_params}, change this value manually to: DECOY_"
                )
                logger.debug(warn)
        with open(f"{comet_params}", "w") as comet_params_write:
            comet_params_write.write(new_comet_params_read)

    # --- PARSE ARGUMENTS FROM .config FILE ----------------------------------------

    config_pwd = pathlib.Path(config).absolute()

    # Get params from yaml
    with open(config_pwd) as stream:
        try:
            config = yaml.load(stream, Loader=yaml.SafeLoader)
        except yaml.YAMLError as exc:
            logger.error(exc)
            log_quit()

    # ------------------------------------------------------------------------------
    # --- CHECK EACH BOX OF VARIABLES ----------------------------------------------
    # is_tool("perl")
    # ------------------------------------------------------------------------------

    # logger.info("IT'S DANGEROUS TO GO ALONE! TAKE THIS.")

    logger = logging.getLogger("AnnotaPipeline")
    logger.info(
        "---------------------------------------------------------------"
    )
    logger.info("Checking parameters in annotapipeline.yaml file")
    check_parameters(config)

    # --- PREPARING SOME VARIABLES -------------------------------------------------
    # \\ Check if user pass protein and gff file -> if it is, redirect variables
    if seq is not None:
        seq_file = pathlib.Path(seq).absolute()
    if protein is not None:
        prot_path = pathlib.Path(protein).absolute()
    if gff is not None:
        gff_path = pathlib.Path(gff).absolute()

    # Parameters from config file, each line is one script/software configuration
    python_exe = config["pipeline"]["python"]
    pipeline_settings = config["pipeline"]
    databases = config["databases"]
    basename = pipeline_settings["basename"]
    keyword_list = pipeline_settings["keywords"]  # join words with '
    augustus_main = config["augustus"]
    augustus_optional = config["augustus-optional"]
    seq_cleaner = config["seq-cleaner"]
    interpro = config["interproscan"]
    hmmscan = config["hmmer"]
    blast = config["local-aligner"]
    rpsblast = config["rpsblast"]
    kallisto = config["kallisto"]
    proteomics = config["proteomics"]

    # AnnotaPipeline main directory
    home_dir = f"AnnotaPipeline_{basename!s}"

    # FileExistsError exception will be ignored
    pathlib.Path(home_dir).mkdir(exist_ok=True)
    os.chdir(home_dir)
    annota_pwd = pathlib.Path(home_dir_pwd / home_dir)

    # --- AUGUSTUS -----------------------------------------------------------------

    augustus_folder = pathlib.Path(
        annota_pwd / str("1_GenePrediction_" + basename)
    )
    pathlib.Path(augustus_folder).mkdir(exist_ok=True)
    os.chdir(augustus_folder)

    # ==============================================================================
    # Run Augustus or start with protein file?
    logger.info(
        "---------------------------------------------------------------"
    )
    logger.info(
        "----------------- Gene Prediction has started -----------------"
    )
    if protein is None:
        augustus_run(augustus_main, augustus_optional, basename)
    else:
        # Copy protein file to AUGUSTUS path and patronized variable to run Annotapipeline after augustus
        shutil.copy2(prot_path, augustus_folder)
        aug_parsing = protein

    # ==============================================================================
    # SEQUENCE CLEANER -------------------------------------------------------------
    logger = logging.getLogger("AnnotaPipeline")
    logger.info("Sequence Cleaner has started")

    # Clean only with min_size
    sequence_cleaner(str(aug_parsing), int(seq_cleaner.get("minsize-seq")))

    # Check if expected file exists
    check_file(f"Clear_{aug_parsing}")

    logger.info(
        f"Sequence Cleaner is finished. Please check Clear_{aug_parsing}"
    )

    os.chdir(annota_pwd)

    # BLAST ------------------------------------------------------------------------

    # Create BLAST directory
    blast_folder = pathlib.Path(
        annota_pwd / str("2_SimilarityAnalysis_" + basename)
    )
    pathlib.Path(blast_folder).mkdir(exist_ok=True)
    os.chdir(blast_folder)

    # BLAST: command line
    logger.info(
        "---------------------------------------------------------------"
    )
    logger.info(
        "--------------- Similarity Analysis has started ---------------"
    )
    logger.info("BLAST execution and parsing has started")

    # Select secondary database from config file
    spdb_path = databases.get("secondary-db")
    # Default custom sep, wich is used only form customdb, so its defined here
    # and will be ignored in blastp_parser.py for other databases
    customsep = "|"
    customfield = 5
    if str(databases.get("secondary-format")).lower() == "eupathdb":
        flag_spdb = "-spdb"
    elif str(databases.get("secondary-format")).lower() == "trembldb":
        flag_spdb = "-trbl"
    elif str(databases.get("secondary-format")).lower() == "nrdb":
        flag_spdb = "-nr"
    elif str(databases.get("secondary-format")).lower() == "custom":
        flag_spdb = "-spdb"
        customsep = databases.get("customsep")
        customfield = databases.get("customcolumn")

    subprocess.run(
        [
            str(python_exe),
            str(str(pipeline_pwd / "blastp_parser.py")),
            "-s",
            str(augustus_folder / str("Clear_" + aug_parsing)),
            "-sp",
            str(databases.get("swissprot-db")),
            "-basename",
            str(basename),
            str(flag_spdb),  # Flag for databse
            str(spdb_path),  # path to database
            "-id",
            str(blast.get("identity")),
            "-pos",
            str(blast.get("positivity")),
            "-cov",
            str(blast.get("coverage")),
            "-kw",
            str(f"{','.join(keyword_list)}"),
            "-t",
            str(pipeline_settings.get("threads")),
            "-max_target_seqs",
            str(blast.get("max_target_seqs")),
            "-evalue",
            str(blast.get("evalue")),
            "-blastp",
            str(blast.get("blast-exe")),
            # Flags used only with customdb
            "-customsep",
            str(customsep),
            "-customcolumn",
            str(customfield),
        ]
    )

    logger.info("BLAST execution and parsing is finished")

    # Check if expected file exists
    check_file(f"{basename!s}_no_hit_products.txt")
    check_file(f"{basename!s}_hypothetical_products.txt")
    check_file(f"{basename!s}_annotated_products.txt")

    os.chdir(annota_pwd)

    # INTERPROSCAN -----------------------------------------------------------------

    # Create INTERPROSCAN directory
    interpro_folder = pathlib.Path(
        annota_pwd / str("3_FunctionalAnnotation_" + basename)
    )
    pathlib.Path(interpro_folder).mkdir(exist_ok=True)
    os.chdir(interpro_folder)
    logger.info(
        "---------------------------------------------------------------"
    )
    logger.info(
        "------------- Functional Annotation has started ---------------"
    )
    # Running interproscan with hypothetical proteins
    interpro_run(
        "hypothetical",
        basename,
        blast_folder,
        augustus_folder,
        aug_parsing,
        interpro,
        interpro.get("interproscan-exe"),
    )

    # Running interproscan with annotated proteins
    interpro_run(
        "annotated",
        basename,
        blast_folder,
        augustus_folder,
        aug_parsing,
        interpro,
        interpro.get("interproscan-exe"),
    )

    # HMMER ------------------------------------------------------------------------
    logger = logging.getLogger("HMMSCAN")
    logger.info("Running HMMSCAN with Hypothetical Proteins")

    # General
    hmmscan_command_line = (
        f"{hmmscan.get('hmmscan-exe')!s} "
        f"--cpu {pipeline_settings.get('threads')!s} "
        f"--tblout {basename!s}_hmmscan_output.txt "
        f"--noali"
    )

    # Optionals
    for param in hmmscan:
        if str(hmmscan.get(param)).lower() == "flag":
            hmmscan_command_line += f" --{param!s}"
        elif param != "hmmscan-exe":
            # These specific arguments are passed through '-'
            if any(arg == str(param) for arg in ("E", "Z", "T")):
                hmmscan_command_line += f" -{param!s} {hmmscan.get(param)!s}"
            # Everyone else are passed through '--'
            else:
                hmmscan_command_line += f" --{param!s} {hmmscan.get(param)!s}"

    # Add Database
    hmmscan_command_line += (
        f" {databases.get('pfam')!s} "
        f"Hypothetical_Products.fasta "
        f"> /dev/null 2> hmmscan.err"
    )

    logger.debug(str(hmmscan_command_line))

    subprocess.getoutput(hmmscan_command_line)

    # Check if expected file exists
    check_file(f"{basename!s}_hmmscan_output.txt")

    logger.info("HMMSCAN is finished")

    # RPSBLAST ---------------------------------------------------------------------
    logger = logging.getLogger("RPSBLAST")
    logger.info("Running RPSBLAST with Hypothetical Proteins")

    # General
    rpsblast_command_line = (
        f"{rpsblast.get('rpsblast-exe')!s} -query Hypothetical_Products.fasta "
        f"-out {basename!s}_rpsblast_output.outfmt6 "
        f"-db {databases.get('cdd-db')!s} "
        f'-outfmt "6 qseqid sseqid sacc bitscore evalue ppos pident qcovs stitle" '
        f"-num_threads {pipeline_settings.get('threads')!s}"
    )

    # Optionals
    for param in rpsblast:
        if str(rpsblast.get(param)).lower() == "flag":
            rpsblast_command_line += f" -{param!s}"
        elif param != "rpsblast-exe":
            rpsblast_command_line += f" -{param!s} {rpsblast.get(param)!s}"

    logger.debug(str(rpsblast_command_line))

    subprocess.getoutput(rpsblast_command_line)

    # Check if expected file exists
    check_file(f"{basename!s}_rpsblast_output.outfmt6")

    logger.info("RPSBLAST is finished")

    # -----------------------------------------------------------------------
    logger = logging.getLogger("AnnotaPipeline")
    logger.info(
        "---------------------------------------------------------------"
    )
    logger.info(
        "---------------- Generating Annotation Files ------------------"
    )
    logger.info("Parsing information from INTERPROSCAN, HMMSCAN and RPSBLAST")

    subprocess.run(
        [
            str(python_exe),
            str(pipeline_pwd / "funcannotation_parser.py"),
            "-ipr_annot",
            str(basename + "_interproscan_annotated_output.gff3"),
            "-ipr_hyp",
            str(basename + "_interproscan_hypothetical_output.gff3"),
            "-hmm",
            str(basename + "_hmmscan_output.txt"),
            "-rpsblast",
            str(basename + "_rpsblast_output.outfmt6"),
            "-basename",
            str(basename),
        ]
    )

    # Cleaning the house
    try:
        os.remove("hmmscan.err")
        shutil.rmtree("temp", ignore_errors=True)
    except Exception as warn:
        logger.warning("Failed to remove HMMSCAN log and temp dir")
        logger.debug(f"code error: {warn}")

    logger.info(
        "INTERPROSCAN, HMMSCAN and RPSBLAST execution and parsing is finished"
    )

    os.chdir(annota_pwd)

    # -----------------------------------------------------------------------

    logger.info("Preparing file for protein annotation")

    subprocess.run(
        [
            str(python_exe),
            str(pipeline_pwd / "info_parser.py"),
            "-ipr1",
            str(
                interpro_folder
                / str(basename + "_interproscan_annotated_output.gff3")
            ),
            "-ipr2",
            str(
                interpro_folder
                / str(basename + "_interproscan_hypothetical_output.gff3")
            ),
            "-a",
            str(blast_folder / str(basename + "_annotated_products.txt")),
            "-hy",
            str(blast_folder / str(basename + "_hypothetical_products.txt")),
            "-nh",
            str(blast_folder / str(basename + "_no_hit_products.txt")),
        ]
    )

    # Check if expected file exists
    check_file("All_annotation_products.txt")

    try:
        # Sort annotations
        os.system(
            "sort -V All_annotation_products.txt -o All_Annotated_Products.txt"
        )
        os.remove("All_annotation_products.txt")
    except Exception as warn:
        logger.warning("Failed to sort All_annotation_products.txt")
        logger.debug(f"code error: {warn}")

    logger.info("All_Annotated_Products.txt file is complete")

    # ------------  Defining what file will be used ---------------------------
    if (
        gff is not None and protein is not None
    ):  # User gave protein file and gff file
        # Run parser to generate fasta_file
        gff_file = str(gff_path)
        logger.info(f"Generating AnnotaPipeline_{basename}_proteins.fasta")
        gfftofasta(str(python_exe))
        logger.info("Generating annotated GFF file")
        run_fastatogff(str(python_exe))
        logger.info(f"GFF file is ready - Check {basename}_Annotated_GFF.gff")
    elif protein is not None and gff is None:  # User gave only protein file
        logger.info(
            "GFF file wasn't given, FASTA file will only have annotations"
        )
        logger.info("Generating simple FASTA")
        subprocess.run(
            [
                str(python_exe),
                str(pipeline_pwd / "fasta_simple.py"),
                "-annot",
                "All_Annotated_Products.txt",
                "-b",
                str(basename),
                "-faf",
                str(augustus_folder / str("Clear_" + aug_parsing)),
                "-org",
                str('"{}"'.format(str(pipeline_settings.get("organism")))),
            ]
        )
        logger.info("GFF file wasn't given, skipping script fastatogff.py")
    else:  # User selected run Augustus
        gff_file = augustus_folder / str("AUGUSTUS_" + str(basename) + ".gff")
        logger.info(f"Generating AnnotaPipeline_{basename}_proteins.fasta")
        gfftofasta(str(python_exe))
        logger.info("Generating annotated GFF file")
        run_fastatogff(str(python_exe))
        logger.info(f"GFF file is ready - Check {basename!s}_Annotated_GFF.gff")
        # Transfer annotation from proteins to Transcript (codingseq) file
        annotate_codingseq(
            annota_pwd / str("AnnotaPipeline_" + basename + "_proteins.fasta"),
            augustus_folder / str("AUGUSTUS_" + basename + ".codingseq"),
            basename,
        )

    logger.info(
        "AnnotaPipeline has annotated the annotations on the annotated file."
    )

    # -----------------------------------------------------------------------
    # ---------------------- Kallisto ---------------------------------------
    # Run kallisto
    # If kalisto_method is empty, lack arguments for kallisto >> skip
    # If proteins were given, lack cdscexon files >> skip
    if kallisto_method is not None and seq is not None:
        logger.info(
            "---------------------------------------------------------------"
        )
        logger.info(
            "---------- Transcript Quantification has started --------------"
        )
        kallisto_output_path = pathlib.Path(
            annota_pwd / str("4_TranscriptQuantification_" + basename)
        )
        pathlib.Path(kallisto_output_path).mkdir(exist_ok=True)
        # Go to /4_TranscriptQuantification_
        os.chdir(kallisto_output_path)

        # Change logger
        logger = logging.getLogger("KALLISTO")

        # Annotated_Products.cdsexon
        kallisto_run(
            str(python_exe),
            kallisto.get("kallisto-exe"),
            kallisto_paired_end,
            kallisto_method,
            basename,
            f"{annota_pwd / f'AnnotaPipeline_{basename}_transcripts.fasta'}",
        )

        logger.info("KALLISTO execution and parsing is finished")

        # Return to AnnotaPipeline basedir
        os.chdir(annota_pwd)

    # ------------------------------------------------------------------------
    # ----------------------- Commet ----------------------------------------
    if proteomics.get("comet-exe") is not None:
        logger.info(
            "---------------------------------------------------------------"
        )
        logger.info(
            "------------ Peptide Identification has started ---------------"
        )
        if kallisto_method is None or seq is None:
            comet_output_path = pathlib.Path(
                annota_pwd / str("4_PeptideIdentification_" + basename)
            )
        else:
            comet_output_path = pathlib.Path(
                annota_pwd / str("5_PeptideIdentification_" + basename)
            )

        logger = logging.getLogger("COMET")
        # Go to /X_PeptideIdentification
        pathlib.Path(comet_output_path).mkdir(exist_ok=True)
        os.chdir(comet_output_path)

        # Check if overwrite parameters will be used
        if use_last_and_first:
            first_last_param = (
                f"-F{proteomics.get('first')} -L{proteomics.get('last')}"
            )
        else:
            first_last_param = ""
        # Mass files location
        mass_path = f"{str(proteomics.get('comet-spectrometry')).rstrip('/')}/"

        # Add percolator output in comet.params
        # change arguments for comet
        modify_comet_params(proteomics.get("comet-params"))
        commet_command = (
            f"{proteomics.get('comet-exe')} -P{proteomics.get('comet-params')} "
            f"-D{annota_pwd / f'AnnotaPipeline_{basename}_proteins.fasta'} "
            f"{first_last_param} {mass_path}*.{proteomics.get('comet-ext')}"
        )

        logger.info("COMET execution has started")
        logger.info(commet_command)
        subprocess.getoutput(commet_command)
        logger.info("COMET execution is finished")
        # ----------- Create Comet Output path --------------------------
        comet_path = pathlib.Path(comet_output_path / "COMET_Output")
        pathlib.Path(comet_path).mkdir(exist_ok=True)
        # ----------------------------------------------------------------
        # Get all output files from mass_path >> default output path
        file_names = pathlib.Path(mass_path).glob("*.pin")
        if not list(file_names):
            logger.error(
                "COMET - Fail trying to run COMET, check version of comet.params"
            )
            log_quit()
        # ----------------------------------------------------------------
        # MUST GET AGAIN BECAUSE LIST METHOD EMPTY VALUES IF IT HAD ONE
        file_names = pathlib.Path(mass_path).glob("*.pin")
        # ----------- Create Percolator Output path ----------------------
        logger = logging.getLogger("PERCOLATOR")
        logger.info("PERCOLATOR execution has started")
        # ----------------------------------------------------------------
        # ----------- Create Percolator Output path ----------------------
        percolator_path_raw = pathlib.Path(comet_output_path / "PERCOLATOR_Raw")
        pathlib.Path(percolator_path_raw).mkdir(exist_ok=True)
        # ----------------------------------------------------------------
        percolator_path_parsed = pathlib.Path(
            comet_output_path / "PERCOLATOR_Parsed"
        )
        pathlib.Path(percolator_path_parsed).mkdir(exist_ok=True)
        # ----------------------------------------------------------------
        for comet_output_file in file_names:
            # -------- check files created by comet to run percolator ---------------
            check_file(f"{comet_output_file}")
            # --------- RUN Percolator inside Percolator RAW path -------------------
            os.chdir(percolator_path_raw)
            percolator_out_basename = re.sub(
                r"\.[0-9].*", "", comet_output_file.stem
            )
            percolator_command = (
                f"{proteomics.get('percolator-exe')} -r {percolator_out_basename}_peptide_output.tsv"
                f" -m {percolator_out_basename}_percolator_output.tsv"
                f" -B {percolator_out_basename}_decoy_output.tsv {comet_output_file}"
            )
            subprocess.getoutput(percolator_command)
            check_file(f"{percolator_out_basename}_percolator_output.tsv")
            # --------- RUN Percolator parser inside Percolator PARSED path -------------------
            os.chdir(percolator_path_parsed)
            parser_percolator_command = (
                f"{python_exe} {pipeline_pwd / 'percolator_parser.py'!s}"
                f" -p {percolator_path_raw}/{percolator_out_basename}_percolator_output.tsv"
                f" -qv {proteomics.get('percolator-qvalue')}"
                f" -b {basename}_{percolator_out_basename}"
            )
            try:
                subprocess.getoutput(parser_percolator_command)
            except Exception as warn:
                # Empty sheets will fall here
                logger.warning(
                    f"Failed trying to parser {percolator_out_basename}_percolator_output.tsv"
                )
                logger.debug(f"code error {warn}")

            # -----------------------------------------
            # Move Comet output files
            try:
                shutil.move(str(comet_output_file), str(comet_path))
            except Exception as warn:
                logger.warning(
                    f"Failed trying move comet output files to {comet_path}"
                )
                logger.debug(f"code error {warn}")

            os.chdir(comet_output_path)

        logger.info("PERCOLATOR execution and parsing is finished")
        logger.info("Creating quantitative report of Spectrum and Peptides")
        quantitative_proteomics(f"{percolator_path_parsed}", basename)
        try:
            # Sort spectrum count
            os.system(
                f"sort -V {basename}_pre_total_Proteomics_Quantification.tsv "
                f"-o  {basename}_Total_Proteomics_Quantification.tsv"
            )
            os.remove(f"{basename}_pre_total_Proteomics_Quantification.tsv")
        except Exception as warn:
            logger.warning(
                f"Failed to sort {basename}_pre_total_Proteomics_Quantification.tsv"
            )
            logger.debug(f"code error: {warn}")

    # Return to AnnotaPipeline basedir
    os.chdir(annota_pwd)
    logger = logging.getLogger("AnnotaPipeline")
    # --------------------------------------------------------------------------------------------------
    # summary PARSER
    summary_parser_command_line = (
        f"{python_exe} {pipeline_pwd / 'summary_parser.py'!s}"
        f" -b {basename}"
        f" -annot {annota_pwd / 'All_Annotated_Products.txt'!s}"
        f" -ipr_hyp {interpro_folder / str(basename + '_interproscan_hypothetical_output.gff3')!s}"
        f" -ipr_annot {interpro_folder / str(basename + '_interproscan_annotated_output.gff3')!s}"
    )

    # Add optional parametes (if kallisto and/or comet were executed)
    if kallisto_method is not None and seq is not None:
        summary_parser_command_line = (
            summary_parser_command_line
            + f" -tr {kallisto_output_path / str(basename + '_Transcript_Quantification.tsv')!s}"
        )
    if proteomics.get("comet-exe") is not None:
        summary_parser_command_line = (
            summary_parser_command_line
            + f" -proteomics {comet_output_path / str(basename + '_Total_Proteomics_Quantification.tsv')!s}"
        )

    # This is the end of the line
    try:
        logger.info("Generating Summary file")
        subprocess.getoutput(summary_parser_command_line)
    except Exception:
        logger.error("Failed to create summary file, check files")
        logger.debug(summary_parser_command_line)

    logger.info("That's all folks")
    # close logger
    logging.shutdown()


if __name__ == "__main__":
    app()
