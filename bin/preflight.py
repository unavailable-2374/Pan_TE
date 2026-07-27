#!/usr/bin/env python3
"""Environment preflight for Pan_TE.

Every dependency the requested run will actually touch is verified in one pass before
any work starts. The pipeline shells out to ~20 external tools and two Perl modules
spread across five stages, so without this a missing dependency surfaces as an opaque
crash minutes or hours in -- which is exactly how a fresh cluster install used to fail.

Checks are gated by what this particular run needs: classification tools are only
required when --model-dir is given, GFA/PAF tools only in those modes. Failures are
collected rather than raised one at a time, so a single run reports everything that is
wrong, and each one names its fix.
"""

import os
import shutil
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

# Tools every run needs, grouped by the stage that invokes them. The stage label is what
# the operator sees, so it points at the failing part of the pipeline directly. The last
# field lists the pipeline steps that need the tool: a tool is only *required* while at
# least one of its steps still has to run. This matters for resubmitting a partially
# completed job -- a checkpointed Step 2 must not be able to block a Step 4 resume just
# because look4ltrs is missing from the environment now.
CORE_TOOLS: List[Tuple[str, str, str, Tuple[str, ...]]] = [
    ('samtools',     'genome processing',   'bioconda::samtools',      ('genome', 'mdl_repeat', 'te_looker')),
    ('bedtools',     'genome processing',   'bioconda::bedtools',      ('genome', 'ltr')),
    ('look4ltrs',    'Step 2 (LTR)',        'bioconda::look4ltrs',     ('ltr',)),
    ('FastGA',       'Step 2 (LTR)',        'bioconda::fastga',        ('ltr',)),
    ('mdl-repeat',   'Step 3 (mdl-repeat)', 'bioconda::mdl-repeat',    ('mdl_repeat',)),
    ('dustmasker',   'Step 3 (mdl-repeat)', 'bioconda::blast',         ('mdl_repeat',)),
    ('trf',          'Step 3 (mdl-repeat)', 'bioconda::trf',           ('mdl_repeat',)),
    ('rmblastn',     'Step 3 (Refiner)',    'bioconda::rmblast',       ('mdl_repeat',)),
    ('blastn',       'Step 3 (Refiner)',    'bioconda::blast',         ('mdl_repeat',)),
    ('makeblastdb',  'Step 3 (Refiner)',    'bioconda::blast',         ('mdl_repeat',)),
    ('mafft',        'Step 3 (Refiner)',    'bioconda::mafft',         ('mdl_repeat',)),
    ('abpoa',        'Step 3 (Refiner)',    'bioconda::abpoa',         ('mdl_repeat',)),
    ('RepeatMasker', 'Step 3 (Refiner)',    'bioconda::repeatmasker',  ('mdl_repeat', 'anno')),
    ('dtr',          'Step 4 (te-looker)',  'bioconda::te-looker',     ('te_looker',)),
    ('hmmbuild',     'Step 4 (te-looker)',  'bioconda::hmmer',         ('te_looker',)),
    ('cd-hit-est',   'Step 5 (Combine)',    'bioconda::cd-hit',        ('te_looker', 'combine')),
]

# Used by optional/rescue paths. Their absence degrades a specific feature rather than
# stopping the run, so they warn instead of failing.
OPTIONAL_TOOLS: List[Tuple[str, str, str, Tuple[str, ...]]] = [
    ('nhmmer',   'te-looker boundary refine',  'bioconda::hmmer',    ('te_looker',)),
    ('sourmash', 'te-looker Track 2',          'bioconda::sourmash', ('te_looker',)),
    ('spoa',     'te-looker consensus refine', 'bioconda::spoa',     ('te_looker',)),
    ('iqtree',   'te-looker family merger',    'bioconda::iqtree',   ('te_looker',)),
    ('TEsorter', 'Refiner_mdl rescue',         'bioconda::tesorter', ('mdl_repeat',)),
    ('barrnap',  'non-TE separation',          'bioconda::barrnap',  ('mdl_repeat',)),
]

CLASSIFIER_TOOLS: List[Tuple[str, str, str, Tuple[str, ...]]] = [
    ('RepeatClassifier', 'Step 5 (classification)', 'bioconda::repeatmodeler', ('combine',)),
]

PERL_MODULES: List[Tuple[str, str, str, Tuple[str, ...]]] = [
    ('Log::Log4perl',         'build_mdl', 'bioconda::perl-log-log4perl',       ('mdl_repeat',)),
    ('Parallel::ForkManager', 'build_mdl', 'bioconda::perl-parallel-forkmanager', ('mdl_repeat',)),
]

PYTHON_MODULES: List[Tuple[str, str, str, Tuple[str, ...]]] = [
    ('Bio',   'Refiner / LTR_detect', 'conda-forge::biopython', ('genome', 'ltr', 'mdl_repeat')),
    ('numpy', 'Refiner',              'conda-forge::numpy',     ('mdl_repeat',)),
    ('scipy', 'Refiner consensus',    'conda-forge::scipy',     ('mdl_repeat',)),
]


class Report:
    """Collected preflight findings. Nothing is raised until every check has run."""

    def __init__(self):
        self.failures: List[str] = []
        self.warnings: List[str] = []

    def fail(self, message: str) -> None:
        self.failures.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    @property
    def ok(self) -> bool:
        return not self.failures


def pending_steps(pipeline) -> set:
    """Pipeline steps that still have to run (checkpoint absent).

    Read directly from the checkpoint filenames rather than via
    PanTEPipeline.check_checkpoint, which logs a line per query.
    """
    steps = getattr(pipeline, 'steps', None)
    if not steps:
        return set()
    legacy = getattr(pipeline, 'legacy_steps', {}) or {}
    output_dir = pipeline.output_dir
    pending = set()
    for step, marker in steps.items():
        names = [marker] + ([legacy[step]] if step in legacy else [])
        if not any(os.path.exists(os.path.join(output_dir, n)) for n in names):
            pending.add(step)
    if getattr(pipeline, 'anno_enabled', False) is False:
        pending.discard('anno')
    return pending


def _check_tools(tools, report: Report, required: bool, pending: set) -> None:
    for name, stage, package, steps in tools:
        # Skip anything only used by steps this run will not execute.
        if pending and not (set(steps) & pending):
            continue
        if shutil.which(name):
            continue
        message = f"{name!r} not found on PATH -- needed by {stage}. Install: {package}"
        (report.fail if required else report.warn)(message)


def _perl_inc_dirs(bin_dir: str) -> List[str]:
    """The extra @INC entries build_mdl adds for RepeatMasker/RepeatModeler modules.

    Kept in step with the BEGIN block in bin/build_mdl: if the two disagree, preflight
    would pass on a setup that then fails at runtime (or the reverse).
    """
    dirs = []
    prefix = os.environ.get('CONDA_PREFIX')
    if prefix:
        import glob as _glob
        dirs += [os.path.join(prefix, 'share', 'RepeatModeler'),
                 os.path.join(prefix, 'share', 'RepeatMasker')]
        dirs += _glob.glob(os.path.join(prefix, 'share', 'RepeatModeler-*'))
        dirs += _glob.glob(os.path.join(prefix, 'share', 'RepeatMasker-*'))
    dirs.append(os.path.join(os.path.dirname(bin_dir), 'share'))
    dirs.append(os.path.join(bin_dir, '..', 'share'))
    return [d for d in dirs if os.path.isdir(d)]


def _check_perl(bin_dir: str, report: Report, pending: set) -> None:
    if pending and 'mdl_repeat' not in pending:
        return
    perl = shutil.which('perl')
    if not perl:
        report.fail("'perl' not found on PATH -- needed by build_mdl. Install: conda-forge::perl")
        return

    inc_args = []
    for directory in [bin_dir] + _perl_inc_dirs(bin_dir):
        inc_args += ['-I', directory]
    for extra in (os.environ.get('PAN_TE_PERL5LIB') or '').split(':'):
        if extra:
            inc_args += ['-I', extra]

    for module, stage, package, _steps in PERL_MODULES:
        proc = subprocess.run([perl] + inc_args + ['-e', f'use {module};'],
                              capture_output=True, encoding='utf-8')
        if proc.returncode != 0:
            report.fail(f"Perl module {module} not available -- needed by {stage}. "
                        f"Install: {package}")

    # RepModelConfig is optional (build_mdl falls back to a PATH lookup for rmblastn),
    # but its absence used to be fatal, so it is worth reporting as a warning.
    proc = subprocess.run([perl] + inc_args + ['-e', 'require RepModelConfig;'],
                          capture_output=True, encoding='utf-8')
    if proc.returncode != 0:
        report.warn("Perl module RepModelConfig not found (RepeatModeler's share/ dir). "
                    "build_mdl will resolve rmblastn from PATH instead. "
                    "Install: bioconda::repeatmodeler")


def _check_python(report: Report, pending: set) -> None:
    import importlib.util
    for module, stage, package, steps in PYTHON_MODULES:
        if pending and not (set(steps) & pending):
            continue
        if importlib.util.find_spec(module) is None:
            report.fail(f"Python module {module!r} not importable -- needed by {stage}. "
                        f"Install: {package}")


def _resolve_classifyte_prefix() -> Tuple[Optional[str], List[str]]:
    """Mirror resolve_classifyte_prefix() in bin/run_Classifier, returning (prefix, tried).

    Kept deliberately identical to the shell version -- a preflight that resolves the
    environment differently from the runtime is worse than none, because it would pass
    on a setup that then dies at the very last step of a multi-hour run. Note the
    ordering: everything that does NOT need the `conda` command comes first, because a
    batch job often has the main env on PATH without conda initialised.
    """
    tried: List[str] = []
    explicit = os.environ.get('CLASSIFYTE_ENV_PREFIX')
    if explicit:
        return explicit, tried

    name = os.environ.get('CLASSIFYTE_ENV_NAME', 'ClassifyTE_env')

    def _ok(prefix: str) -> bool:
        return os.access(os.path.join(prefix, 'bin', 'python'), os.X_OK)

    # 1. Sibling of the active environment -- no conda invocation needed.
    active = os.environ.get('CONDA_PREFIX')
    if active:
        candidate = os.path.join(os.path.dirname(active.rstrip(os.sep)), name)
        tried.append(f"{candidate} (sibling of $CONDA_PREFIX)")
        if _ok(candidate):
            return candidate, tried

    # 2. Configured env directories (clusters often relocate these off $HOME).
    for var in ('CONDA_ENVS_PATH', 'CONDA_ENVS_DIRS'):
        for directory in (os.environ.get(var) or '').split(os.pathsep):
            if not directory:
                continue
            candidate = os.path.join(directory, name)
            tried.append(f"{candidate} (from ${var})")
            if _ok(candidate):
                return candidate, tried

    # 3. Base derived from $CONDA_EXE's location, without running it.
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        candidate = os.path.join(os.path.dirname(os.path.dirname(conda_exe)), 'envs', name)
        tried.append(f"{candidate} (derived from $CONDA_EXE)")
        if _ok(candidate):
            return candidate, tried

    # 4. Only now ask conda, which may not exist in a batch shell at all.
    conda = conda_exe if (conda_exe and os.access(conda_exe, os.X_OK)) else shutil.which('conda')
    if not conda:
        tried.append("(conda command unavailable -- CONDA_EXE unset and 'conda' not on PATH)")
        return None, tried
    try:
        base = subprocess.run([conda, 'info', '--base'], capture_output=True,
                              encoding='utf-8', timeout=120)
        if base.returncode == 0 and base.stdout.strip():
            candidate = os.path.join(base.stdout.strip(), 'envs', name)
            tried.append(f"{candidate} (conda info --base)")
            if _ok(candidate):
                return candidate, tried
        tried.append(f"conda env list (searched for an env named {name})")
        listing = subprocess.run([conda, 'env', 'list'], capture_output=True,
                                 encoding='utf-8', timeout=120)
        for line in (listing.stdout or '').splitlines():
            parts = line.split()
            if parts and parts[0] == name:
                return parts[-1], tried
    except (OSError, subprocess.SubprocessError) as exc:
        tried.append(f"(conda invocation failed: {exc})")
    return None, tried


# What run_Classifier symlinks in and executes; missing any of these is fatal there.
CLASSIFYTE_REQUIRED_DIRS = ('features', 'nodes', 'models')
CLASSIFYTE_REQUIRED_SCRIPTS = ('generate_feature_file_parallel.py', 'evaluate.py')


def _check_classifier(pipeline, report: Report, pending: set) -> None:
    if not getattr(pipeline, 'model_dir', None):
        return
    if pending and 'combine' not in pending:
        return
    _check_tools(CLASSIFIER_TOOLS, report, required=True, pending=pending)

    # The isolated ClassifyTE_env. This is the last thing the pipeline touches, so an
    # absent environment used to surface only after every other stage had completed.
    prefix, tried = _resolve_classifyte_prefix()
    env_name = os.environ.get('CLASSIFYTE_ENV_NAME', 'ClassifyTE_env')
    if not prefix:
        searched = ''.join(f"\n        - {t}" for t in tried)
        report.fail(
            f"Conda environment {env_name!r} not found -- run_Classifier (the final step) "
            f"needs it. Searched:{searched}\n"
            f"    Create it with `conda env create -f environment-classifyte.yml`, or set "
            f"CLASSIFYTE_ENV_PREFIX to its absolute prefix (the directory containing "
            f"bin/python).")
    elif not os.access(os.path.join(prefix, 'bin', 'python'), os.X_OK):
        report.fail(
            f"ClassifyTE environment resolved to {prefix!r} but {prefix}/bin/python is "
            f"not executable. Recreate it with "
            f"`conda env create -f environment-classifyte.yml`, or fix "
            f"CLASSIFYTE_ENV_PREFIX.")

    model_dir = pipeline.model_dir
    if not os.path.isdir(model_dir):
        report.fail(f"--model-dir {model_dir!r} is not a directory. Fetch the ClassifyTE "
                    f"assets with `pan_te-setup-data --classifyte-dir <dir>`")
        return

    # run_Classifier auto-corrects one level up when features/ or models/ is missing,
    # so accept whichever of the two layouts actually holds the assets.
    effective = model_dir
    if not all(os.path.isdir(os.path.join(model_dir, d)) for d in ('features', 'models')):
        parent = os.path.dirname(model_dir.rstrip(os.sep))
        if all(os.path.isdir(os.path.join(parent, d)) for d in ('features', 'models')):
            effective = parent

    missing = [d for d in CLASSIFYTE_REQUIRED_DIRS
               if not os.path.exists(os.path.join(effective, d))]
    if missing:
        report.fail(f"--model-dir {effective!r} is missing {', '.join(missing)}/ -- "
                    f"run_Classifier symlinks these and fails without them. Re-fetch with "
                    f"`pan_te-setup-data --classifyte-dir <dir>`")
    missing_scripts = [s for s in CLASSIFYTE_REQUIRED_SCRIPTS
                       if not os.path.isfile(os.path.join(effective, s))]
    if missing_scripts:
        report.fail(f"--model-dir {effective!r} is missing {', '.join(missing_scripts)} -- "
                    f"this must be a ClassifyTE source checkout, not just its models. "
                    f"Re-fetch with `pan_te-setup-data --classifyte-dir <dir>`")


def _estimate_genome_bp(pipeline) -> Tuple[int, str]:
    """Genome size for the memory check: the .fai when it exists, else the FASTA size.

    Preflight runs before Step 1, so the processed genome's index usually does not exist
    yet. An uncompressed FASTA's byte count is within a couple of percent of its base
    count (headers plus newlines), which is ample for a feasibility gate; the exact check
    runs again from the .fai immediately before dtr.
    """
    exact = pipeline.get_genome_size_bp()
    if exact > 0:
        return exact, 'processed genome .fai'

    source = getattr(pipeline, 'original_genome_file', None)
    if source and os.path.isfile(source):
        if source.endswith(('.gz', '.bz2')):
            return 0, 'compressed input (not estimated)'
        return int(os.path.getsize(source) * 0.98), 'input FASTA size (approximate)'
    return 0, 'unknown'


def _check_memory(pipeline, report: Report, pending: set) -> None:
    # Only Step 4 has this requirement; on a resume that already has families.fasta the
    # expensive discovery is skipped entirely and the estimate does not apply.
    if pending and 'te_looker' not in pending:
        return
    run_dir = os.path.join(pipeline.output_dir, 'te-looker', 'run', 'families.fasta')
    legacy_run = os.path.join(pipeline.output_dir, 'RepGraph', 'run', 'families.fasta')
    for path in (run_dir, legacy_run):
        if os.path.exists(path) and os.path.getsize(path) > 0:
            return

    main = sys.modules.get('__main__')
    estimate = getattr(main, 'estimate_te_looker_peak_bytes', None)
    detect = getattr(main, 'detect_memory_limit', None)
    human = getattr(main, 'human_bytes', None)
    if not (estimate and detect and human):
        return

    genome_bp, origin = _estimate_genome_bp(pipeline)
    if genome_bp <= 0:
        report.warn(f"Genome size could not be determined ({origin}); the te-looker "
                    f"memory check will run at Step 4 instead")
        return

    w = getattr(pipeline, 'te_looker_w', 1)
    # Assume the worst case (a dtr that holds two k-mer buffers) unless the binary is
    # present to probe -- an optimistic guess here would defeat the point of the gate.
    required = estimate(genome_bp, w, True)
    limit = detect(getattr(pipeline, 'memory_limit_mb', 0))
    if limit is None:
        report.warn(f"te-looker needs ~{human(required)} for a {genome_bp/1e9:.2f} Gb "
                    f"genome, but no memory ceiling could be detected")
        return
    if required > limit['bytes']:
        report.fail(
            f"te-looker (Step 4) needs ~{human(required)} for this "
            f"{genome_bp/1e9:.2f} Gb genome at w={w}, but only {human(limit['bytes'])} "
            f"is available (from {limit['source']}). Give the job more memory, or opt "
            f"into k-mer subsampling with --te-looker-w N (which changes the seed set). "
            f"Genome size from: {origin}")


def run_preflight(pipeline, check_only: bool = False) -> bool:
    """Verify the environment for this run. Returns True when it is usable.

    Raises EnvironmentError outside --check-env mode so the pipeline stops before
    creating any output; in --check-env mode the caller decides the exit status.
    """
    import logging
    logger = logging.getLogger('PanTE')

    bin_dir = os.path.dirname(os.path.abspath(__file__))
    report = Report()

    pending = pending_steps(pipeline)
    done = sorted(set(getattr(pipeline, 'steps', {})) - pending - {'anno'})
    if done:
        logger.info(f"Preflight: {', '.join(done)} already checkpointed; checking only "
                    f"what this run still has to do ({', '.join(sorted(pending)) or 'nothing'})")
    else:
        logger.info("Preflight: verifying tools, modules and resources")

    _check_tools(CORE_TOOLS, report, required=True, pending=pending)
    _check_tools(OPTIONAL_TOOLS, report, required=False, pending=pending)
    _check_perl(bin_dir, report, pending)
    _check_python(report, pending)
    _check_classifier(pipeline, report, pending)
    _check_memory(pipeline, report, pending)

    for message in report.warnings:
        logger.warning(f"Preflight: {message}")

    if report.ok:
        logger.info(f"Preflight passed ({len(report.warnings)} warning(s))")
        return True

    detail = '\n'.join(f"  - {message}" for message in report.failures)
    summary = (f"Preflight found {len(report.failures)} blocking problem(s):\n{detail}\n"
               f"Nothing has been run and no output directories were created.")
    if check_only:
        logger.error(summary)
        return False
    raise EnvironmentError(summary)
