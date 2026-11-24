#!/usr/bin/env python3
"""Compare the output of a Python API test using the runtest rules."""
import os
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import yaml

from py123d.scripts.prescriptions.local_run import LocalRun
from py123d.scripts.yamlc.yaml_config import ConfigCase
from py123d.loggers import printf


def main(test_dir, script, prepare=False):
    root = Path(test_dir).resolve()
    name = Path(script).name
    with (root / 'py_config.yaml').open(encoding='utf-8') as config_file:
        config = yaml.safe_load(config_file)
    common = config.get('common_config') or {}
    case_config = next((case for case in config['test_cases']
                        if name in case['files']), None)
    if case_config is None:
        raise ValueError(f'No configuration for {name}')
    options = dict(common)
    options.update(case_config)
    options.update(files=name, proc=1)
    case = ConfigCase(options, SimpleNamespace(root=str(root)))
    if prepare:
        if Path(case.fs.output).exists():
            shutil.rmtree(case.fs.output)
        return 0

    os.environ.setdefault('FLOW123D_DIR', str(Path(__file__).resolve().parents[3]))
    printf.init_logger()
    comparisons = LocalRun(case).create_comparisons()
    comparisons.progress = False
    comparisons.start()
    comparisons.join()
    if (not comparisons.is_over() or any(not thread.is_over() for thread in comparisons.threads)
            or comparisons.returncode not in (0, None)):
        print(f'Comparison failed; see {case.fs.ndiff_log}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1], sys.argv[2], '--prepare' in sys.argv[3:]))
    except Exception as error:
        sys.exit(f'Comparison error: {error}')
