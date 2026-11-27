"""Run PPR_model_v9_10.py with the input folder and the copy folder given on the command line.

    python code/run_model.py --data-root /path/to/inputs [--copy-dir /path/to/copies] [model options ...]

--data-root   folder with the input layers (see INPUT_DATA.md); environment variable PPR_DATA_ROOT also works
--copy-dir    folder that receives a copy of the tables and figures of each run (default: ./copies/v9_10)
Every other option goes to the model unchanged (python code/PPR_model_v9_10.py --help lists them).
Results are written to ./results/ppr_v9_10_review<suffix>/ under the folder you run from.
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def _pop(argv, flag, env):
    val = os.environ.get(env)
    if flag in argv:
        i = argv.index(flag)
        val = argv[i + 1]
        del argv[i:i + 2]
    return val


def main(argv):
    argv = list(argv)
    data = _pop(argv, "--data-root", "PPR_DATA_ROOT")
    copy = _pop(argv, "--copy-dir", "PPR_COPY_DIR")
    import PPR_model_v9_10 as M
    C = M.Config
    C.CWD = Path.cwd()
    C.OUT_DIR = C.CWD / "results" / "ppr_v9_10_review"
    if data:
        C.DATA_ROOT = Path(data).expanduser().resolve()
        C.ENV_DIR = C.DATA_ROOT / "PPR_env"
    elif not Path(C.DATA_ROOT).exists():
        sys.exit("give the input folder with --data-root (or set PPR_DATA_ROOT); see INPUT_DATA.md")
    C.DL_DIR = Path(copy).expanduser().resolve() if copy else C.CWD / "copies" / "v9_10"
    M.v910_cli(argv)


if __name__ == "__main__":
    main(sys.argv[1:])
