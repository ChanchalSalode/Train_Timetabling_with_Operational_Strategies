# main.py

import argparse
import os

from config import OUTPUT_DIR, LOG_DIR

from data_adapters.adapter import load_model_inputs
from model_builder import build_model
from solver import solve_model
from results import print_results

# =========================
# CREATE OUTPUT FOLDERS
# =========================

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

def parse_args():

    parser = argparse.ArgumentParser(
        description="Run simple model with Excel or GTFS/AFC-style inputs."
    )

    parser.add_argument(
        "--input_type",
        choices=["excel", "gtfs_afc"],
        default="excel",
        help="Input adapter to use. Default keeps the original Excel workflow."
    )

    parser.add_argument(
        "--input_path",
        default=None,
        help="Excel workbook path for input_type=excel, or GTFS folder for input_type=gtfs_afc."
    )

    parser.add_argument(
        "--afc_path",
        default=None,
        help="Optional AFC CSV path for input_type=gtfs_afc. Defaults to input_path/afc_taps.csv."
    )

    return parser.parse_args()


def main():

    args = parse_args()

    # =========================
    # LOAD DATA
    # =========================

    data = load_model_inputs(
        input_type=args.input_type,
        input_path=args.input_path,
        afc_path=args.afc_path
    )

    params = {
        'r': data.r,
        'e': data.e,
        'p': data.p
    }

    # =========================
    # BUILD MODEL
    # =========================

    mdl, vars, sets = build_model(params)

    # =========================
    # SOLVE MODEL
    # =========================

    solve_model(mdl)

    # =========================
    # PRINT RESULTS
    # =========================

    print_results(
        mdl,
        vars,
        sets
    )


if __name__ == "__main__":

    main()
