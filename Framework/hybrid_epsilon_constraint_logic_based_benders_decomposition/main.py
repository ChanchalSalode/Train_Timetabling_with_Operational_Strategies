# main.py

import argparse

from config import *

from data_adapters.adapter import load_model_inputs

from sets import create_sets

from model_builder import build_model

from epsilon_controller import run_epsilon_loop


def parse_args():

    parser = argparse.ArgumentParser(
        description="Run hybrid epsilon model with Excel or GTFS/AFC-style inputs."
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


# =========================================
# MAIN FUNCTION
# =========================================

def main():
    args = parse_args()

    print("\n=================================")
    print("HYBRID EPSILON + BENDERS MODEL")
    print("=================================")

    # =====================================
    # Load Input Data
    # =====================================

    print("\nLoading input data...")

    model_input = load_model_inputs(
        input_type=args.input_type,
        input_path=args.input_path,
        afc_path=args.afc_path
    )

    data = {
        'pr': model_input.pr,
        'ac': model_input.ac,
        'dc': model_input.dc,
        'p': model_input.p,
        'e': model_input.e
    }

    print("Data loaded successfully.")

    # =====================================
    # Parameters Dictionary
    # =====================================

    params = PARAMETERS
    # =====================================
    # Create Sets
    # =====================================

    print("\nCreating sets...")

    sets = create_sets(params)

    print("Sets created successfully.")
    # =====================================
    # Build Optimization Model
    # =====================================

    print("\nBuilding optimization model...")

    mdl, vars = build_model(
        sets=sets,
        data=data,
        params=params
    )
    mdl._params_data = params
    mdl._input_data = data
    print("Model built successfully.")

    # =====================================
    # Run Hybrid ε-Constraint Controller
    # =====================================

    print("\nStarting ε-constraint loop...")

    run_epsilon_loop(
        mdl,
        vars,
        sets,
        params,
        data
    )

    # =====================================
    # Finished
    # =====================================

    print("\n=================================")
    print("PROGRAM COMPLETED")
    print("=================================")


# =========================================
# DRIVER
# =========================================

if __name__ == "__main__":

    main()
