import argparse
import os
import sys

from config import PARAMETERS

from data_adapters.adapter import load_model_inputs

from sets import create_sets

from model_builder import build_model

from solver import solve_model


class Tee:

    def __init__(
        self,
        *streams
    ):

        self.streams = streams

    def write(
        self,
        data
    ):

        for stream in self.streams:

            stream.write(data)

            stream.flush()

    def flush(self):

        for stream in self.streams:

            stream.flush()


def main():
    args = parse_args()

    os.makedirs(
        PARAMETERS['OUTPUT_DIR'],
        exist_ok=True
    )

    os.makedirs(
        PARAMETERS['LOG_DIR'],
        exist_ok=True
    )

    open(
        PARAMETERS['GUROBI_LOG_FILE'],
        'w',
        encoding='utf-8'
    ).close()

    with open(
        PARAMETERS['TERMINAL_OUTPUT_FILE'],
        'w',
        encoding='utf-8'
    ) as terminal_log:

        original_stdout = sys.stdout

        original_stderr = sys.stderr

        sys.stdout = Tee(
            original_stdout,
            terminal_log
        )

        sys.stderr = Tee(
            original_stderr,
            terminal_log
        )

        try:

            run(args)

        finally:

            sys.stdout = original_stdout

            sys.stderr = original_stderr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run classical epsilon model with Excel or GTFS/AFC-style inputs."
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


def run(args):

    # =========================================
    # Input File
    # =========================================

    print("\nLoading input data...")

    model_input = load_model_inputs(
        input_type=args.input_type,
        input_path=args.input_path,
        afc_path=args.afc_path
    )

    data = {
        'r': model_input.r,
        'e': model_input.e,
        'p': model_input.p
    }

    print("Data loaded successfully.")

    # =========================================
    # Create Sets
    # =========================================

    print("\nCreating sets...")

    sets = create_sets(
        PARAMETERS
    )

    print("Sets created successfully.")

    # =========================================
    # Build Model
    # =========================================

    print("\nBuilding optimization model...")

    mdl, vars = build_model(

        sets=sets,

        params=PARAMETERS,

        data=data

    )

    print("Model built successfully.")

    # =========================================
    # Solve Model
    # =========================================

    print("\nOptimizing model...")

    solve_model(

        mdl=mdl,

        vars=vars,

        sets=sets,

        params=PARAMETERS

    )

    print("\n=================================")
    print("PROGRAM COMPLETED")
    print("=================================")


# =============================================
# Run Main
# =============================================

if __name__ == "__main__":

    main()
