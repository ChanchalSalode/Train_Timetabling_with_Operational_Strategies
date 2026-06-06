# config.py

import os

# =========================
# MODEL PARAMETERS
# =========================

S_U = 8
S_D = 16

K_U = 6
K_D = 12

h_max = 360
h_min = 90

delta_min = 135

RS = 14

C = 250

M_time = 100000.0
M_flow = 5000.0

# Big-M
M = M_flow

# =========================
# SOLVER PARAMETERS
# =========================

TIME_LIMIT = 14400

OUTPUT_FLAG = 1

# =========================
# BASE DIRECTORY
# =========================

BASE_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

# =========================
# DATA FILE PATH
# =========================

INPUT_FILE = os.path.join(
    BASE_DIR,
    "data",
    "Input_Parameters.xlsx"
)

# =========================
# OUTPUT DIRECTORIES
# =========================

OUTPUT_DIR = os.path.join(
    BASE_DIR,
    "outputs"
)

LOG_DIR = os.path.join(
    OUTPUT_DIR,
    "logs"
)

# =========================
# OUTPUT FILES
# =========================

SOLUTION_FILE = os.path.join(
    OUTPUT_DIR,
    "solution.txt"
)

TERMINAL_OUTPUT_FILE = os.path.join(
    OUTPUT_DIR,
    "terminal_writeup.txt"
)

GUROBI_LOG_FILE = os.path.join(
    LOG_DIR,
    "gurobi.log"
)

# =========================
# EXCEL SHEET NAMES
# =========================

RUNNING_TIME_SHEET = "Running_time"

DWELLING_TIME_SHEET = "Dwelling_time"

DEMAND_SHEET = "Demand_30(7.30-8.00am)mint"

# =========================
# MODEL NAME
# =========================

MODEL_NAME = "TRB_Santiago_16station"

# =========================
# PARAMETER DICTIONARY
# =========================

PARAMETERS = {

    'S_U': S_U,
    'S_D': S_D,

    'K_U': K_U,
    'K_D': K_D,

    'h_max': h_max,
    'h_min': h_min,

    'delta_min': delta_min,

    'RS': RS,

    'C': C,

    'M_time': M_time,
    'M_flow': M_flow,

    'M': M,

    'TIME_LIMIT': TIME_LIMIT,

    'OUTPUT_FLAG': OUTPUT_FLAG,

    'BASE_DIR': BASE_DIR,

    'INPUT_FILE': INPUT_FILE,

    'OUTPUT_DIR': OUTPUT_DIR,

    'LOG_DIR': LOG_DIR,

    'SOLUTION_FILE': SOLUTION_FILE,

    'TERMINAL_OUTPUT_FILE': TERMINAL_OUTPUT_FILE,

    'GUROBI_LOG_FILE': GUROBI_LOG_FILE,

    'LOG_TO_CONSOLE': 0,

    'RUNNING_TIME_SHEET': RUNNING_TIME_SHEET,

    'DWELLING_TIME_SHEET': DWELLING_TIME_SHEET,

    'DEMAND_SHEET': DEMAND_SHEET,

    'MODEL_NAME': MODEL_NAME

}
