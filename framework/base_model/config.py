# config.py

import os

# =========================
# MODEL PARAMETERS
# =========================

S_up = 8
S_dn = 16

K_up = 6
K_dn = 12

h_max = 360
h_min = 90

delta_min = 135

RS = 5
C = 250

M_time = 100000.0
M_flow = 5000.0

TIME_LIMIT = 14400

# =========================
# FILE PATHS
# =========================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

DATA_FILE = os.path.join(
    BASE_DIR,
    "data",
    "Input_Parameters.xlsx"
)

OUTPUT_DIR = os.path.join(
    BASE_DIR,
    "outputs"
)

LOG_DIR = os.path.join(
    OUTPUT_DIR,
    "logs"
)

SOLUTION_FILE = os.path.join(
    OUTPUT_DIR,
    "solution.txt"
)

GUROBI_LOG_FILE = os.path.join(
    LOG_DIR,
    "gurobi.log"
)
