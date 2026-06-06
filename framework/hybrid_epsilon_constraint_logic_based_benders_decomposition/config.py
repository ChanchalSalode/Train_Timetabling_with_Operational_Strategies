# config.py

import os

# =========================================
# MODEL PARAMETERS
# =========================================

S_U = 8
S_D = 16

K_U = 9
K_D = 17

h_max = 360
h_min = 90

delta_min = 135

RS = 14

C = 250

M_time = 100000.0
M_flow = 5000.0

# Big-M
M = M_flow

# =========================================
# SOLVER PARAMETERS
# =========================================

TIME_LIMIT = 14400

OUTPUT_FLAG = 1

LAZY_CONSTRAINTS = 1

MAX_Y = 12

# =========================================
# PATH SETTINGS
# =========================================

BASE_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

DATA_DIR = os.path.join(
    BASE_DIR,
    "data"
)

OUTPUT_DIR = os.path.join(
    BASE_DIR,
    "outputs"
)

LOG_DIR = os.path.join(
    OUTPUT_DIR,
    "logs"
)

# =========================================
# INPUT FILE
# =========================================

DATA_FILE = os.path.join(
    DATA_DIR,
    "Input_Parameters.xlsx"
)

# =========================================
# OUTPUT FILES
# =========================================

SOLUTION_FILE = os.path.join(
    OUTPUT_DIR,
    "solution.txt"
)

GUROBI_LOG_FILE = os.path.join(
    LOG_DIR,
    "gurobi.log"
)

# =========================================
# EXCEL SHEETS
# =========================================

PURE_RUNNING_TIME_SHEET = "Pure_running_time"

ACCELERATION_TIME_SHEET = "Acceleration_time"

DECELERATION_TIME_SHEET = "Deceleration_time"

DWELLING_TIME_SHEET = "Dwelling_time"

DEMAND_SHEET = "Demand_30(7.30-8.00am)mint_peak"

# =========================================
# MODEL NAME
# =========================================

MODEL_NAME = "TRB_Santiago_16station"

# =========================================
# CREATE OUTPUT DIRECTORIES
# =========================================

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)

os.makedirs(
    LOG_DIR,
    exist_ok=True
)

# =========================================
# PARAMETERS DICTIONARY
# =========================================

PARAMETERS = {

    # Station Parameters
    'S_U': S_U,
    'S_D': S_D,

    # Train Parameters
    'K_U': K_U,
    'K_D': K_D,

    # Headway Parameters
    'h_max': h_max,
    'h_min': h_min,

    # Transfer Parameters
    'delta_min': delta_min,

    # Rolling Stock
    'RS': RS,

    # Capacity
    'C': C,

    # Big-M Parameters
    'M_time': M_time,
    'M_flow': M_flow,
    'M': M,

    # Solver Parameters
    'TIME_LIMIT': TIME_LIMIT,
    'OUTPUT_FLAG': OUTPUT_FLAG,
    'LAZY_CONSTRAINTS': LAZY_CONSTRAINTS,

    # Epsilon Parameter
    'MAX_Y': MAX_Y,

    # File Paths
    'DATA_FILE': DATA_FILE,
    'OUTPUT_DIR': OUTPUT_DIR,
    'LOG_DIR': LOG_DIR,
    'SOLUTION_FILE': SOLUTION_FILE,
    'GUROBI_LOG_FILE': GUROBI_LOG_FILE,

    # Excel Sheets
    'PURE_RUNNING_TIME_SHEET': PURE_RUNNING_TIME_SHEET,
    'ACCELERATION_TIME_SHEET': ACCELERATION_TIME_SHEET,
    'DECELERATION_TIME_SHEET': DECELERATION_TIME_SHEET,
    'DWELLING_TIME_SHEET': DWELLING_TIME_SHEET,
    'DEMAND_SHEET': DEMAND_SHEET,

    # Model Name
    'MODEL_NAME': MODEL_NAME

}
