# config.py

import os

# =========================================
# NETWORK CONFIGURATION
# =========================================

S_U = 8
S_D = 16

K_U = 6
K_D = 12


# =========================================
# OPERATIONAL PARAMETERS
# =========================================

h_max = 360
h_min = 90

delta_min = 135

RS = 14

C = 250


# =========================================
# BIG-M PARAMETERS
# =========================================

M_time = 100000.0
M_flow = 5000.0

M = M_flow


# =========================================
# TIME SETTINGS
# =========================================

H = 27000                 # 07:30:00 in seconds

T_WINDOW = 1800           # 30 minutes

BIN = 60                  # 1-minute bins


# =========================================
# DISCRETIZATION
# =========================================

_times = list(
    range(
        0,
        T_WINDOW + 1,
        BIN
    )
)

_n_bins = len(_times) - 1


# =========================================
# NUMERICAL TOLERANCE
# =========================================

eps = 1e-9


# =========================================
# PREPROCESSOR SETTINGS
# =========================================

PREPROC_HISTORY = []

PREPROC_MAX_ITERS = 3

PREPROC_TOL = 1e-3

PREPROC_MODE = 'soft'

PREPROC_ALPHA = 0.1

PREPROC_MIN_PENALTY = 0.1

PREPROC_BETA = 1.0


# =========================================
# QUEUE MODEL PARAMETERS
# =========================================

MU_BOARD = 0.5

DWELL_BASE_SEC = 20.0

DWELL_PER_PAX_SEC = 0.25


# =========================================
# FILE PATHS
# =========================================

DATA_DIR = 'data'


def resolve_data_file(filename):
    """Prefer data/filename, but keep root-level files backward-compatible."""

    data_path = os.path.join(DATA_DIR, filename)

    if os.path.exists(data_path):

        return data_path

    if os.path.exists(filename):

        return filename

    return data_path


MAIN_DATA_FILE = resolve_data_file(
    'Input_Parameters.xlsx'
)

TIME_DEPENDENT_FILE = resolve_data_file(
    'time_dependent_16stations_30mint.xlsx'
)


# =========================================
# OUTPUT FILES
# =========================================

METRICS_OUTPUT_FILE = (
    'Qmetrics_comparison_'
    'QTRB_30_mint_Model_1_'
    'Off_peak_Morning_14_trains.txt'
)
