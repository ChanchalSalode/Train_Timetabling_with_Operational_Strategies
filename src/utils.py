# src/utils.py
import pandas as pd
import yaml
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def load_params(path: str = None):
    if path is None:
        path = ROOT / "data" / "parameters.yaml"
    with open(path, "r") as f:
        return yaml.safe_load(f)

def load_network(path: str = None):
    if path is None:
        path = ROOT / "data" / "network_topology.csv"
    return pd.read_csv(path)

def load_od_static(path: str = None):
    if path is None:
        path = ROOT / "data" / "OD_static.csv"
    return pd.read_csv(path)

def load_od_time_dependent(path: str = None):
    if path is None:
        path = ROOT / "data" / "OD_time_dependent.csv"
    return pd.read_csv(path)

def load_station_profile(station_csv_path):
    return pd.read_csv(station_csv_path)
