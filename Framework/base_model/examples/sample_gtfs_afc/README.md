# Sample GTFS/AFC Input

Run from `simple_model_final`:

```bash
python main.py --input_type gtfs_afc --input_path examples/sample_gtfs_afc
```

The adapter reads GTFS files plus `afc_taps.csv`, then maps them to the internal schema: `r`, `e`, `p`, `phi`, `P_backlog`, `station_map`, `line_direction`, and `service_candidates`.
