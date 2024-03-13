import os
import pytest
import yaml
import numpy as np
from slimflow import GridRuns


@pytest.fixture
def load_config():
    config_path = os.path.join(os.path.dirname(__file__), "data",
                               "test_config.yml")
    with open(config_path, "r") as file:
        return yaml.safe_load(file)


def test_slim_runs_entry(load_config):
    slim_runs = GridRuns(load_config)

    assert slim_runs.name == "bgs"
    assert slim_runs.seed == 42
    assert isinstance(slim_runs.rng, np.random.Generator)
    assert len(slim_runs.variables) == 4


def test_generate_targets_basic(load_config):
    slim_runs = GridRuns(load_config)
    targets = slim_runs.generate_targets()

    # 4 params x 2 reps x 2 suffices = 16
    assert targets.shape[0] == 16


def test_directory_creation(load_config, tmp_path):
    """
    Test that directories are created when `create=True`.
    Uses pytest's `tmp_path` fixture to avoid creating directories in the
    actual project structure.
    """
    load_config['dir'] = str(tmp_path)  # Use temporary path
    slim_runs = GridRuns(load_config)
    _ = list(slim_runs.generate_targets())

    # Check if directories were created
    expected_dir = tmp_path / "bgs"
    assert expected_dir.exists()
