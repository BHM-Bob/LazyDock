'''
Date: 2026-08-28
LastEditors: BHM-Bob 2262029386@qq.com

Core logic migrated from BindFlow (https://github.com/IFMIFMIF/BindFlow)
'''
from pathlib import Path
import inspect


def _data_dir(dataDir=None) -> Path:
    """return the data directory of the abfe subpackage (or a subdirectory)."""
    abfe_dir = Path(inspect.getfile(_data_dir)).resolve().parent
    if dataDir:
        return abfe_dir / 'data' / dataDir
    return abfe_dir / 'data'


def gmx_ff_data_dir() -> Path:
    return _data_dir('gmx_ff')


def gmx_water_models_data_dir() -> Path:
    return _data_dir('gmx_water_models')


if __name__ == '__main__':
    print(gmx_ff_data_dir())
    print(gmx_water_models_data_dir())
