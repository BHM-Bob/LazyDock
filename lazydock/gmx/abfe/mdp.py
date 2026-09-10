'''
Date: 2026-08-28
LastEditors: BHM-Bob 2262029386@qq.com
Description: MDP helpers for ABFE, migrated from BindFlow bindflow.mdp.mdp.
             Core logic migrated from BindFlow (https://github.com/IFMIFMIF/BindFlow)
'''
import json
from pathlib import Path
from typing import List, Union
import os

PathLike = Union[str, Path]

_MDP_PARAM_DEFAULT = {
    "integrator": "steep",
    "emtol": "1000.0",
    "nsteps": "5000",
    "nstlist": "200",
    "cutoff-scheme": "Verlet",
    "rlist": "1.0",
    "vdwtype": "Cut-off",
    "vdw-modifier": "Potential-shift-Verlet",
    "rvdw-switch": "0",
    "rvdw": "1.0",
    "coulombtype": "pme",
    "rcoulomb": "1.0",
    "epsilon-r": "1",
    "epsilon-rf": "1",
    "constraints": "h-bonds",
    "constraint-algorithm": "LINCS"
}


def list_if_file(path: PathLike, ext: str = None) -> list:
    """List files in path, optionally filtered by extension."""
    path = Path(path)
    if ext:
        return [p for p in path.iterdir() if p.is_file() and p.suffix == ext]
    return [p for p in path.iterdir() if p.is_file()]


class MDP:
    """Base class to work with MDP files."""

    def __init__(self, **kwargs):
        self.parameters = dict()
        self._set_default_parameters()
        self.set_parameters(**kwargs)

    def _set_default_parameters(self):
        self.parameters = dict(_MDP_PARAM_DEFAULT)

    def set_parameters(self, **kwargs):
        kwargs = {key.replace('_', '-'): value for key, value in kwargs.items()}
        self.parameters.update(kwargs)

    def from_file(self, template_filename, clean_current_parameters=True):
        with open(template_filename, 'r') as f:
            lines = f.readlines()
            if clean_current_parameters:
                self.parameters = {}
            for line in lines:
                if line.startswith(';') or line.startswith('#'):
                    continue
                tokens = line.strip().split('=', 1)
                if len(tokens) != 2:
                    continue
                parameter_name = tokens[0].strip().replace('_', '-')
                parameter_value = tokens[1].strip()
                self.parameters[parameter_name] = parameter_value
        return self

    def to_string(self):
        s = ''
        for parameter_name, parameter_value in self.parameters.items():
            s += f'{parameter_name:<40} = {parameter_value}\n'
        return s

    def write(self, filename: str):
        with open(filename, 'w') as f:
            f.write(self.to_string())

    def __repr__(self):
        return f"{self.__class__.__name__}({json.dumps(self.parameters, indent=4)})"


class StepMDP(MDP):
    """StepMDP works with the mdp templates under abfe/mdp_templates.
    It defines ``set_new_step`` to switch between the step mdp files."""

    def __init__(self, step: str = None, step_path: PathLike = None, **kwargs):
        super().__init__(**kwargs)
        self.step = step
        self.step_path = Path(step_path)
        if self.step:
            self._from_archive()

    def set_new_step(self, step):
        self._from_archive(explicit_step=step)
        return self

    def _from_archive(self, explicit_step: str = None):
        if explicit_step:
            self.step = explicit_step
        valid_steps = [step.stem for step in list_if_file(self.step_path, ext='.mdp')]
        if self.step not in valid_steps:
            raise ValueError(f"name = {self.step} is not a valid step mdp, must be one of: {valid_steps}")
        self.from_file(self.step_path / f"{self.step}.mdp")


def make_fep_dir_structure(sim_dir: PathLike, template_dir: PathLike, lambda_values: List[float],
                           lambda_type: str, sys_type: str, dt_max: float, mdp_extra_kwargs: dict = None,
                           couple_intramol: str = 'yes', fep_rlist: float = None,
                           fep_cutoff: float = None):
    """Create the simulation directory structure:
    ``{sim_dir}/simulation/{lambda_type}.{i}/{step}/{step}.mdp``

    Where i is the init-lambda-state and step the simulation step name.
    Migrated from BindFlow's make_fep_dir_structure.

    ``couple_intramol``: 'yes' (default, small molecules) or 'no' (flexible
    ligands such as peptides).  With 'no', the intramolecular non-bonded
    interactions of the ligand are kept fully on during the decoupling
    (they cancel in the thermodynamic cycle), which prevents long flexible
    chains from unfolding in the fully-decoupled (lambda=0 vdw / lambda=0
    coul) windows.

    GROMACS requires that all perturbed, excluded (1-4) pairs stay within
    rlist.  A large flexible ligand (peptide diameter can exceed 2 nm) thus
    needs ``rlist`` and ``table-limit`` covering its maximum internal
    distance; the default BindFlow rlist=1.2 nm is fine for small molecules
    but would make mdrun abort with "perturbed, excluded non-bonded pair
    interactions beyond the pair-list cut-off".  For this reason the
    flexible-peptide mode overrides rlist/table-limit to cover the ligand
    (see the engine's ``couple_intramol='no'`` path).

    ``fep_rlist``: optional explicit rlist value (nm) for peptide FEP
    windows; if None the built-in default (3.1) is used.  Larger values
    cost performance, smaller values may trigger the excluded-pair fatal.
    """
    sim_dir = Path(sim_dir)
    template_dir = Path(template_dir)
    valid_lambda_types = ["vdw", "coul", "bonded"]
    valid_sys_types = ['ligand', 'complex']
    if lambda_type not in valid_lambda_types:
        raise ValueError(f"Non valid lambda_type = {lambda_type}. Must be one of {valid_lambda_types}")
    if sys_type not in valid_sys_types:
        raise ValueError(f"Non valid sys_type = {sys_type}. Must be one of {valid_sys_types}")

    input_mdp = [step.name for step in list_if_file(template_dir / f"{lambda_type}", ext='.mdp')]
    lambda_range_str = " ".join(map(str, lambda_values))
    mdp_template = StepMDP(step_path=template_dir / lambda_type)
    for mdp_file in input_mdp:
        step = Path(mdp_file).stem
        mdp_template.set_new_step(step)
        if 'dt' in mdp_template.parameters:
            if float(mdp_template.parameters['dt'].split(';')[0]) > dt_max:
                mdp_template.set_parameters(dt=dt_max)
        mdp_template.set_parameters(couple_intramol=couple_intramol)
        # NOTE: nstlist 已直接在模板中设为 200 (00_min 保持 1)。用户 --mdp-extra 可覆盖 (下方 mdp_extra_kwargs)。
        if couple_intramol == 'no':
            # couple-intramol=no: 配体分子内非键保持全开(避免解耦窗口链展开),
            # 但 GROMACS 2021+ 要求所有被扰动的排除对(1-4/跨环)距离 <= rlist,
            # 因此 rlist 必须 >= 配体最大内部距离.
            #   - 实测(GMX 2026.3, 16肽solv-d=2.0盒, a=6.494 nm):
            #     * couple-intramol=no 时 tpr 中 LIG 内 238 个原子全部互为排除对
            #       (excls 56644 元素 = 238^2), exclusionchecker 对它们逐一检查
            #     * 02_npt.gro (posres 压实) 配体最大内部距离 = 2.567 nm
            #     * 去 posres 自由涨落后, rlist=2.75 余量仅 0.18 nm, 首个
            #       neighbor search 即 fatal "1 perturbed excluded pair"
            #   - rlist=3.1: 仍 < 半最短盒矢 3.247 nm (dodecahedron 盒),
            #     给热涨落 ~0.53 nm 余量；可用 fep_rlist 覆盖
            # - verlet-buffer-tolerance=-1: 禁用自动 Verlet 缓冲, 固定 rlist
            #   (否则 grompp 会把 rlist 重算为 rcoulomb+buffer, 导致 rlist < 配体直径,
            #    mdrun 报 "perturbed excluded pairs beyond pair-list cutoff")
            # - rcoulomb/rvdw 保持标准值 1.0: rlist 只负责邻居列表范围,
            #   实际相互作用截断仍由 rcoulomb/rvdw 决定(PME 实空间), 两者可独立
            rlist = 3.1 if fep_rlist is None else float(fep_rlist)
            # NOTE : couple-intramol=no 时, 分子内非键全程保持开启
            # (non-perturbed), grompp 会要求所有 1-4/环状排除对距离 <= max(rlist,
            # rvdw, rcoulomb)。肽/环肽柔性链在平衡/生产阶段 1-4 排除对可伸展超过
            # 1.0 nm (实测XX肽 6TYR-16TYR 排除对 1.074 nm), 而 rlist 只负责
            # 邻居列表, 不影响 grompp 的排除对检查 -> 必须同步提高 rcoulomb/rvdw。
            # 默认 1.4 (覆盖典型伸长 ~1.1 nm + 涨落余量), 可用 fep_cutoff 覆盖。
            cutoff = 1.4 if fep_cutoff is None else float(fep_cutoff)
            mdp_template.set_parameters(rlist=rlist,
                                        rcoulomb=cutoff,
                                        rvdw=cutoff,
                                        verlet_buffer_tolerance=-1)
        if mdp_extra_kwargs:
            try:
                mdp_template.set_parameters(**mdp_extra_kwargs[lambda_type][step])
            except KeyError:
                pass
        mdp_template.set_parameters(**{f"{lambda_type}-lambdas": lambda_range_str})
        # Set to 1 all the bonded-lambdas in case of vdw and coul for the complex
        if sys_type.lower() == 'complex' and lambda_type in ['vdw', 'coul']:
            mdp_template.set_parameters(**{"bonded-lambdas": " ".join(map(str, len(lambda_values) * [1]))})
        for i in range(len(lambda_values)):
            (sim_dir / f"simulation/{lambda_type}.{i}/{step}").mkdir(exist_ok=True, parents=True)
            mdp_template.set_parameters(**{"init-lambda-state": i})
            mdp_template.write(sim_dir / f"simulation/{lambda_type}.{i}/{step}/{step}.mdp")


def get_number_of_frames(input_mdp: PathLike):
    loaded_mdp_params = MDP().from_file(input_mdp).parameters
    return -(int(loaded_mdp_params['nsteps'].split(';')[0]) // -int(loaded_mdp_params['nstxout-compressed'].split(';')[0]))


if __name__ == '__main__':
    pass
