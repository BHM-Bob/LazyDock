'''
Date: 2026-08-28
Description: mdp template paths for ABFE, migrated from BindFlow bindflow.mdp.templates.
             Core logic migrated from BindFlow (https://github.com/IFMIFMIF/BindFlow)
'''
from pathlib import Path

_TEMPLATE_ROOT = Path(__file__).resolve().parent / 'mdp_templates'


class _TemplatePath(object):
    def __init__(self, path: Path):
        self.path = path
        self.equi = path / 'equi'
        self.fep = path / 'fep'

    def __str__(self):
        return str(self.path)


class _LigandPath(_TemplatePath):
    pass


class _ComplexMembranePath(_TemplatePath):
    pass


class _ComplexSolublePath(_TemplatePath):
    pass


class _ComplexPath(object):
    def __init__(self):
        self.membrane = _ComplexMembranePath(_TEMPLATE_ROOT / 'complex' / 'membrane')
        self.soluble = _ComplexSolublePath(_TEMPLATE_ROOT / 'complex' / 'soluble')


class TemplatePath(object):
    ligand = _LigandPath(_TEMPLATE_ROOT / 'ligand')
    complex = _ComplexPath()


if __name__ == '__main__':
    print(TemplatePath.ligand.equi)
    print(TemplatePath.complex.soluble.fep)
