# most of the code is from PepGLAD
import io
import logging
import os
import time
from typing import List, Optional, Union

import openmm
import pdbfixer
from openmm import app as openmm_app
from openmm import unit

ENERGY = unit.kilojoules_per_mole
LENGTH = unit.nanometer  # pyright: ignore[reportAttributeAccessIssue]


class ForceFieldMinimizer(object):

    def __init__(self, stiffness=10.0, max_iterations=0, tolerance=10, platform='CUDA', constraints=openmm_app.HBonds):
        super().__init__()
        self.stiffness = stiffness * ENERGY/(LENGTH ** 2) # coefficient for restraints force, unit: kJ*mol^-1*nm^-2  # pyright: ignore[reportOperatorIssue]
        self.max_iterations = max_iterations
        self.tolerance = tolerance * ENERGY/LENGTH  # pyright: ignore[reportOperatorIssue]
        assert platform in ('CUDA', 'CPU')
        self.platform = platform
        self.constraints = constraints

    def _fix(self, pdb_str):
        fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms(seed=0)
        fixer.addMissingHydrogens()

        out_handle = io.StringIO()
        openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
        return out_handle.getvalue()

    def _get_pdb_string(self, topology, positions):
        with io.StringIO() as f:
            openmm_app.PDBFile.writeFile(topology, positions, f, keepIds=True)
            return f.getvalue()
        
    def _minimize(self, pdb_str: str, restrain_chain: Optional[Union[str, List[str]]] = None,
                  restrain_backbone: bool = False):
        pdb = openmm_app.PDBFile(io.StringIO(pdb_str))

        force_field = openmm_app.ForceField("charmm36.xml") # referring to http://docs.openmm.org/latest/userguide/application/02_running_sims.html
        system = force_field.createSystem(pdb.topology, constraints=self.constraints)

        # Add constraints to restrain_chain
        restrain_chain = restrain_chain or []
        restrain_name = ['N', 'CA', 'C', 'O'] if restrain_backbone else []
        if restrain_chain:
            force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
            force.addGlobalParameter("k", self.stiffness)
            for p in ["x0", "y0", "z0"]:
                force.addPerParticleParameter(p)
            
            for i, a in enumerate(pdb.topology.atoms()):
                if a.residue.chain.id in restrain_chain and (not restrain_backbone or a.name in restrain_name):
                    force.addParticle(i, pdb.positions[i])

            system.addForce(force)

        # Set up the integrator and simulation
        integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
        platform = openmm.Platform.getPlatformByName(self.platform)
        simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)

        # Perform minimization
        ret = {}
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        ret["einit"] = state.getPotentialEnergy().value_in_unit(ENERGY)
        ret["posinit"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)

        simulation.minimizeEnergy(maxIterations=self.max_iterations, tolerance=self.tolerance)

        state = simulation.context.getState(getEnergy=True, getPositions=True)
        ret["efinal"] = state.getPotentialEnergy().value_in_unit(ENERGY)
        ret["pos"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
        ret["min_pdb"] = self._get_pdb_string(simulation.topology, state.getPositions())

        return ret['min_pdb'], ret
    
    def _add_energy_remarks(self, pdb_str, ret):
        pdb_lines = pdb_str.splitlines()
        pdb_lines.insert(1, "REMARK   1  FINAL ENERGY:   {:.3f} KCAL/MOL".format(ret['efinal']))
        pdb_lines.insert(1, "REMARK   1  INITIAL ENERGY: {:.3f} KCAL/MOL".format(ret['einit']))
        return "\n".join(pdb_lines)

    def __call__(self, pdb_str, out_path, restrain_chain: Optional[Union[str, List[str]]] = None, restrain_backbone: bool = False, return_info=True):
        if '\n' not in pdb_str and pdb_str.lower().endswith(".pdb"):
            with open(pdb_str) as f:
                pdb_str = f.read()

        pdb_fixed = self._fix(pdb_str)
        pdb_min, ret = self._minimize(pdb_fixed, restrain_chain, restrain_backbone)
        pdb_min = self._add_energy_remarks(pdb_min, ret)
        if out_path and os.path.exists(out_path):
            with open(out_path, 'w') as f:
                f.write(pdb_min)
        if return_info:
            return pdb_min, ret
        else:
            return pdb_min


if __name__ == '__main__':
    import sys
    force_field = ForceFieldMinimizer(stiffness=10**8)
    # force_field(sys.argv[1], sys.argv[2], restrain_chain=[sys.argv[3]])
    force_field('data_tmp/docking/CB1R_0.pdb', 'data_tmp/docking/CB1R_0_openmm_relax_test.pdb',
                restrain_chain=['A'], restrain_backbone=True)
