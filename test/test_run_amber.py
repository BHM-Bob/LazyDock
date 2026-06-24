import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from lazydock.scripts.run_amber import gamd_md


class TestGaMDRestartRequirements(unittest.TestCase):
    def setUp(self):
        os.environ.setdefault('AMBERHOME', '/tmp/amber')
        self.tempdir = tempfile.TemporaryDirectory()
        self.working_dir = Path(self.tempdir.name)
        self.prmtop = self.working_dir / 'system.prmtop'
        self.inpcrd = self.working_dir / 'step05_npt2.rst'
        self.prmtop.write_text('prmtop')
        self.inpcrd.write_text('rst')
        self.args = type('Args', (), {
            'batch_dir': [str(self.working_dir)],
            'prmtop_name': '*.prmtop',
            'inpcrd_name': '*.rst',
            'eq_nstlim': 3000000,
            'prod_nstlim': 100000000,
            'ntcmd': 1000000,
            'ieqdtd': 1000000,
            'nteb': None,
            'ntave': None,
            'ntcmdprep': None,
            'ntebprep': None,
            'igamd': 3,
            'sigma0P': 6.0,
            'sigma0D': 6.0,
            'ntp': 1,
            'barostat': 2,
            'cut': 10.0,
            'dt': 0.002,
            'prod_temp0': 300.0,
            'prod_ntpr': 50000,
            'prod_ntwx': 50000,
            'prod_ntwr': 5000000,
            'start_phase': 'equil',
            'end_phase': 'prod'
        })

    def tearDown(self):
        self.tempdir.cleanup()

    @patch('lazydock.scripts.run_amber.put_err')
    def test_equil_requires_velocity_restart(self, mock_put_err):
        cmd = gamd_md(self.args)
        with patch.object(cmd, '_check_velocities_with_pmemd', return_value=(False, 'could not find enough velocities')):
            result = cmd.run_phase('equil', self.prmtop, self.inpcrd, self.working_dir)
        self.assertFalse(result)
        mock_put_err.assert_called()

    @patch('lazydock.scripts.run_amber.put_err')
    def test_prod_requires_equil_restart(self, mock_put_err):
        cmd = gamd_md(self.args)
        with patch.object(cmd, '_check_velocities_with_pmemd', return_value=(True, '')):
            result = cmd.run_phase('prod', self.prmtop, self.inpcrd, self.working_dir)
        self.assertFalse(result)
        mock_put_err.assert_called()
