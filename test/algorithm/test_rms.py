import unittest
from pathlib import Path

import numpy as np
import torch
from lazydock.algorithm.rms import (batch_calc_rmsd_rotational_matrix, batch_fast_calc_rmsd,
                                    batch_inner_product, batch_rmsd,
                                    calc_rms_rotational_matrix,
                                    fast_calc_rmsd_rotation, fit_to,
                                    inner_product, rmsd)
from mbapy_lite.file import opts_file
from MDAnalysis import Universe
from MDAnalysis.analysis.align import _fit_to as fit_to_mda
from MDAnalysis.analysis.align import rotation_matrix as rotation_matrix_mda
from MDAnalysis.analysis.rms import rmsd as rmsd_mda


class TestRMSDAPIs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Load test trajectory data once for all tests"""
        # Modify paths according to your test data location
        paths = opts_file(Path(__file__).parent.parent.parent / 'data_tmp/test_config.json', way='json')['test_rms']
        cls.u = Universe(paths['top_path'], paths['traj_path'])
        cls.ag = cls.u.select_atoms('protein and name CA')
        cls.ref_coords = []
        cls.mobile_coords = []
        
        # Load first 10 frames
        for ts in cls.u.trajectory[:100]:
            coords = cls.ag.positions.copy().astype(np.float64)
            if isinstance(cls.ref_coords, list):
                cls.ref_coords = coords  # First frame as reference
            else:
                cls.mobile_coords.append(coords)
        
        # Precompute MDAnalysis results for comparison
        cls.mda_rmsds = []
        cls.mda_rotations = []
        for mobile in cls.mobile_coords:
            _, mda_rmsd = rotation_matrix_mda(cls.ref_coords, mobile)
            cls.mda_rmsds.append(mda_rmsd)
            rotation = rotation_matrix_mda(cls.ref_coords, mobile)[0]
            cls.mda_rotations.append(rotation.ravel())
        cls.mda_rmsds = np.array(cls.mda_rmsds)
        cls.mda_rotations = np.array(cls.mda_rotations)

    def _convert_to_backend(self, data, backend):
        """Convert numpy array to backend tensor"""
        if backend == 'numpy':
            return data.astype(np.float64)
        elif backend == 'torch':
            return torch.tensor(data, dtype=torch.float64)
        elif backend == 'cuda':
            return torch.tensor(data, dtype=torch.float64).cuda()

    def test_inner_product(self):
        """Test inner_product with different backends"""
        backends = ['numpy', 'torch', 'cuda'] if torch.cuda.is_available() else ['numpy', 'torch']
        ref = self.ref_coords
        mobile = self.mobile_coords[0]
        
        for backend in backends:
            with self.subTest(backend=backend):
                if backend == 'cuda' and not torch.cuda.is_available():
                    self.skipTest("CUDA not available")
                
                # Compute with current backend
                ref_t = self._convert_to_backend(ref, backend)
                mobile_t = self._convert_to_backend(mobile, backend)
                A, E0 = inner_product(ref_t, mobile_t, backend=backend)
                
                # Convert to numpy for comparison
                if backend != 'numpy':
                    A = A.cpu().numpy() if hasattr(A, 'cpu') else A
                    E0 = E0.cpu().numpy() if hasattr(E0, 'cpu') else E0
                
                # Compare with numpy implementation
                A_np, E0_np = inner_product(ref, mobile, backend='numpy')
                np.testing.assert_allclose(A, A_np, atol=1e-6)
                np.testing.assert_allclose(E0, E0_np, atol=1e-6)

    def test_rotation_calculation(self):
        """Test rotation matrix calculation against MDAnalysis"""
        backends = ['numpy', 'torch', 'cuda'] if torch.cuda.is_available() else ['numpy', 'torch']
        
        for idx, mobile in enumerate(self.mobile_coords):
            for backend in backends:
                with self.subTest(frame=idx, backend=backend):
                    if backend == 'cuda' and not torch.cuda.is_available():
                        self.skipTest("CUDA not available")
                    
                    # Compute rotation matrix
                    ref_t = self._convert_to_backend(self.ref_coords, backend)
                    mobile_t = self._convert_to_backend(mobile, backend)
                    rot = np.zeros(9) if backend == 'numpy' else torch.zeros(9)
                    
                    rmsd_val, rot = calc_rms_rotational_matrix(
                        ref_t, mobile_t, rot=rot, backend=backend
                    )
                    
                    # Convert results to numpy
                    if backend != 'numpy':
                        rmsd_val = rmsd_val.cpu().numpy()
                        rot = rot.cpu().numpy()
                    
                    # Verify against MDAnalysis
                    np.testing.assert_allclose(rmsd_val, self.mda_rmsds[idx], atol=1e-5)
                    np.testing.assert_allclose(rot, self.mda_rotations[idx], atol=1e-5)

    def test_batch_processing(self):
        """Verify batch processing matches single frame results"""
        backends = ['numpy', 'torch', 'cuda'] if torch.cuda.is_available() else ['numpy', 'torch']
        batch_ref = np.array([self.ref_coords]*len(self.mobile_coords))
        batch_mobile = np.array(self.mobile_coords)
        
        for backend in backends:
            with self.subTest(backend=backend):
                if backend == 'cuda' and not torch.cuda.is_available():
                    self.skipTest("CUDA not available")
                
                # Batch computation
                batch_ref_t = self._convert_to_backend(batch_ref, backend)
                batch_mobile_t = self._convert_to_backend(batch_mobile, backend)
                batch_rmsds = batch_calc_rmsd_rotational_matrix(batch_ref_t, batch_mobile_t, backend=backend)
                
                # Single frame computations
                single_rmsds = []
                for mobile in self.mobile_coords:
                    ref_t = self._convert_to_backend(self.ref_coords, backend)
                    mobile_t = self._convert_to_backend(mobile, backend)
                    single_rmsds.append(
                        calc_rms_rotational_matrix(ref_t, mobile_t, backend=backend)
                    )
                
                # Convert to numpy for comparison
                if backend != 'numpy':
                    batch_rmsds = batch_rmsds.cpu().numpy()
                    single_rmsds = [r.cpu().numpy() for r in single_rmsds]
                
                err_tol = 1e-5 if backend == 'numpy' else 0.005  # Different tolerance for 'cuda' and 'torch'
                np.testing.assert_allclose(single_rmsds, self.mda_rmsds, atol=err_tol)
                np.testing.assert_allclose(batch_rmsds, self.mda_rmsds, atol=err_tol)
                np.testing.assert_allclose(batch_rmsds, single_rmsds, atol=err_tol)

    def test_rmsd_function(self):
        """Test rmsd function with different parameters"""
        backends = ['numpy', 'torch', 'cuda'] if torch.cuda.is_available() else ['numpy', 'torch']
        test_cases = [
            {'center': False, 'superposition': False},
            {'center': True, 'superposition': False},
            {'center': True, 'superposition': True}
        ]
        
        for case in test_cases:
            for backend in backends:
                with self.subTest(backend=backend, params=case):
                    if backend == 'cuda' and not torch.cuda.is_available():
                        self.skipTest("CUDA not available")
                    
                    ref = self.ref_coords
                    mobile = self.mobile_coords[0]
                    
                    # Compute RMSD
                    computed = rmsd(
                        self._convert_to_backend(ref, backend),
                        self._convert_to_backend(mobile, backend),
                        **case,
                        backend=backend
                    )
                    
                    # Convert to numpy
                    if backend != 'numpy':
                        computed = computed.cpu().numpy()
                    
                    # Compute reference with MDAnalysis
                    expected = rmsd_mda(ref, mobile, **case)
                    
                    np.testing.assert_allclose(computed, expected, atol=1e-5)

if __name__ == '__main__':
    unittest.main()