import unittest
from pathlib import Path

import numpy as np
import torch
from lazydock.algorithm.rms import (batch_calc_rmsd_rotational_matrix, batch_rmsd,
                                    calc_rms_rotational_matrix, fit_to, batch_fit_to,
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
        
        # Load first 100 frames
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
            rotation, mda_rmsd = rotation_matrix_mda(cls.ref_coords, mobile)
            cls.mda_rmsds.append(mda_rmsd)
            cls.mda_rotations.append(rotation.ravel())
        cls.mda_rmsds = np.array(cls.mda_rmsds)
        cls.mda_rotations = np.array(cls.mda_rotations)

    @staticmethod
    def _convert_to_backend(data, backend):
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

    def test_batch_rmsd_rotational_matrix(self):
        """Verify batch rmsd matches single frame results"""
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
                
                np.testing.assert_allclose(single_rmsds, self.mda_rmsds, atol=1e-5)
                np.testing.assert_allclose(batch_rmsds, self.mda_rmsds, atol=1e-5)
                np.testing.assert_allclose(batch_rmsds, single_rmsds, atol=1e-5)

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
                    
                    for mobile_t in self.mobile_coords:                        
                        # Compute RMSD
                        computed = rmsd(
                            self._convert_to_backend(ref, backend),
                            self._convert_to_backend(mobile_t, backend),
                            **case,
                            backend=backend
                        )
                        
                        # Convert to numpy
                        if backend != 'numpy':
                            computed = computed.cpu().numpy()
                        
                        # Compute reference with MDAnalysis
                        expected = rmsd_mda(ref, mobile_t, **case)
                        
                        np.testing.assert_allclose(computed, expected, atol=1e-5)
                    
    def test_batch_rmsd(self):
        """Verify batch rmsd matches single frame results"""
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
                batch_rmsds = batch_rmsd(batch_ref_t, batch_mobile_t, backend=backend)
                
                # Single frame computations via MDAnalysis with center=True and superposition=True
                mda_rmsds = []
                for mobile in self.mobile_coords:
                    mda_rmsds.append(
                        rmsd_mda(self.ref_coords, mobile, center=True, superposition=True)
                    )
                
                # Convert to numpy for comparison
                if backend != 'numpy':
                    batch_rmsds = batch_rmsds.cpu().numpy()
                
                np.testing.assert_allclose(batch_rmsds, mda_rmsds, atol=1e-5)


class TestFitAPIs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Load test trajectory data once for all tests"""
        # Modify paths according to your test data location
        paths = opts_file(Path(__file__).parent.parent.parent / 'data_tmp/test_config.json', way='json')['test_rms']
        cls.u = Universe(paths['top_path'], paths['traj_path'])
        cls.ag = cls.u.select_atoms('protein and name CA')
        cls.ref_com = cls.ag.center_of_geometry()
        cls.ref_coords = cls.ag.positions.copy().astype(np.float64)

        # Load first 100 frames and precompute MDAnalysis results for comparison
        cls.mda_rmsds = []
        cls.mda_positions = []
        ref_com = cls.ag.center_of_geometry()
        for _ in cls.u.trajectory[:100]:
            new_pos, mds_rmsd = fit_to_mda(cls.ag.positions, cls.ref_coords,
                                           cls.ag, cls.ag.center_of_geometry(), ref_com)
            cls.mda_rmsds.append(mds_rmsd)
            cls.mda_positions.append(new_pos.positions.copy().astype(np.float64))
        cls.mda_rmsds = np.array(cls.mda_rmsds)
        cls.mda_positions = np.array(cls.mda_positions)
    
    def test_fit_to(self):
        """Test fit_to function with different backends"""
        backends = ['numpy', 'torch', 'cuda'] if torch.cuda.is_available() else ['numpy', 'torch']

        for backend in backends:
            with self.subTest(backend=backend):
                if backend == 'cuda' and not torch.cuda.is_available():
                    self.skipTest("CUDA not available")
                
                ref_coords = TestRMSDAPIs._convert_to_backend(self.ref_coords, backend)
                ref_com = TestRMSDAPIs._convert_to_backend(self.ref_com, backend)
                for i, _ in enumerate(self.u.trajectory[:100]):
                    # Load frame
                    mobile_t = self.ag.positions.copy().astype(np.float64)
                    mobile_com = self.ag.center_of_geometry()
                    
                    mobile_t = TestRMSDAPIs._convert_to_backend(mobile_t, backend)
                    mobile_com = TestRMSDAPIs._convert_to_backend(mobile_com, backend)
                    
                    # Compute with current backend
                    new_pos, np_rmsd = fit_to(mobile_t, ref_coords, mobile_com, ref_com, backend=backend)
                    
                    if backend!= 'numpy':
                        new_pos = new_pos.cpu().numpy()
                        np_rmsd = np_rmsd.cpu().numpy()
                    
                    np.testing.assert_allclose(self.mda_rmsds[i], np_rmsd, atol=1e-5)
                    np.testing.assert_allclose(self.mda_positions[i], new_pos, atol=1e-5)

    def test_batch_fit_to(self):
        """Verify batch fit_to matches single frame results"""
        backends = ['numpy', 'torch', 'cuda'] if torch.cuda.is_available() else ['numpy', 'torch']
        batch_ref = np.array([self.ref_coords]*len(self.mda_rmsds))

        for backend in backends:
            with self.subTest(backend=backend):
                if backend == 'cuda' and not torch.cuda.is_available():
                    self.skipTest("CUDA not available")
                
                # gether coords
                batch_mobile = []
                for i, _ in enumerate(self.u.trajectory[:100]):
                    batch_mobile.append(self.ag.positions.copy().astype(np.float64))
                batch_mobile = np.array(batch_mobile)
                
                # Batch computation
                batch_ref_t = TestRMSDAPIs._convert_to_backend(batch_ref, backend)
                batch_mobile_t = TestRMSDAPIs._convert_to_backend(batch_mobile, backend)
                batch_positions, batch_rmsds = batch_fit_to(batch_mobile_t, batch_ref_t, backend=backend)
                
                # Convert to numpy for comparison
                if backend!= 'numpy':
                    batch_positions = batch_positions.cpu().numpy()
                    batch_rmsds = batch_rmsds.cpu().numpy()

                # Compare with MDAnalysis
                np.testing.assert_allclose(batch_positions, self.mda_positions, atol=1e-5)
                np.testing.assert_allclose(batch_rmsds, self.mda_rmsds, atol=1e-5)


if __name__ == '__main__':
    unittest.main()