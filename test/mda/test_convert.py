import unittest

import numpy as np
from MDAnalysis import Universe

from lazydock.gmx.mda.convert import FakeAtomGroup, PDBConverter


class TestFakeAtomGroupIndexing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.u = Universe(
            'data_tmp/MD/traj1/pull.tpr',
            'data_tmp/MD/traj1/pull.xtc'
        )
        cls.ag_protein = cls.u.select_atoms('protein')
        cls.fake = FakeAtomGroup(cls.ag_protein)
        cls.n_atoms = len(cls.ag_protein)

    def test_len(self):
        self.assertEqual(len(self.fake), self.n_atoms)
        self.assertEqual(len(self.fake), len(self.ag_protein))

    def test_single_integer_index_returns_length_one(self):
        result = self.fake[0]
        self.assertIsInstance(result, FakeAtomGroup)
        self.assertEqual(len(result), 1)
        np.testing.assert_array_equal(result.ids, self.fake.ids[[0]])
        np.testing.assert_array_equal(result.positions, self.fake.positions[[0]])

    def test_single_integer_negative_index(self):
        result = self.fake[-1]
        self.assertIsInstance(result, FakeAtomGroup)
        self.assertEqual(len(result), 1)
        np.testing.assert_array_equal(result.ids, self.fake.ids[[-1]])
        np.testing.assert_array_equal(result.resnames, self.fake.resnames[[-1]])

    def test_bool_array_mask(self):
        mask = np.zeros(self.n_atoms, dtype=bool)
        mask[0] = True
        mask[5] = True
        mask[10] = True
        result = self.fake[mask]
        self.assertIsInstance(result, FakeAtomGroup)
        self.assertEqual(len(result), 3)
        expected_positions = self.fake.positions[mask]
        np.testing.assert_array_equal(result.positions, expected_positions)
        np.testing.assert_array_equal(result.ids, self.fake.ids[mask])
        np.testing.assert_array_equal(result.resnames, self.fake.resnames[mask])
        np.testing.assert_array_equal(result.chainIDs, self.fake.chainIDs[mask])
        np.testing.assert_array_equal(result.resids, self.fake.resids[mask])
        np.testing.assert_array_equal(result.names, self.fake.names[mask])
        np.testing.assert_array_equal(result.elements, self.fake.elements[mask])
        np.testing.assert_array_equal(result.altLocs, self.fake.altLocs[mask])
        np.testing.assert_array_equal(result.occupancies, self.fake.occupancies[mask])
        np.testing.assert_array_equal(result.tempfactors, self.fake.tempfactors[mask])
        np.testing.assert_array_equal(result.segids, self.fake.segids[mask])
        np.testing.assert_array_equal(result.charges, self.fake.charges[mask])
        np.testing.assert_array_equal(result.icodes, self.fake.icodes[mask])

    def test_bool_array_all_true(self):
        mask = np.ones(self.n_atoms, dtype=bool)
        result = self.fake[mask]
        self.assertEqual(len(result), self.n_atoms)
        np.testing.assert_array_equal(result.ids, self.fake.ids)
        np.testing.assert_array_equal(result.positions, self.fake.positions)

    def test_bool_array_all_false(self):
        mask = np.zeros(self.n_atoms, dtype=bool)
        result = self.fake[mask]
        self.assertEqual(len(result), 0)
        self.assertEqual(result.positions.shape, (0, 3))
        self.assertEqual(len(result.ids), 0)

    def test_integer_array_index(self):
        idx = np.array([0, 3, 7, 2])
        result = self.fake[idx]
        self.assertIsInstance(result, FakeAtomGroup)
        self.assertEqual(len(result), 4)
        np.testing.assert_array_equal(result.positions, self.fake.positions[idx])
        np.testing.assert_array_equal(result.ids, self.fake.ids[idx])
        np.testing.assert_array_equal(result.resids, self.fake.resids[idx])

    def test_integer_list_index(self):
        idx = [1, 5, 9]
        result = self.fake[idx]
        self.assertIsInstance(result, FakeAtomGroup)
        self.assertEqual(len(result), 3)
        np.testing.assert_array_equal(result.positions, self.fake.positions[idx])
        np.testing.assert_array_equal(result.ids, self.fake.ids[idx])

    def test_slice_index(self):
        result = self.fake[10:20]
        self.assertIsInstance(result, FakeAtomGroup)
        self.assertEqual(len(result), 10)
        np.testing.assert_array_equal(result.positions, self.fake.positions[10:20])
        np.testing.assert_array_equal(result.ids, self.fake.ids[10:20])
        np.testing.assert_array_equal(result.resnames, self.fake.resnames[10:20])

    def test_slice_step(self):
        result = self.fake[::100]
        expected_len = (self.n_atoms + 99) // 100
        self.assertEqual(len(result), expected_len)
        np.testing.assert_array_equal(result.positions, self.fake.positions[::100])

    def test_chained_indexing(self):
        first = self.fake[10:100]
        second = first[np.array([0, 5, 10], dtype=int)]
        expected = self.fake[[10, 15, 20]]
        self.assertEqual(len(second), 3)
        np.testing.assert_array_equal(second.ids, expected.ids)
        np.testing.assert_array_equal(second.positions, expected.positions)

    def test_indexed_result_works_with_pdb_converter(self):
        mask = np.zeros(self.n_atoms, dtype=bool)
        mask[:50] = True
        result = self.fake[mask]
        converter = PDBConverter(result)
        pdb_str = converter.fast_convert()
        self.assertIsInstance(pdb_str, str)
        self.assertTrue(pdb_str.startswith('ATOM') or pdb_str.startswith('HETATM'))
        n_lines = len([l for l in pdb_str.strip().split('\n')
                       if l.startswith('ATOM') or l.startswith('HETATM')])
        self.assertEqual(n_lines, 50)

    def test_ids_not_reindexed_after_bool_mask(self):
        mask = np.zeros(self.n_atoms, dtype=bool)
        mask[5] = True
        mask[17] = True
        result = self.fake[mask]
        np.testing.assert_array_equal(result.ids, self.fake.ids[[5, 17]])
        self.assertNotEqual(list(result.ids), [1, 2])

    def test_single_integer_ids_not_reindexed(self):
        k = 42
        if k < self.n_atoms:
            result = self.fake[k]
            self.assertEqual(int(result.ids[0]), int(self.fake.ids[k]))


if __name__ == '__main__':
    unittest.main()
