import unittest

from lazydock.pyrt.pose_utils import _Pose


class TestPose(unittest.TestCase):

    def setUp(self):
        self.pose = _Pose(seq='A'*10)

    def test_add_resi_before(self):
        resi_name3 = "ASP"
        insert_pos = 1
        pos = 'before'
        self.pose.add_resi(resi_name3, insert_pos, pos)
        # Add assertions here

    def test_add_resi_after(self):
        resi_name3 = "GLU"
        insert_pos = 5
        pos = 'after'
        self.pose.add_resi(resi_name3, insert_pos, pos)
        # Add assertions here

    def test_delete_resi(self):
        resi_id = 3
        self.pose.delete_resi(resi_id)
        # Add assertions here

    def test_mutate_resi(self):
        resi_id = 2
        new_name1 = "C"
        pack_radius = 0.5
        self.pose.mutate_resi(resi_id, new_name1, pack_radius)
        # Add assertions here

    def test_add_chain(self):
        chain = _Pose(seq='A'*10)
        result = self.pose.add_chain(chain)
        # Add assertions here

    def test_del_region(self):
        start_resi_id = 1
        stop_resi_id = 3
        self.pose.del_region(start_resi_id, stop_resi_id)
        # Add assertions here

    def test_merge_chain(self):
        chain = _Pose(seq='C'*10)
        result = self.pose.merge_chain(chain)
        # Add assertions here

    def test_delete_chain(self):
        chain_id = 1
        self.pose.delete_chain(chain_id)
        # Add assertions here

    def test_swap_chain(self):
        order = ['A', 'B']
        self.pose.merge_chain(_Pose(seq='WWW'), True)
        print(self.pose.pose.pdb_info())
        # result = self.pose.swap_chain(order)
        # Add assertions here

    def test_split_chain(self):
        result = self.pose.split_chain()
        # Add assertions here

    def test_get_resi_id_in_pose_via_pdb(self):
        chain_name = "A"
        resi_id = 1
        result = self.pose.get_resi_id_in_pose_via_pdb(chain_name, resi_id)
        # Add assertions here

    def test_get_resi_id_in_pdb_via_pose(self):
        resi_id = 2
        result = self.pose.get_resi_id_in_pdb_via_pose(resi_id)
        # Add assertions here

    def test_get_chain_id(self):
        chain_name = "A"
        result = self.pose.get_chain_id(chain_name)
        # Add assertions here

    def test_get_chain_name(self):
        chain_id = 1
        result = self.pose.get_chain_name(chain_id)
        # Add assertions here

    def test_get_chain_seq(self):
        chain_id = 1
        result = self.pose.get_chain_seq(chain_id)
        # Add assertions here

if __name__ == '__main__':
    from pyrosetta import init
    init()
    
    unittest.main()