import unittest
from typing import List

from lazydock.pml.autodock_utils import DlgFile, merge_dlg


def sum_pose_num(dlg_lst: List[DlgFile]):
    return sum(map(lambda x: len(x.pose_lst), dlg_lst))


class TestMergeDlg(unittest.TestCase):

    def test_merge_0_dlg(self):
        # Test case with less than 2 dlg files
        result = merge_dlg([])
        self.assertEqual(result, None)

    def test_merge_2_dlg(self):
        # Test case with 2 dlg files
        dlg_lst = [DlgFile(path='data_tmp/dlg/1000run.dlg') for _ in range(2)]
        result = merge_dlg(dlg_lst)
        self.assertEqual(len(result.pose_lst), sum_pose_num(dlg_lst))

    def test_merge_4_dlg(self):
        # Test case with more than 2 dlg files
        dlg_lst = [DlgFile(path='data_tmp/dlg/1000run.dlg') for _ in range(4)]
        result = merge_dlg(dlg_lst)
        self.assertEqual(len(result.pose_lst), sum_pose_num(dlg_lst))


if __name__ == '__main__':
    unittest.main()