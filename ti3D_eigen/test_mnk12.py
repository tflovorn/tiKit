import unittest
import numpy as np
import mnk12

class TestBlockTridiagonal(unittest.TestCase):
    # H(k) may include top-left and bottom-right 4x4 blocks if periodic but
    # otherwise should be tridiagonal in 4x4 blocks.
    # Assume non-periodic for now.
    def test_btd(self):
        # arbitrary k (may want to pick several)
        k = np.array([0.05, 0.1, 0.03])
        # generate some H(k)'s
        H_list = []
        for numLayers in range(3, 6):
            H_fn = mnk12.HamiltonianFn("mnk12", numLayers)
            H_list.append(H_fn(k))
        # call _check_btd on each
        for H in H_list:
            self._check_btd(H)

    def _check_btd(self, H):
        for i in range(0, len(H), 4):
            for j in range(0, len(H), 4):
                if j != i and j != i-4 and j != i+4:
                    self._check_block_zero(H, i, j)

    # Fails if the 4x4 block with upper-left corner (i, j) contains
    # nonzero elements.
    def _check_block_zero(self, H, i, j):
        for r in range(4):
            for s in range(4):
                self.assertEqual(H[i+r, j+s], 0)

if __name__ == "__main__":
    unittest.main()
