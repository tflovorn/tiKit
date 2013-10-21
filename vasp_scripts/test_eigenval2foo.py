import unittest
import eigenval2foo

class FindSymmetryPoints(unittest.TestCase):
    def setUp(self):
        self.e = eigenval2foo.EIGENVAL("TEST_EIGENVAL")

    def test_known_points(self):
        expectedIndices = [0, 24, 48, 72, 96]
        expectedPoints = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5), (0.5, 0.5, 0.0),
                    (0.0, 0.0, 0.0), (0.5, 0.0, 0.0)]
        self.assertEqual(expectedIndices, self.e.symIndices)
        for i in range(len(expectedIndices)):
            kExpected = expectedPoints[i]
            kVal = self.e.kpoints[expectedIndices[i]]
            self.assertEqual(kExpected, kVal)

if __name__ == "__main__":
    unittest.main()
