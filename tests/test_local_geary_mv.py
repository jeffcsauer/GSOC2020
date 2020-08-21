import unittest
import libpysal
from libpysal.common import pandas, RTOL, ATOL
from esda import local_geary
import numpy as np

PANDAS_EXTINCT = pandas is None

from ..local_geary_mv import Local_Geary_MV

class Local_Geary_MV_Tester(unittest.TestCase):
    def setUp(self):
        np.random.seed(10)
        self.w = libpysal.io.open(libpysal.examples.get_path("stl.gal")).read()
        f = libpysal.io.open(libpysal.examples.get_path("stl_hom.txt"))
        self.y1 = np.array(f.by_col['HR8893'])
        self.y2 = np.array(f.by_col['HC8488'])

    def test_local_geary_mv(self):
        lG_mv = Local_Geary_MV(connectivity=self.w).fit([self.y1, self.y2])
        self.assertAlmostEqual(lG_mv.localG[0], 0.4096931479581422)
        self.assertAlmostEqual(lG_mv.p_sim[0], 0.32)
        
suite = unittest.TestSuite()
test_classes = [
    Local_Geary_MV_Tester
]
for i in test_classes:
    a = unittest.TestLoader().loadTestsFromTestCase(i)
    suite.addTest(a)

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite)