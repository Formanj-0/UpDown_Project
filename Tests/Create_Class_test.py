#import inc_dec    # The code to test
import unittest   # The test framework
import numpy as np

class Test_TestIncrementDecrement(unittest.TestCase):
    def test_load_data(self):
        self.assertEqual(True, np.arange(12).reshape((3, 4)))

    def test_decrement(self):
        self.assertEqual(inc_dec.decrement(3), 4)

if __name__ == '__main__':
    unittest.main()