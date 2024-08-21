import unittest
import platform
import os
import tempfile
import shutil
import sys
import json
import warnings
import multiprocessing
import hashlib
import gzip


class UnitTest(unittest.TestCase):
    TEST_DATA = os.path.join(os.path.dirname(__file__), 'test_data')

    # This enables running tests from any directory.
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
    SRC_DIR = os.path.join(os.path.dirname(__file__), '..', 'src')

    def setUp(self):
        # Remember working directory before test case is launched.
        self.working_dir = os.getcwd()

        # Suppress import warnings (Python 3.10 on virtual Linux).
        warnings.simplefilter('ignore', category=ImportWarning)

        # Create and delete a temporary directory for each test case.
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Delete the temporary directory.
        shutil.rmtree(self.temp_dir)