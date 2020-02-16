import os
import sys

# # PACKAGE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
#
# sys.path.append(PACKAGE_DIR)
TEST_DIR = os.path.dirname(__file__)

TEST_DATA_DIR = os.path.join(TEST_DIR, 'unit_test_data/')
EXAMPLES_DIR = TEST_DIR[:-5] + 'examples/'
