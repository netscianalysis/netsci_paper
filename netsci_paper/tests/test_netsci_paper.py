"""
Unit and regression test for the netsci_paper package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import netsci_paper


def test_netsci_paper_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "netsci_paper" in sys.modules
