from ..imports import load_dekkers
from ..emceeDeconv import getEmceeTrace


def test_import():
    """Test that we can import the data."""
    data = load_dekkers()


def test_trace():
    """Fit the model."""
    trace = getEmceeTrace()
