"""Configuration file for pytest tests."""
import pytest


@pytest.fixture
def chrome_options(chrome_options):
    """By default use these options with Chrome."""
    # https://pytest-selenium.readthedocs.io/en/latest/user_guide.html#id2
    chrome_options.add_argument("--headless")
    return chrome_options
