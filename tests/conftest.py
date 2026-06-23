"""
holds command line options for pytest

"""


def pytest_addoption(parser):
    parser.addoption("--gpu-available", action="store_true", dest="gpu_available")
    parser.addoption("--nvidia", action="store_true", dest="nvidia")
    parser.addoption("--threads", action="store", default=1, dest="threads")
    parser.addoption("--euks", action="store_true", dest="euks")
