import pytest

@pytest.fixture(scope="session", autouse=True)
def set_up_overall(request):
    request.addfinalizer(tear_down)

def tear_down():
    import os
    import glob

    patterns = ['*pdb', '*dat']
    pytest_scratches = []
    for pat in patterns:
        pytest_scratches.extend(glob.glob(pat))
    for fl in pytest_scratches:
        os.unlink(fl)
