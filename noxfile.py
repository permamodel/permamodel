import os
import pathlib
import shutil
from itertools import chain

import nox

PROJECT = "permamodel"
HERE = pathlib.Path(__file__)
ROOT = HERE.parent
PATHS = [PROJECT + "/components", PROJECT + "/tests", HERE.name]
PYTHON_VERSIONS = ["3.9", "3.10", "3.11"]


@nox.session(python=PYTHON_VERSIONS)
def test(session: nox.Session) -> None:
    """Run the tests."""
    session.install(".[testing]")
    args = session.posargs or ["--cov", "--cov-report=term", "-vvv"]
    session.run("pytest", *args)


@nox.session(name="test-bmi", python=PYTHON_VERSIONS, venv_backend="conda")
def test_bmi(session: nox.Session) -> None:
    """Test the Basic Model Interface."""
    session.conda_install("bmi-tester", "pymt>=1.3")
    session.install(".")
    fn_tmpdir = session.create_tmp()
    session.run(
        "bmi-test",
        "permamodel.components.bmi_frost_number:BmiFrostnumberMethod",
        "--config-file",
        "./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg",
        "--root-dir",
        fn_tmpdir,
        "-vvv",
    )
    ku_tmpdir = session.create_tmp()
    session.run(
        "bmi-test",
        "permamodel.components.bmi_ku_component:BmiKuMethod",
        "--config-file",
        "./permamodel/examples/Ku_method.cfg",
        "--root-dir",
        ku_tmpdir,
        "-vvv",
    )


@nox.session
def format(session: nox.Session) -> None:
    """Clean lint and assert style."""
    session.install(".[dev]")

    if session.posargs:
        black_args = session.posargs
    else:
        black_args = []

    session.run("black", *black_args, *PATHS)
    session.run("isort", *PATHS)
    session.run("flake8", *PATHS)


@nox.session
def build(session: nox.Session) -> None:
    """Build source and binary distributions."""
    session.install(".[build]")
    session.run("python", "-m", "build")


@nox.session
def release(session):
    """Tag, build, and publish a new release to PyPI."""
    session.install(".[build]")
    session.run("fullrelease")


@nox.session(name="testpypi")
def publish_testpypi(session):
    """Upload package to TestPyPI."""
    build(session)
    session.run("twine", "check", "dist/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "--repository-url",
        "https://test.pypi.org/legacy/",
        "dist/*",
    )


@nox.session(name="pypi")
def publish_pypi(session):
    """Upload package to PyPI."""
    build(session)
    session.run("twine", "check", "dist/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "--repository-url",
        "https://upload.pypi.org/legacy/",
        "dist/*",
    )


@nox.session(python=False)
def clean(session):
    """Remove virtual environments, build files, and caches."""
    shutil.rmtree("build", ignore_errors=True)
    shutil.rmtree("dist", ignore_errors=True)
    shutil.rmtree(f"{PROJECT}.egg-info", ignore_errors=True)
    shutil.rmtree(".pytest_cache", ignore_errors=True)
    shutil.rmtree(".venv", ignore_errors=True)
    if os.path.exists(".coverage"):
        os.remove(".coverage")
    for p in chain(ROOT.rglob("*.py[co]"), ROOT.rglob("__pycache__")):
        if p.is_dir():
            p.rmdir()
        else:
            p.unlink()


@nox.session(python=False)
def nuke(session):
    """Clean and also remove the .nox directory."""
    clean(session)
    shutil.rmtree(".nox", ignore_errors=True)
