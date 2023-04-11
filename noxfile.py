import os
import pathlib
import shutil
from itertools import chain

import nox
import toml

PROJECT = "permamodel"
HERE = pathlib.Path(__file__)
ROOT = HERE.parent
PATHS = [PROJECT + "/components", PROJECT + "/tests", HERE.name]
PYTHON_VERSIONS = ["3.9", "3.10", "3.11"]
BUILD_DIR = ROOT / "build"


@nox.session(python=PYTHON_VERSIONS)
def test(session: nox.Session) -> None:
    """Run the tests."""
    session.install(".[testing]")

    args = [
        "--cov",
        PROJECT,
        "-vvv",
    ] + session.posargs

    if "CI" in os.environ:
        args.append(f"--cov-report=xml:{ROOT.absolute()!s}/coverage.xml")
    session.run("pytest", *args)

    if "CI" not in os.environ:
        session.run("coverage", "report", "--ignore-errors", "--show-missing")


@nox.session(name="test-bmi", python=PYTHON_VERSIONS, venv_backend="conda")
def test_bmi(session: nox.Session) -> None:
    """Test the Basic Model Interface."""
    session.conda_install("bmi-tester", "pymt>=1.3")
    session.install(".")

    bmi_test_dir = BUILD_DIR / "bmi_test"
    with bmi_test_setup(bmi_test_dir):
        session.run(
            "bmi-test",
            "permamodel.components.bmi_frost_number:BmiFrostnumberMethod",
            "--config-file",
            "./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg",
            "--root-dir",
            bmi_test_dir,
            "-vvv",
        )
    with bmi_test_setup(bmi_test_dir):
        session.run(
            "bmi-test",
            "permamodel.components.bmi_Ku_component:BmiKuMethod",
            "--config-file",
            "./permamodel/examples/Ku_method.cfg",
            "--root-dir",
            bmi_test_dir,
            "-vvv",
        )
    with bmi_test_setup(bmi_test_dir):
        cfg_file = "./permamodel/examples/Ku_bmi_example_config.toml"
        _set_absolute_path_in_config(cfg_file)
        session.run(
            "bmi-test",
            "permamodel.components.bmi_Ku:BmiKuModel",
            "--config-file",
            cfg_file,
            "--root-dir",
            bmi_test_dir,
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
    shutil.rmtree(BUILD_DIR, ignore_errors=True)
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


class bmi_test_setup:
    def __init__(self, dir="."):
        self._dir = dir

    def __enter__(self):
        os.makedirs(self._dir, exist_ok=True)

    def __exit__(self, type_, value, traceback):
        shutil.rmtree(self._dir)


def _set_absolute_path_in_config(cfg_file):
    with open(cfg_file, "r") as f:
        config = toml.load(f)

    if not pathlib.Path(config["directories"]["inputs_dir"]).is_absolute():
        config["directories"][
            "inputs_dir"
        ] = f'{ROOT / config["directories"]["inputs_dir"]}{os.sep}'
        config["directories"][
            "outputs_dir"
        ] = f'{ROOT / config["directories"]["outputs_dir"]}{os.sep}'

        with open(cfg_file, "w") as f:
            toml.dump(config, f)
