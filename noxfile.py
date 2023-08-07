"""Nox sessions."""
import tempfile
from typing import Any

import nox
from nox.sessions import Session


package = "covsonar"
nox.options.sessions = "lint", "tests"
locations = "src", "tests", "./noxfile.py"


def install_with_constraints(session: Session, *args: str, **kwargs: Any) -> None:
    """Install packages constrained by Poetry's lock file.
    This function is a wrapper for nox.sessions.Session.install. It
    invokes pip to install packages inside of the session's virtualenv.
    Additionally, pip is passed a constraints file generated from
    Poetry's lock file, to ensure that the packages are pinned to the
    versions specified in poetry.lock. This allows you to manage the
    packages as Poetry development dependencies.
    Arguments:
        session: The Session object.
        args: Command-line arguments for pip.
        kwargs: Additional keyword arguments for Session.install.
    """
    with tempfile.NamedTemporaryFile() as requirements:
        session.run(
            "poetry",
            "export",
            "--with",
            "dev",
            "--format=requirements.txt",
            "--without-hashes",
            f"--output={requirements.name}",
            external=True,
        )
        session.install(f"--requirement={requirements.name}", *args, **kwargs)


@nox.session(python="3.9")
def black(session: Session) -> None:
    """Run black code formatter."""
    args = session.posargs or locations
    install_with_constraints(session, "black")
    session.run("black", *args)


@nox.session(python="3.9")
def zimports(session: Session) -> None:
    """Run zimports import formatter."""
    args = session.posargs or locations
    install_with_constraints(session, "zimports")
    session.run("zimports", "-m", "covsonar,tests", "-W", "2", *args)


@nox.session(python="3.9")
def lint(session: Session) -> None:
    """Lint using flake8."""
    args = session.posargs or locations
    install_with_constraints(
        session,
        "flake8",
        "flake8-annotations",
        "flake8-bandit",
        "flake8-black",
        "flake8-bugbear",
        "flake8-docstrings",
        "flake8-import-order",
        "darglint",
    )
    session.run("flake8", *args)


@nox.session(python="3.9", venv_backend="mamba")
def tests(session):
    args = session.posargs or ["--cov"]
    session.conda_install("parasail-python==1.3.4", "libiconv", channel="bioconda")
    session.run("poetry", "install", "--only", "main", external=True)
    install_with_constraints(
        session, "coverage[toml]", "nox", "pytest", "pytest-cov", "pytest-sugar"
    )
    session.run("pytest", "--doctest-modules", *args)
