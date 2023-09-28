from typing import Optional

# Determine version using pyproject.toml file
try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:  # pragma: no cover
    from importlib_metadata import version, PackageNotFoundError


def get_package_version(package_name: str) -> Optional[str]:
    """Fetch the package version from the installed packages information.

    Args:
        package_name (str): Name of the package for which version is to be fetched.

    Returns:
        str: Version of the package. Returns "unknown" if package version is not found.
    """
    try:
        return version(package_name)
    except PackageNotFoundError:  # pragma: no cover
        return "unknown"


__version__ = get_package_version(__name__)
