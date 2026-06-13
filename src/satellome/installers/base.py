"""
Base utilities for installers
"""

import os
import json
import shutil
import hashlib
import platform
import logging
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)

# Bumped when the provenance-manifest schema changes so old manifests can be
# detected and migrated rather than silently misread.
MANIFEST_VERSION = 1


def detect_platform() -> Tuple[str, str]:
    """
    Detect the current platform and architecture.

    Returns:
        Tuple[str, str]: (platform_name, architecture)
        platform_name: 'linux', 'darwin' (macOS), or 'unknown'
        architecture: 'x86_64', 'arm64', or 'unknown'
    """
    system = platform.system().lower()
    machine = platform.machine().lower()

    # Normalize platform names
    if system == 'linux':
        platform_name = 'linux'
    elif system == 'darwin':
        platform_name = 'darwin'
    else:
        platform_name = 'unknown'

    # Normalize architecture
    if machine in ['x86_64', 'amd64']:
        arch = 'x86_64'
    elif machine in ['arm64', 'aarch64']:
        arch = 'arm64'
    else:
        arch = 'unknown'

    return platform_name, arch


def get_satellome_bin_dir() -> Path:
    """
    Get or create the Satellome binary directory.

    Priority:
    1. <site-packages>/satellome/bin/ (primary location, cleaner)
    2. ~/.satellome/bin/ (fallback if no write permissions)

    Returns:
        Path: Path to binary directory
    """
    # Try to use package directory first (cleaner, no pollution of user's home)
    try:
        import satellome
        package_dir = Path(satellome.__file__).parent
        bin_dir = package_dir / 'bin'

        # Test if we can write to this directory
        bin_dir.mkdir(parents=True, exist_ok=True)
        test_file = bin_dir / '.write_test'
        try:
            test_file.touch()
            test_file.unlink()
            logger.debug(f"Using package binary directory: {bin_dir}")
            return bin_dir
        except (PermissionError, OSError):
            logger.debug(f"No write permission to {bin_dir}, falling back to ~/.satellome/bin/")
    except Exception as e:
        logger.debug(f"Could not use package directory: {e}")

    # Fallback to user home directory
    bin_dir = Path.home() / '.satellome' / 'bin'
    bin_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Using home directory: {bin_dir}")
    return bin_dir


def check_command_exists(command: str) -> bool:
    """
    Check if a command is available in the system PATH.

    Args:
        command: Command name to check

    Returns:
        bool: True if command exists, False otherwise
    """
    return shutil.which(command) is not None


def check_binary_exists(binary_name: str, check_system_path: bool = True) -> Optional[str]:
    """
    Check if a binary exists in Satellome bin directory or system PATH.

    Args:
        binary_name: Name of the binary to check
        check_system_path: Whether to check system PATH in addition to ~/.satellome/bin/

    Returns:
        Optional[str]: Full path to binary if found, None otherwise
    """
    # Check ~/.satellome/bin/ first (higher priority)
    satellome_bin = get_satellome_bin_dir() / binary_name
    if satellome_bin.exists() and os.access(satellome_bin, os.X_OK):
        return str(satellome_bin)

    # Check system PATH if requested
    if check_system_path:
        system_path = shutil.which(binary_name)
        if system_path:
            return system_path

    return None


def verify_installation(binary_name: str, test_command: Optional[str] = None) -> bool:
    """
    Verify that a binary is properly installed and executable.

    Args:
        binary_name: Name of the binary to verify
        test_command: Optional command to test (e.g., "fastan --help")

    Returns:
        bool: True if verification successful, False otherwise
    """
    binary_path = check_binary_exists(binary_name)

    if not binary_path:
        logger.error(f"{binary_name} not found in PATH or ~/.satellome/bin/")
        return False

    logger.info(f"Found {binary_name} at: {binary_path}")

    # If test command provided, try to run it
    if test_command:
        import subprocess
        try:
            subprocess.run(
                test_command.split(),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=5,
                check=False  # Some programs return non-zero for --help
            )
            logger.info(f"{binary_name} test command executed successfully")
            return True
        except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError) as e:
            logger.error(f"Failed to run test command for {binary_name}: {e}")
            return False

    return True


def find_system_gcc() -> Optional[str]:
    """
    Find system gcc (not conda's gcc which may have linking issues).

    Returns:
        Optional[str]: Path to system gcc, or None if not found
    """
    for gcc_path in ['/usr/bin/gcc', '/usr/local/bin/gcc']:
        if os.path.exists(gcc_path):
            return gcc_path
    return None


def run_make_with_fallback(
    cwd: Path,
    timeout: int = 300,
    clean_command: Optional[list] = None
) -> Tuple[bool, str]:
    """
    Run make with fallback to system gcc if conda gcc fails with -lz error.

    Args:
        cwd: Working directory for make
        timeout: Timeout in seconds
        clean_command: Command to clean failed build (default: ['make', 'clean'])

    Returns:
        Tuple[bool, str]: (success, error_message)
    """
    import subprocess

    if clean_command is None:
        clean_command = ['make', 'clean']

    compile_env = os.environ.copy()

    # First attempt: use default compiler
    result = subprocess.run(
        ['make'],
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        env=compile_env
    )

    if result.returncode == 0:
        return True, ""

    stderr_output = result.stderr.decode()

    # Check if error is related to conda gcc and -lz linking issue
    if 'cannot find -lz' in stderr_output and 'conda' in os.environ.get('PATH', ''):
        logger.warning("Compilation failed with conda gcc (cannot find -lz)")
        logger.info("Retrying with system gcc...")

        # Clean up failed build
        subprocess.run(clean_command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Try with system gcc
        system_gcc = find_system_gcc()
        if system_gcc:
            compile_env['CC'] = system_gcc
            logger.info(f"Using system compiler: {system_gcc}")

            result = subprocess.run(
                ['make'],
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
                env=compile_env
            )

            if result.returncode == 0:
                return True, ""

            return False, f"Failed to compile with system gcc:\n{result.stderr.decode()}"

        return False, f"System gcc not found. Original error:\n{stderr_output}"

    return False, f"Compilation failed:\n{stderr_output}"


def check_build_dependencies() -> Tuple[bool, str]:
    """
    Check if required build tools are available.

    Returns:
        Tuple[bool, str]: (success, error_message)
    """
    missing = []

    # Check for git
    if not check_command_exists('git'):
        missing.append('git')

    # Check for make
    if not check_command_exists('make'):
        missing.append('make')

    # Check for C compiler (gcc or clang)
    has_compiler = check_command_exists('gcc') or check_command_exists('clang') or check_command_exists('cc')
    if not has_compiler:
        missing.append('gcc or clang')

    if missing:
        error_msg = f"Missing required build tools: {', '.join(missing)}\n"
        system_name, _ = detect_platform()

        if system_name == 'darwin':
            error_msg += "On macOS, install Xcode Command Line Tools: xcode-select --install"
        elif system_name == 'linux':
            error_msg += "On Ubuntu/Debian: sudo apt-get install build-essential git\n"
            error_msg += "On CentOS/RHEL: sudo yum groupinstall 'Development Tools' && sudo yum install git"

        return False, error_msg

    return True, ""


# ---------------------------------------------------------------------------
# Binary provenance & integrity manifests
#
# Each managed binary gets a sibling "<binary>.manifest.json" recorded at
# install time (git commit of the source, SHA-256 of the built artifact, the
# self-reported build signature, platform). At runtime the pipeline recomputes
# the SHA-256 and compares it. A mismatch (corrupt / replaced / silently stale
# binary) is a user-visible error, never a silent fallback.
# ---------------------------------------------------------------------------


def manifest_path_for(binary_path) -> Path:
    """Return the manifest path that sits next to a binary."""
    binary_path = Path(binary_path)
    return binary_path.with_name(binary_path.name + ".manifest.json")


def compute_sha256(path) -> str:
    """Compute the SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_build_signature(binary_path) -> Optional[str]:
    """
    Best-effort read of a binary's self-reported build signature.

    Tries ``<binary> -V`` first (FasTAN fork prints its signature to stdout and
    exits 0). Falls back to scanning the binary image for the embedded
    ``FASTAN_SATELLOME_SIGNATURE=...`` marker so the build can still be
    identified when it cannot be executed (e.g. cross-arch).

    Returns the signature string, or None if neither mechanism yields one.
    """
    binary_path = str(binary_path)

    try:
        result = subprocess.run(
            [binary_path, "-V"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5,
            check=False,
        )
        out = result.stdout.decode(errors="replace").strip()
        if result.returncode == 0 and ("satellome-fork" in out or out.startswith("FasTAN")):
            return out
    except (OSError, subprocess.SubprocessError) as e:
        logger.debug(f"Could not run '{binary_path} -V': {e}")

    # Fall back to the embedded greppable marker.
    marker = b"FASTAN_SATELLOME_SIGNATURE="
    try:
        with open(binary_path, "rb") as f:
            data = f.read()
        idx = data.find(marker)
        if idx != -1:
            end = idx
            while end < len(data) and 32 <= data[end] < 127:
                end += 1
            return data[idx:end].decode("ascii", errors="replace")
    except OSError as e:
        logger.debug(f"Could not scan {binary_path} for signature marker: {e}")

    return None


def get_git_commit(repo_dir) -> Optional[str]:
    """Return the full HEAD commit SHA of a git checkout, or None."""
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10,
            check=False,
        )
        if result.returncode == 0:
            sha = result.stdout.decode().strip()
            return sha or None
    except (OSError, subprocess.SubprocessError) as e:
        logger.debug(f"Could not read git commit for {repo_dir}: {e}")
    return None


def write_binary_manifest(
    binary_path,
    *,
    tool: str,
    repo: Optional[str] = None,
    git_sha: Optional[str] = None,
    source_version: Optional[str] = None,
    extra: Optional[dict] = None,
) -> Optional[Path]:
    """
    Write a provenance/integrity manifest next to an installed binary.

    Records the source provenance (repo + git commit + source version), the
    SHA-256 of the built artifact, the binary's self-reported signature, and
    the build platform. Returns the manifest path, or None if writing failed
    (a failed manifest write is surfaced as a warning so install does not
    silently produce an unverifiable binary).
    """
    binary_path = Path(binary_path)
    platform_name, arch = detect_platform()
    try:
        manifest = {
            "manifest_version": MANIFEST_VERSION,
            "tool": tool,
            "repo": repo,
            "git_sha": git_sha,
            "source_version": source_version,
            "signature": read_build_signature(binary_path),
            "sha256": compute_sha256(binary_path),
            "size_bytes": binary_path.stat().st_size,
            "platform": f"{platform_name}/{arch}",
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        if extra:
            manifest.update(extra)

        path = manifest_path_for(binary_path)
        with open(path, "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
            f.write("\n")
        logger.info(f"Wrote provenance manifest: {path} (git {git_sha or 'unknown'})")
        return path
    except (OSError, ValueError) as e:
        # Surface, do not swallow: an unverifiable binary is a real problem.
        logger.warning(
            f"Could not write provenance manifest for {binary_path}: {e}. "
            f"The binary will install but its provenance cannot be verified at runtime."
        )
        return None


def verify_binary_manifest(binary_path) -> Tuple[str, str]:
    """
    Verify an installed binary against its provenance manifest.

    Returns ``(status, message)`` where status is one of:
      - "ok"         : manifest present and SHA-256 matches.
      - "unverified" : no manifest beside the binary (e.g. a user-supplied
                       binary on PATH). Caller should surface a visible warning.
      - "mismatch"   : manifest present but corrupt, or SHA-256 differs from the
                       recorded value. Caller MUST surface a visible error and
                       refuse to run the binary silently.

    The message is human-readable and names the concrete discrepancy.
    """
    binary_path = Path(binary_path)
    path = manifest_path_for(binary_path)

    if not path.exists():
        return (
            "unverified",
            f"no provenance manifest next to {binary_path.name} "
            f"({path.name} missing) — cannot confirm this is the satellome build",
        )

    try:
        with open(path) as f:
            manifest = json.load(f)
    except (OSError, ValueError) as e:
        return (
            "mismatch",
            f"provenance manifest {path.name} is unreadable/corrupt: {e}",
        )

    expected = manifest.get("sha256")
    if not expected:
        return (
            "mismatch",
            f"provenance manifest {path.name} has no sha256 field (corrupt or wrong schema)",
        )

    try:
        actual = compute_sha256(binary_path)
    except OSError as e:
        return ("mismatch", f"cannot read {binary_path} to verify checksum: {e}")

    if actual != expected:
        return (
            "mismatch",
            f"{binary_path.name} sha256 mismatch — manifest={expected[:12]}… "
            f"actual={actual[:12]}… (binary was replaced, corrupted, or rebuilt "
            f"out of band)",
        )

    sig = manifest.get("signature") or "n/a"
    git_sha = manifest.get("git_sha") or "unknown"
    return ("ok", f"{manifest.get('tool', binary_path.name)} verified: {sig} (git {git_sha})")


def resolve_binary(name: str) -> Optional[str]:
    """
    Resolve a managed binary to a single absolute path, deterministically.

    Resolution order — managed dir first, then system PATH — so that the
    satellome-installed (manifest-backed) binary always wins over an unrelated
    one that happens to be on PATH. This unifies what used to be two divergent
    resolvers (``shutil.which`` first in main.py vs managed-first in
    ``check_binary_exists``), which could otherwise pick different binaries.

    Returns the absolute path, or None if the binary is not found anywhere.
    """
    managed = get_satellome_bin_dir() / name
    if managed.exists() and os.access(managed, os.X_OK):
        return str(managed)

    system_path = shutil.which(name)
    if system_path:
        return system_path

    return None
