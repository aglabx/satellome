"""Installation / PATH self-diagnosis for Satellome.

Motivation
----------
``pip install --user satellome`` prints::

    WARNING: The script satellome is installed in '/home/user/.local/bin'
    which is not on PATH.

and then the user types ``satellome`` and gets ``command not found`` — the
package is installed correctly, only the shell cannot see the launcher. The
same trap hits every companion console script Satellome shells out to
(``arraysplitter`` is pip-installed at runtime): ``shutil.which`` returns
``None`` even though the executable exists in this interpreter's scripts
directory, so a pipeline step is silently skipped or "installed but not found
in PATH" is reported forever.

Silent-error policy (see project guidelines):
- A tool that is genuinely *not installed* is an expected absence -> optional.
- A tool that IS installed but invisible because of PATH is a misconfiguration
  the user must be told about: we resolve it anyway (so the run does not
  degrade) and we say so, with the exact command that fixes the shell.

Nothing in this module executes the tools it inspects; it only looks at the
filesystem, ``sysconfig`` and the shebang line, so it is cheap enough to run at
startup of every pipeline run.
"""

import logging
import os
import shutil
import sys
import sysconfig
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

# Any non-empty value silences the startup PATH warning (--doctor still works).
ENV_DISABLE = "SATELLOME_NO_ENV_CHECK"

CONSOLE_SCRIPT = "satellome"

# Companion console scripts / binaries Satellome shells out to. ``pip`` marks
# the ones installed as Python console scripts, i.e. the ones exposed to the
# very same "--user install is not on PATH" failure.
PIP_COMPANIONS = ("arraysplitter",)
MANAGED_BINARIES = ("fastan", "tanbed", "trf")
BUNDLED_BINARIES = ("sat-family", "telomere-check", "find-gaps", "bed-extract", "genome-size")


@dataclass
class ToolLocation:
    """Where an executable actually is, and whether the shell can see it."""

    name: str
    path: Optional[str] = None
    # "path"        - found through $PATH (the normal case)
    # "scripts-dir" - found in this interpreter's scripts dir, NOT on $PATH
    # "managed"     - found in the satellome-managed bin dir
    # "missing"     - not found anywhere
    source: str = "missing"
    on_path: bool = False
    dir_off_path: Optional[str] = None

    @property
    def found(self) -> bool:
        return self.path is not None

    @property
    def hidden(self) -> bool:
        """Installed, but the user's shell cannot see it."""
        return self.found and not self.on_path and self.source == "scripts-dir"


@dataclass
class EntrypointReport:
    """State of the ``satellome`` launcher itself."""

    version: str = ""
    package_dir: str = ""
    interpreter: str = ""
    script_on_path: Optional[str] = None      # what `satellome` in the shell runs
    script_for_interpreter: Optional[str] = None  # launcher belonging to THIS python
    path_script_interpreter: Optional[str] = None  # shebang of script_on_path
    problems: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.problems


def _script_names(name: str) -> List[str]:
    if os.name == "nt":
        return [name + ".exe", name + "-script.py", name]
    return [name]


def _user_scheme() -> Optional[str]:
    """Name of the sysconfig scheme used by ``pip install --user``."""
    try:
        return sysconfig.get_preferred_scheme("user")  # Python >= 3.10
    except (AttributeError, KeyError, ValueError):
        pass
    if os.name == "nt":
        return "nt_user"
    if sys.platform == "darwin" and getattr(sys, "_framework", None):
        return "osx_framework_user"
    return "posix_user"


def interpreter_script_dirs() -> List[Path]:
    """Directories where pip puts console scripts for the running interpreter.

    Covers the three layouts that matter: a venv/conda ``bin`` next to the
    interpreter, the interpreter's own scripts scheme, and the ``--user``
    scheme (``~/.local/bin``) that produced the warning this module exists for.
    """
    dirs: List[Path] = []

    def add(value) -> None:
        if not value:
            return
        candidate = Path(value)
        if candidate not in dirs:
            dirs.append(candidate)

    try:
        add(Path(sys.executable).resolve().parent)
    except OSError:  # unreadable /proc-style symlink; not fatal for a hint
        add(Path(sys.executable).parent)

    for scheme in (None, _user_scheme()):
        try:
            add(sysconfig.get_path("scripts") if scheme is None
                else sysconfig.get_path("scripts", scheme=scheme))
        except (KeyError, ValueError, OSError):
            continue

    try:
        import site

        user_base = site.getuserbase()
    except (ImportError, AttributeError, OSError):
        user_base = None
    if user_base:
        add(Path(user_base) / ("Scripts" if os.name == "nt" else "bin"))

    return dirs


def path_entries() -> List[str]:
    """Normalised $PATH directories (realpath'd, empty entries dropped)."""
    entries = []
    for raw in (os.environ.get("PATH") or "").split(os.pathsep):
        if not raw:
            continue
        try:
            entries.append(os.path.realpath(os.path.expanduser(raw)))
        except OSError:
            entries.append(raw)
    return entries


def dir_on_path(directory) -> bool:
    """True if ``directory`` is one of the $PATH entries."""
    if directory is None:
        return False
    try:
        target = os.path.realpath(os.path.expanduser(str(directory)))
    except OSError:
        target = str(directory)
    return target in path_entries()


def _executable(candidate: Path) -> bool:
    return candidate.is_file() and os.access(str(candidate), os.X_OK)


def find_console_script(name: str) -> ToolLocation:
    """Locate a Python console script, PATH first, interpreter scripts dir second.

    Unlike bare ``shutil.which`` this still finds a script that pip installed
    into ``~/.local/bin`` on a machine whose PATH lacks that directory — the
    exact case pip warns about at install time. The returned location records
    that the shell cannot see it so callers can surface it instead of behaving
    as if the tool were missing.
    """
    found = shutil.which(name)
    if found:
        return ToolLocation(name=name, path=os.path.abspath(found), source="path", on_path=True)

    for directory in interpreter_script_dirs():
        for script_name in _script_names(name):
            candidate = directory / script_name
            if _executable(candidate):
                return ToolLocation(
                    name=name,
                    path=str(candidate),
                    source="scripts-dir",
                    on_path=False,
                    dir_off_path=str(directory),
                )

    return ToolLocation(name=name)


def script_shebang_interpreter(script_path) -> Optional[str]:
    """Interpreter a console script will run under, read from its shebang.

    Returns ``None`` for a binary launcher (Windows ``.exe``), an unreadable
    file or a script without a shebang — all cases where we simply cannot tell,
    and must not pretend a mismatch.
    """
    try:
        with open(str(script_path), "rb") as handle:
            first_line = handle.readline(512)
    except OSError as error:
        logger.debug("cannot read shebang of %s: %s", script_path, error)
        return None

    if not first_line.startswith(b"#!"):
        return None

    parts = first_line[2:].decode("utf-8", "replace").strip().split()
    if not parts:
        return None
    # "#!/usr/bin/env python3" -> the interpreter is the second token
    if os.path.basename(parts[0]) == "env" and len(parts) > 1:
        return parts[1]
    return parts[0]


def _same_file(left, right) -> bool:
    if not left or not right:
        return False
    try:
        return os.path.samefile(str(left), str(right))
    except OSError:
        return os.path.realpath(str(left)) == os.path.realpath(str(right))


def _shell_fix_lines(directory) -> List[str]:
    """Copy-pasteable ways out, most convenient first."""
    home = os.path.expanduser("~")
    shown = str(directory)
    if shown.startswith(home + os.sep):
        shown = "$HOME" + shown[len(home):]
    return [
        f'    export PATH="{shown}:$PATH"      # this shell only',
        f"    echo 'export PATH=\"{shown}:$PATH\"' >> ~/.bashrc   # permanent",
        f"    python -m {CONSOLE_SCRIPT} ...   # no PATH change needed",
    ]


def check_satellome_entrypoint() -> EntrypointReport:
    """Diagnose the ``satellome`` launcher: present? visible? the right one?"""
    import satellome

    report = EntrypointReport(
        version=getattr(satellome, "__version__", "unknown"),
        package_dir=str(Path(satellome.__file__).resolve().parent),
        interpreter=sys.executable,
    )

    on_path = shutil.which(CONSOLE_SCRIPT)
    report.script_on_path = os.path.abspath(on_path) if on_path else None

    for directory in interpreter_script_dirs():
        for script_name in _script_names(CONSOLE_SCRIPT):
            candidate = directory / script_name
            if _executable(candidate):
                report.script_for_interpreter = str(candidate)
                break
        if report.script_for_interpreter:
            break

    if report.script_on_path:
        report.path_script_interpreter = script_shebang_interpreter(report.script_on_path)

    own = report.script_for_interpreter

    if own and not report.script_on_path:
        # The pip warning, one step later: installed, invisible to the shell.
        directory = os.path.dirname(own)
        report.problems.append(
            f"The '{CONSOLE_SCRIPT}' command is installed at {own} but that "
            f"directory is not on your PATH, so typing '{CONSOLE_SCRIPT}' will "
            "fail with 'command not found'. Fix with any one of:\n"
            + "\n".join(_shell_fix_lines(directory))
        )
    elif own and report.script_on_path and not _same_file(own, report.script_on_path):
        # Two installs. The one the user types is not the one running now.
        report.problems.append(
            f"Your shell resolves '{CONSOLE_SCRIPT}' to {report.script_on_path}, "
            f"but this interpreter's launcher is {own}. Two installations are "
            "competing and 'satellome' may run a different version than "
            f"'{os.path.basename(sys.executable)} -m {CONSOLE_SCRIPT}'. Remove the "
            f"stale one, or always call: python -m {CONSOLE_SCRIPT} ..."
        )
    elif not own and not report.script_on_path:
        report.notes.append(
            f"No '{CONSOLE_SCRIPT}' launcher found — running from a source "
            f"checkout or via 'python -m {CONSOLE_SCRIPT}'. That works; only the "
            "bare 'satellome' command is unavailable."
        )

    if (
        report.script_on_path
        and report.path_script_interpreter
        and not _same_file(report.path_script_interpreter, sys.executable)
    ):
        detail = (
            f"The '{CONSOLE_SCRIPT}' on your PATH runs "
            f"{report.path_script_interpreter}, not the interpreter you are using "
            f"now ({sys.executable}). Its dependencies and version can differ. "
            f"Use 'python -m {CONSOLE_SCRIPT}' to always follow the active "
            "environment."
        )
        # When shadowing was already reported these are two faces of one
        # misconfiguration; keep it as a single problem so the fix is one action.
        if report.problems and "competing" in report.problems[-1]:
            report.problems[-1] += "\n    " + detail
        else:
            report.problems.append(detail)

    return report


def check_companion_tools() -> List[ToolLocation]:
    """Locate every external tool Satellome may shell out to."""
    from satellome.installers.base import get_satellome_bin_dir

    locations: List[ToolLocation] = []

    try:
        managed_dir = get_satellome_bin_dir()
    except OSError as error:
        # Unwritable home/package dir: report it rather than pretending the
        # managed binaries are simply absent.
        logger.warning("cannot access the satellome bin directory: %s", error)
        managed_dir = None

    for name in MANAGED_BINARIES + BUNDLED_BINARIES:
        if managed_dir is not None:
            candidate = Path(managed_dir) / name
            if _executable(candidate):
                locations.append(
                    ToolLocation(
                        name=name,
                        path=str(candidate),
                        source="managed",
                        on_path=dir_on_path(managed_dir),
                        dir_off_path=None if dir_on_path(managed_dir) else str(managed_dir),
                    )
                )
                continue
        locations.append(find_console_script(name))

    for name in PIP_COMPANIONS:
        locations.append(find_console_script(name))

    return locations


def hidden_tool_warning(location: ToolLocation) -> str:
    """One line telling the user a tool exists but their shell cannot see it."""
    return (
        f"{location.name}: found at {location.path}, but {location.dir_off_path} "
        "is not on your PATH. Satellome will use it by absolute path; your shell "
        "will not find it."
    )


def warn_if_entrypoint_misconfigured(log=None) -> EntrypointReport:
    """Cheap startup check: surface PATH problems, never abort the run.

    Returns the report so callers can record it; set ``SATELLOME_NO_ENV_CHECK``
    to keep the pipeline output quiet on machines where this is known.
    """
    log = log or logger
    report = check_satellome_entrypoint()
    if os.environ.get(ENV_DISABLE):
        return report
    for problem in report.problems:
        for line in problem.splitlines():
            log.warning(line)
        log.warning(f"(silence this check with {ENV_DISABLE}=1, or run 'satellome --doctor')")
    return report


def format_doctor_report(report: EntrypointReport, tools: List[ToolLocation]) -> List[str]:
    """Render the full ``--doctor`` output as log lines."""
    lines = [
        "SATELLOME ENVIRONMENT",
        f"  version        : {report.version}",
        f"  package        : {report.package_dir}",
        f"  interpreter    : {sys.executable}",
    ]

    if report.script_on_path:
        interp = report.path_script_interpreter or "unknown interpreter"
        lines.append(f"  'satellome' cmd: {report.script_on_path} (runs {interp})")
    else:
        lines.append("  'satellome' cmd: NOT on PATH")
    if report.script_for_interpreter:
        marker = "" if _same_file(report.script_for_interpreter, report.script_on_path) else "  <- not the one on PATH"
        lines.append(f"  launcher here  : {report.script_for_interpreter}{marker}")
    lines.append(f"  module fallback: {os.path.basename(sys.executable)} -m {CONSOLE_SCRIPT} (always works)")

    lines.append("")
    lines.append("SCRIPT DIRECTORIES")
    for directory in interpreter_script_dirs():
        state = "on PATH" if dir_on_path(directory) else "NOT on PATH"
        exists = "" if directory.is_dir() else " (does not exist)"
        lines.append(f"  {state:<11}: {directory}{exists}")

    lines.append("")
    lines.append("EXTERNAL TOOLS")
    for tool in tools:
        if not tool.found:
            lines.append(f"  {tool.name:<14}: not installed")
        elif tool.hidden:
            lines.append(f"  {tool.name:<14}: {tool.path}  [off PATH — usable by satellome only]")
        else:
            lines.append(f"  {tool.name:<14}: {tool.path}  [{tool.source}]")

    hidden = [t for t in tools if t.hidden]
    lines.append("")
    if report.problems or hidden:
        lines.append(f"PROBLEMS ({len(report.problems) + len(hidden)})")
        for problem in report.problems:
            first, *rest = problem.splitlines()
            lines.append(f"  ! {first}")
            lines.extend(f"  {line}" for line in rest)
        for tool in hidden:
            lines.append(f"  ! {hidden_tool_warning(tool)}")
            lines.extend(_shell_fix_lines(tool.dir_off_path)[:2])
    else:
        lines.append("No problems found.")

    for note in report.notes:
        lines.append(f"  note: {note}")

    return lines


def run_doctor(log=None) -> bool:
    """Print the environment report. Returns True when nothing is wrong."""
    log = log or logger
    report = check_satellome_entrypoint()
    tools = check_companion_tools()
    healthy = report.ok and not any(t.hidden for t in tools)
    emit = log.info if healthy else log.warning
    for line in format_doctor_report(report, tools):
        emit(line)
    return healthy
