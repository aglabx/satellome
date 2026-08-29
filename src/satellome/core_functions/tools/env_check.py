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
# Any non-empty value stops satellome from editing shell startup files on its own.
ENV_NO_FIX = "SATELLOME_NO_PATH_FIX"

# Marker around the block we own in a shell startup file. Its presence is what
# makes the edit idempotent: we append exactly once and never touch it again.
PATH_BLOCK_START = "# >>> satellome path >>>"
PATH_BLOCK_END = "# <<< satellome path <<<"

CONSOLE_SCRIPT = "satellome"
SEPARATOR = "-" * 70

# What every external tool costs when it is absent. Reporting a tool as merely
# "not installed" is useless to someone who then has a run produce fewer
# results than they expected: the consequence is the part that matters, and it
# differs sharply per tool. Three classes:
#   "required" - the analysis cannot run at all
#   "results"  - the run succeeds but a step is SKIPPED, so output is missing
#   "speedup"  - a Python fallback produces the same result, more slowly
#   "optional" - only needed for a non-default mode
RUST_TOOLS_INSTALL = "satellome --install-rust-tools"


@dataclass(frozen=True)
class ToolSpec:
    name: str
    kind: str          # "managed" (installer-provided) or "script" (pip console script)
    impact_class: str  # required | results | speedup | optional
    impact: str        # what the user loses without it
    install: str       # the command that provides it

    @property
    def breaks_results(self) -> bool:
        """True when a missing tool silently changes what a run produces."""
        return self.impact_class in ("required", "results")


TOOL_SPECS = (
    ToolSpec("fastan", "managed", "required",
             "the default tandem-repeat search cannot run",
             "satellome --install-fastan"),
    ToolSpec("tanbed", "managed", "required",
             "monomer/array extraction cannot run",
             "satellome --install-tanbed"),
    ToolSpec("trf", "managed", "optional",
             "only needed with --run-trf (TRF is off by default)",
             "satellome --install-trf"),
    ToolSpec("arraysplitter", "script", "required",
             "consensus/HOR calculation fails and the FasTAN step aborts",
             "pip install arraysplitter"),
    ToolSpec("sat-family", "managed", "results",
             "satellite family clustering is SKIPPED - no families output",
             RUST_TOOLS_INSTALL),
    ToolSpec("telomere-check", "managed", "results",
             "the telomere check is SKIPPED - no telomere output",
             RUST_TOOLS_INSTALL),
    ToolSpec("find-gaps", "managed", "speedup",
             "gap detection falls back to Python (same result, slower)",
             RUST_TOOLS_INSTALL),
    ToolSpec("bed-extract", "managed", "speedup",
             "sequence extraction falls back to Python (same result, slower)",
             RUST_TOOLS_INSTALL),
    ToolSpec("genome-size", "managed", "speedup",
             "genome size falls back to Python (same result, slower)",
             RUST_TOOLS_INSTALL),
)

TOOL_SPECS_BY_NAME = {spec.name: spec for spec in TOOL_SPECS}

# Kept for callers that only need the names.
PIP_COMPANIONS = tuple(s.name for s in TOOL_SPECS if s.kind == "script")
MANAGED_BINARIES = tuple(s.name for s in TOOL_SPECS if s.kind == "managed")


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
    def spec(self) -> Optional["ToolSpec"]:
        return TOOL_SPECS_BY_NAME.get(self.name)

    @property
    def found(self) -> bool:
        return self.path is not None

    @property
    def built_for(self) -> Optional[str]:
        return binary_platform(self.path) if self.found else None

    @property
    def wrong_arch(self) -> bool:
        """Present and executable, but cannot run on this machine."""
        return self.found and runnable_here(self.path) is False

    @property
    def usable(self) -> bool:
        return self.found and not self.wrong_arch

    @property
    def breaks_results(self) -> bool:
        """Unusable, and that changes what a run produces.

        A binary built for another platform counts as missing here: it is worse
        than missing, since it is reported as installed and fails only when the
        pipeline finally tries to execute it.
        """
        spec = self.spec
        return not self.usable and spec is not None and spec.breaks_results

    @property
    def hidden(self) -> bool:
        """Installed, but the user's shell cannot see it."""
        return self.found and not self.on_path and self.source == "scripts-dir"


@dataclass
class PathFixResult:
    """Outcome of trying to put a directory on the user's PATH for good."""

    status: str = "no-target"
    # "added"              - the block was appended to a startup file just now
    # "already-configured" - a previous run (or the user) already added it
    # "already-on-path"    - nothing to do
    # "no-target"          - no shell startup file could be determined
    # "error"              - the startup file exists but could not be read/written
    directory: Optional[str] = None
    rc_file: Optional[str] = None
    detail: str = ""

    @property
    def changed(self) -> bool:
        return self.status == "added"

    @property
    def ok(self) -> bool:
        return self.status in ("added", "already-configured", "already-on-path")

    @property
    def activate_command(self) -> Optional[str]:
        """What the *current* shell needs, since we cannot change its PATH."""
        if self.rc_file and self.status in ("added", "already-configured"):
            return f"source {_tilde(self.rc_file)}"
        return None


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
    path_fix: Optional["PathFixResult"] = None

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


# Executable-format magic. A file can be present, executable and completely
# unable to run: five of the eight binaries on a production Linux box turned
# out to be macOS arm64 builds copied from a developer's laptop, and satellome
# reported them as installed right up until the run died 65 minutes in with
# "OSError: Exec format error".
ELF_MAGIC = b"\x7fELF"
MACHO_MAGICS = (b"\xfe\xed\xfa\xce", b"\xfe\xed\xfa\xcf",
                b"\xce\xfa\xed\xfe", b"\xcf\xfa\xed\xfe")
MACHO_UNIVERSAL = (b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca")
PE_MAGIC = b"MZ"

# ELF e_machine values we care about.
ELF_MACHINES = {0x3E: "x86_64", 0xB7: "arm64", 0x28: "arm", 0x03: "x86"}
# Mach-O cpu types (low 32 bits, ABI64 bit stripped).
MACHO_CPUS = {7: "x86_64", 12: "arm64"}


def binary_platform(path) -> Optional[str]:
    """What this executable was built for, e.g. "linux-x86_64", "macos-arm64".

    Returns None when the format is unknown or the file cannot be read - which
    includes shell scripts and other interpreted entry points, where the
    question does not apply.
    """
    try:
        with open(str(path), "rb") as handle:
            head = handle.read(24)
    except OSError:
        return None

    if len(head) < 4:
        return None

    if head.startswith(b"#!"):
        return None  # a script: its interpreter decides, not the file

    if head.startswith(ELF_MAGIC):
        if len(head) < 20:
            return "linux"
        machine = int.from_bytes(head[18:20], "little" if head[5] == 1 else "big")
        return f"linux-{ELF_MACHINES.get(machine, f'machine{machine:#x}')}"

    if head.startswith(MACHO_UNIVERSAL):
        return "macos-universal"

    if head[:4] in MACHO_MAGICS:
        big_endian = head[:4] in (b"\xfe\xed\xfa\xce", b"\xfe\xed\xfa\xcf")
        cpu = int.from_bytes(head[4:8], "big" if big_endian else "little")
        return f"macos-{MACHO_CPUS.get(cpu & 0x00FFFFFF, f'cpu{cpu:#x}')}"

    if head.startswith(PE_MAGIC):
        return "windows"

    return None


def host_platform() -> str:
    """This machine, in the same vocabulary as ``binary_platform``."""
    from satellome.installers.base import detect_platform

    name, arch = detect_platform()
    return f"{'macos' if name == 'darwin' else name}-{arch}"


def runnable_here(path) -> Optional[bool]:
    """True/False if the file's format settles it, None if it cannot be told."""
    built_for = binary_platform(path)
    if built_for is None:
        return None
    if built_for == "macos-universal":
        return host_platform().startswith("macos")
    host = host_platform()
    if "unknown" in host:
        return None
    # Compare OS and architecture separately so an unrecognised machine value
    # does not masquerade as a mismatch.
    return built_for == host


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


def _tilde(path, dollar: bool = False) -> str:
    """Render a path under $HOME as ``~/x`` (or ``$HOME/x``) for display."""
    home = os.path.expanduser("~")
    text = str(path)
    if text == home or text.startswith(home + os.sep):
        return ("$HOME" if dollar else "~") + text[len(home):]
    return text


def _shell_fix_lines(directory) -> List[str]:
    """Copy-pasteable ways out, most convenient first."""
    shown = _tilde(directory, dollar=True)
    return [
        f'    export PATH="{shown}:$PATH"      # this shell only',
        f"    echo 'export PATH=\"{shown}:$PATH\"' >> ~/.bashrc   # permanent",
        f"    python -m {CONSOLE_SCRIPT} ...   # no PATH change needed",
    ]


def shell_rc_files() -> List[Path]:
    """Startup files for the user's login shell, most specific first.

    The first entry is what we append to; the rest are only inspected so an
    entry the user already added by hand is not duplicated.
    """
    shell = os.path.basename(os.environ.get("SHELL") or "")
    home = Path.home()

    if shell == "zsh":
        base = Path(os.environ.get("ZDOTDIR") or home)
        return [base / ".zshrc"]
    if shell == "fish":
        return [home / ".config" / "fish" / "config.fish"]
    if shell == "bash":
        # A macOS Terminal tab is a login shell and reads .bash_profile; most
        # Linux terminals read .bashrc. Write to the one the shell will read.
        if sys.platform == "darwin":
            return [home / ".bash_profile", home / ".bashrc", home / ".profile"]
        return [home / ".bashrc", home / ".bash_profile", home / ".profile"]
    return [home / ".profile"]


def _path_block(directory: str) -> str:
    """The lines we append. Marked on both sides so it can be found and removed."""
    shown = _tilde(directory, dollar=True)
    if os.path.basename(os.environ.get("SHELL") or "") == "fish":
        body = f"fish_add_path {shown}"
    else:
        body = f'export PATH="{shown}:$PATH"'
    return f"\n{PATH_BLOCK_START}\n{body}\n{PATH_BLOCK_END}\n"


def _already_configured(text: str, directory: str) -> bool:
    """True if this file already puts ``directory`` on PATH, ours or the user's."""
    if PATH_BLOCK_START in text:
        return True
    forms = {directory, _tilde(directory), _tilde(directory, dollar=True)}
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if "PATH" not in stripped and "fish_add_path" not in stripped:
            continue
        if any(form in stripped for form in forms):
            return True
    return False


def ensure_path_entry(directory, rc_file=None) -> PathFixResult:
    """Put ``directory`` on the PATH of future shells, exactly once.

    This is deliberately not a silent convenience: the caller reports what was
    written and where. A process cannot change its parent shell's environment,
    so the current shell still needs ``source <rc>`` — the result carries that
    command rather than leaving the user to wonder why nothing changed.
    """
    directory = str(directory)
    if dir_on_path(directory):
        return PathFixResult(status="already-on-path", directory=directory)

    targets = [Path(rc_file)] if rc_file else shell_rc_files()
    if not targets:
        return PathFixResult(
            status="no-target",
            directory=directory,
            detail="could not determine a shell startup file for $SHELL",
        )

    # An existing entry anywhere wins: never write a second one.
    existing = []
    for candidate in targets:
        try:
            text = candidate.read_text(encoding="utf-8", errors="replace")
        except FileNotFoundError:
            continue
        except OSError as error:
            # The file is there but unreadable — that is a real problem the user
            # must see, not a reason to quietly append to a different file.
            return PathFixResult(
                status="error",
                directory=directory,
                rc_file=str(candidate),
                detail=f"cannot read {candidate}: {error}",
            )
        existing.append(candidate)
        if _already_configured(text, directory):
            return PathFixResult(
                status="already-configured", directory=directory, rc_file=str(candidate)
            )

    # Append to a file the shell already reads. Creating a startup file that did
    # not exist can change which files the shell reads at all — on macOS a new
    # ~/.bash_profile makes a login shell stop reading ~/.profile — so a file
    # that is already there is always the safer target.
    target = existing[0] if existing else targets[0]
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(str(target), "a", encoding="utf-8") as handle:
            handle.write(_path_block(directory))
    except OSError as error:
        return PathFixResult(
            status="error",
            directory=directory,
            rc_file=str(target),
            detail=f"cannot write {target}: {error}",
        )

    return PathFixResult(status="added", directory=directory, rc_file=str(target))


def fix_entrypoint_path(report=None, log=None) -> Optional[PathFixResult]:
    """Add this interpreter's scripts directory to PATH if it is missing.

    Returns ``None`` when there is nothing to fix (no launcher, or already
    visible). Adding a directory cannot repair a *shadowed* launcher — two
    competing installs need one of them removed — so that case is left alone.
    """
    log = log or logger
    report = report or check_satellome_entrypoint()

    if not report.script_for_interpreter or report.script_on_path:
        return None

    directory = os.path.dirname(report.script_for_interpreter)
    result = ensure_path_entry(directory)

    if result.status == "added":
        log.warning(f"Added {_tilde(directory)} to your PATH in {_tilde(result.rc_file)}")
        log.warning(f"It applies to new shells. For this one: {result.activate_command}")
    elif result.status == "already-configured":
        log.warning(
            f"{_tilde(directory)} is already configured in {_tilde(result.rc_file)}, "
            "but this shell has not picked it up yet. Run: "
            f"{result.activate_command}   (or open a new shell)"
        )
    elif result.status == "error":
        log.error(f"Could not put {_tilde(directory)} on your PATH: {result.detail}")
        for line in _shell_fix_lines(directory):
            log.error(line)
    elif result.status == "no-target":
        log.warning(f"Could not put {_tilde(directory)} on your PATH: {result.detail}")
        for line in _shell_fix_lines(directory):
            log.warning(line)

    return result


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
            "fail with 'command not found'. Equivalent manual fixes:\n"
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

    for spec in TOOL_SPECS:
        if spec.kind == "managed" and managed_dir is not None:
            candidate = Path(managed_dir) / spec.name
            if _executable(candidate):
                on_path = dir_on_path(managed_dir)
                locations.append(
                    ToolLocation(
                        name=spec.name,
                        path=str(candidate),
                        source="managed",
                        on_path=on_path,
                        dir_off_path=None if on_path else str(managed_dir),
                    )
                )
                continue
        locations.append(find_console_script(spec.name))

    return locations


def missing_tools(tools: List[ToolLocation]) -> List[ToolLocation]:
    """Absent tools whose absence changes results, worst class first."""
    order = {"required": 0, "results": 1}
    broken = [t for t in tools if t.breaks_results]
    return sorted(broken, key=lambda t: order.get(t.spec.impact_class, 9))


def degraded_tools(tools: List[ToolLocation]) -> List[ToolLocation]:
    """Absent tools that only cost speed; results are unchanged."""
    return [
        t for t in tools
        if not t.usable and t.spec is not None and t.spec.impact_class == "speedup"
    ]


def hidden_tool_warning(location: ToolLocation) -> str:
    """One line telling the user a tool exists but their shell cannot see it."""
    return (
        f"{location.name}: found at {location.path}, but {location.dir_off_path} "
        "is not on your PATH. Satellome will use it by absolute path; your shell "
        "will not find it."
    )


def warn_if_entrypoint_misconfigured(log=None, fix=True) -> EntrypointReport:
    """Cheap startup check: surface PATH problems and repair the fixable one.

    A launcher that is merely invisible is repaired in place — the scripts
    directory is appended, once, to the user's shell startup file — because
    telling someone to edit ``~/.bashrc`` on every run of a genome pipeline is
    a worse answer than doing it and saying so. A *shadowed* launcher is only
    reported: choosing which of two installs to delete is not ours to make.

    ``SATELLOME_NO_ENV_CHECK`` silences the whole check;
    ``SATELLOME_NO_PATH_FIX`` keeps the warning but never touches a file.
    """
    log = log or logger
    report = check_satellome_entrypoint()
    if os.environ.get(ENV_DISABLE):
        return report

    for problem in report.problems:
        for line in problem.splitlines():
            log.warning(line)

    if report.problems:
        if fix and not os.environ.get(ENV_NO_FIX):
            report.path_fix = fix_entrypoint_path(report=report, log=log)
        log.warning(
            f"(silence this check with {ENV_DISABLE}=1; "
            f"keep it but never edit shell files with {ENV_NO_FIX}=1; "
            "details: 'satellome --doctor')"
        )
    return report


def warn_about_missing_tools(log=None) -> List[ToolLocation]:
    """Say up front which results this run will not produce.

    Called before the pipeline starts, because a genome run takes hours and
    "sat-family binary not found, skipping" scrolling past mid-run is how a
    user ends up with an output directory that quietly lacks family clustering.
    Returns the missing result-affecting tools.
    """
    log = log or logger
    if os.environ.get(ENV_DISABLE):
        return []

    tools = check_companion_tools()
    broken = missing_tools(tools)
    if not broken:
        return []

    log.warning(SEPARATOR)
    log.warning(f"{len(broken)} tool(s) unusable - THIS RUN WILL PRODUCE LESS THAN A COMPLETE ONE:")
    for tool in broken:
        if tool.wrong_arch:
            log.warning(
                f"  {tool.name}: built for {tool.built_for}, this machine is "
                f"{host_platform()} - {tool.spec.impact}"
            )
        else:
            log.warning(f"  {tool.name}: {tool.spec.impact}")
    for command in sorted({t.spec.install for t in broken}):
        log.warning(f"  install with: {command}")
    log.warning("Full report: satellome --doctor")
    log.warning(SEPARATOR)
    return broken


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
        spec = tool.spec
        if not tool.found:
            impact = f" -> {spec.impact}" if spec else ""
            label = "MISSING" if (spec and spec.breaks_results) else "missing"
            lines.append(f"  {tool.name:<14}: {label}{impact}")
        elif tool.wrong_arch:
            lines.append(
                f"  {tool.name:<14}: WRONG ARCHITECTURE - built for {tool.built_for}, "
                f"this machine is {host_platform()}"
            )
            lines.append(f"  {'':<14}  {tool.path}")
        elif tool.hidden:
            lines.append(f"  {tool.name:<14}: {tool.path}  [off PATH - usable by satellome only]")
        else:
            lines.append(f"  {tool.name:<14}: {tool.path}  [{tool.source}]")

    broken = missing_tools(tools)
    slow = degraded_tools(tools)
    hidden = [t for t in tools if t.hidden]

    if slow:
        lines.append("")
        lines.append("SLOWER WITHOUT (results are identical)")
        for tool in slow:
            lines.append(f"  {tool.name:<14}: {tool.spec.impact}")
        installs = sorted({t.spec.install for t in slow})
        for command in installs:
            lines.append(f"  install with: {command}")

    problem_count = len(report.problems) + len(hidden) + len(broken)
    lines.append("")
    if problem_count:
        lines.append(f"PROBLEMS ({problem_count})")
        # Missing tools first: they are what makes a run produce less than the
        # user thinks it did, which no amount of PATH advice will explain.
        for tool in broken:
            if tool.wrong_arch:
                lines.append(
                    f"  ! {tool.name} cannot run here: built for {tool.built_for}, "
                    f"this machine is {host_platform()} - {tool.spec.impact}"
                )
            else:
                lines.append(f"  ! {tool.name} is not installed - {tool.spec.impact}")
        if broken:
            for command in sorted({t.spec.install for t in broken}):
                lines.append(f"    install with: {command}")
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


def run_doctor(log=None, fix=False) -> bool:
    """Print the environment report. Returns True when nothing is wrong.

    With ``fix=True`` the repairable problem (a launcher directory missing from
    PATH) is also written to the shell startup file before the report is
    rendered, so the report reflects what was done.
    """
    log = log or logger
    report = check_satellome_entrypoint()
    applied = fix_entrypoint_path(report=report, log=log) if fix else None
    tools = check_companion_tools()
    # A run with sat-family missing finishes "successfully" and silently lacks
    # family clustering. Calling that healthy is the failure this report exists
    # to prevent, so a missing result-affecting tool is a non-zero exit.
    healthy = report.ok and not any(t.hidden for t in tools) and not missing_tools(tools)
    emit = log.info if healthy else log.warning
    for line in format_doctor_report(report, tools):
        emit(line)
    if applied and applied.activate_command:
        log.warning("")
        log.warning(
            f"PATH updated in {_tilde(applied.rc_file)}. New shells get it "
            f"automatically; for this one run: {applied.activate_command}"
        )
    return healthy


def run_fix_path(log=None) -> bool:
    """Repair the PATH for future shells. Returns True when PATH is settled."""
    log = log or logger
    report = check_satellome_entrypoint()

    result = fix_entrypoint_path(report=report, log=log)
    if result is None:
        if not report.script_for_interpreter:
            log.info(
                f"No '{CONSOLE_SCRIPT}' launcher belongs to this interpreter "
                f"({sys.executable}) — nothing to add to PATH. This is normal "
                "for a source checkout; use 'python -m satellome'."
            )
            return True
        log.info(
            f"'{CONSOLE_SCRIPT}' is already on your PATH "
            f"({report.script_on_path}) — nothing to do."
        )

    # Hidden companion tools live in the same directories; fix those too.
    for tool in check_companion_tools():
        if tool.hidden:
            outcome = ensure_path_entry(tool.dir_off_path)
            if outcome.status == "added":
                log.warning(
                    f"Added {_tilde(tool.dir_off_path)} to your PATH in "
                    f"{_tilde(outcome.rc_file)} (needed by {tool.name})"
                )
            elif outcome.status == "error":
                log.error(
                    f"Could not put {_tilde(tool.dir_off_path)} on your PATH "
                    f"({tool.name}): {outcome.detail}"
                )

    if report.problems and not (result and result.ok):
        return False
    # A shadowed launcher is a problem PATH edits cannot solve.
    return report.ok or bool(result and result.changed)
