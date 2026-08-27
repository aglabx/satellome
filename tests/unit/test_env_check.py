"""Tests for the installation / PATH self-diagnosis.

Each test pins a concrete regression path rather than "returns something":
the failure this module exists for is a tool that IS installed staying
invisible, so the assertions check that the tool is still resolved AND that the
misconfiguration is reported in words the user can act on.
"""

import os
import stat
import sys

import pytest

from satellome.core_functions.tools import env_check
from satellome.core_functions.tools.env_check import (
    ToolLocation,
    check_satellome_entrypoint,
    dir_on_path,
    find_console_script,
    format_doctor_report,
    hidden_tool_warning,
    interpreter_script_dirs,
    script_shebang_interpreter,
    warn_if_entrypoint_misconfigured,
)


def make_script(directory, name, shebang="#!/usr/bin/env python3\n"):
    """Create an executable console-script lookalike."""
    path = directory / name
    path.write_text(shebang + "print('hi')\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


class TestFindConsoleScript:
    def test_found_on_path(self, tmp_path, monkeypatch):
        script = make_script(tmp_path, "toolx")
        monkeypatch.setenv("PATH", str(tmp_path))

        location = find_console_script("toolx")

        assert location.path == str(script)
        assert location.source == "path"
        assert location.on_path is True
        assert location.hidden is False

    def test_found_in_scripts_dir_when_not_on_path(self, tmp_path, monkeypatch):
        """The pip '--user install is not on PATH' case: must still resolve."""
        scripts_dir = tmp_path / "local_bin"
        scripts_dir.mkdir()
        script = make_script(scripts_dir, "toolx")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [scripts_dir])

        location = find_console_script("toolx")

        assert location.path == str(script), "installed tool must not read as missing"
        assert location.source == "scripts-dir"
        assert location.on_path is False
        assert location.hidden is True, "must be flagged so the user is told, not silently used"
        assert location.dir_off_path == str(scripts_dir)

    def test_missing_everywhere(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PATH", str(tmp_path))
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [tmp_path])

        location = find_console_script("definitely-not-installed")

        assert location.found is False
        assert location.hidden is False
        assert location.source == "missing"

    def test_non_executable_file_is_not_a_tool(self, tmp_path, monkeypatch):
        (tmp_path / "toolx").write_text("not executable")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [tmp_path])

        assert find_console_script("toolx").found is False


class TestHiddenToolWarning:
    def test_message_names_tool_path_and_directory(self, tmp_path):
        location = ToolLocation(
            name="arraysplitter",
            path=str(tmp_path / "arraysplitter"),
            source="scripts-dir",
            on_path=False,
            dir_off_path=str(tmp_path),
        )

        message = hidden_tool_warning(location)

        assert "arraysplitter" in message
        assert str(tmp_path) in message
        assert "not on your PATH" in message


class TestPathHelpers:
    def test_dir_on_path_true_and_false(self, tmp_path, monkeypatch):
        other = tmp_path / "other"
        other.mkdir()
        monkeypatch.setenv("PATH", str(tmp_path))

        assert dir_on_path(tmp_path) is True
        assert dir_on_path(other) is False

    def test_dir_on_path_ignores_trailing_separator(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PATH", str(tmp_path) + os.sep)
        assert dir_on_path(tmp_path) is True

    def test_interpreter_script_dirs_includes_interpreter_bin(self):
        dirs = [str(d) for d in interpreter_script_dirs()]
        assert str(os.path.dirname(os.path.realpath(sys.executable))) in dirs
        assert len(dirs) == len(set(dirs)), "no duplicate directories"


class TestShebang:
    def test_env_shebang_resolves_to_second_token(self, tmp_path):
        script = make_script(tmp_path, "toolx", "#!/usr/bin/env python3.11\n")
        assert script_shebang_interpreter(script) == "python3.11"

    def test_absolute_shebang(self, tmp_path):
        script = make_script(tmp_path, "toolx", "#!/opt/py/bin/python\n")
        assert script_shebang_interpreter(script) == "/opt/py/bin/python"

    def test_no_shebang_returns_none_not_a_guess(self, tmp_path):
        script = make_script(tmp_path, "toolx", "")
        assert script_shebang_interpreter(script) is None

    def test_unreadable_file_returns_none(self, tmp_path):
        assert script_shebang_interpreter(tmp_path / "does-not-exist") is None


class TestEntrypointReport:
    def test_launcher_off_path_is_a_problem_with_a_fix(self, tmp_path, monkeypatch):
        scripts_dir = tmp_path / "local_bin"
        scripts_dir.mkdir()
        make_script(scripts_dir, "satellome")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [scripts_dir])

        report = check_satellome_entrypoint()

        assert report.ok is False
        problem = "\n".join(report.problems)
        assert str(scripts_dir) in problem
        assert "not on your PATH" in problem
        assert "export PATH=" in problem, "must hand the user the fix, not just the diagnosis"
        assert "python -m satellome" in problem

    def test_two_installs_are_reported_as_shadowing(self, tmp_path, monkeypatch):
        path_dir = tmp_path / "usr_bin"
        own_dir = tmp_path / "venv_bin"
        path_dir.mkdir()
        own_dir.mkdir()
        make_script(path_dir, "satellome", f"#!{sys.executable}\n")
        make_script(own_dir, "satellome", f"#!{sys.executable}\n")
        monkeypatch.setenv("PATH", str(path_dir))
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [own_dir])

        report = check_satellome_entrypoint()

        assert report.ok is False
        problem = "\n".join(report.problems)
        assert str(path_dir) in problem and str(own_dir) in problem
        assert "competing" in problem

    def test_wrong_interpreter_shebang_is_reported(self, tmp_path, monkeypatch):
        path_dir = tmp_path / "usr_bin"
        path_dir.mkdir()
        make_script(path_dir, "satellome", "#!/some/other/python3\n")
        monkeypatch.setenv("PATH", str(path_dir))
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [path_dir])

        report = check_satellome_entrypoint()

        assert report.ok is False
        assert any("/some/other/python3" in p for p in report.problems)

    def test_healthy_install_has_no_problems(self, tmp_path, monkeypatch):
        scripts_dir = tmp_path / "bin"
        scripts_dir.mkdir()
        make_script(scripts_dir, "satellome", f"#!{sys.executable}\n")
        monkeypatch.setenv("PATH", str(scripts_dir))
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [scripts_dir])

        report = check_satellome_entrypoint()

        assert report.problems == []
        assert report.ok is True

    def test_source_checkout_is_a_note_not_a_problem(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PATH", str(tmp_path))
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [tmp_path])

        report = check_satellome_entrypoint()

        assert report.ok is True, "running from a checkout is legitimate"
        assert any("source" in note for note in report.notes)


class TestStartupWarning:
    def test_problem_reaches_the_log(self, tmp_path, monkeypatch, caplog):
        scripts_dir = tmp_path / "local_bin"
        scripts_dir.mkdir()
        make_script(scripts_dir, "satellome")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.delenv(env_check.ENV_DISABLE, raising=False)
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [scripts_dir])

        with caplog.at_level("WARNING"):
            report = warn_if_entrypoint_misconfigured()

        assert report.ok is False
        assert "not on your PATH" in caplog.text
        assert str(scripts_dir) in caplog.text

    def test_opt_out_silences_the_warning_but_not_the_report(self, tmp_path, monkeypatch, caplog):
        scripts_dir = tmp_path / "local_bin"
        scripts_dir.mkdir()
        make_script(scripts_dir, "satellome")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setenv(env_check.ENV_DISABLE, "1")
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [scripts_dir])

        with caplog.at_level("WARNING"):
            report = warn_if_entrypoint_misconfigured()

        assert caplog.text == ""
        assert report.ok is False, "opting out of the message must not fake a healthy install"


class TestDoctorReport:
    def test_report_lists_problem_and_hidden_tool(self, tmp_path, monkeypatch):
        scripts_dir = tmp_path / "local_bin"
        scripts_dir.mkdir()
        make_script(scripts_dir, "satellome")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [scripts_dir])

        report = check_satellome_entrypoint()
        hidden = ToolLocation(
            name="arraysplitter",
            path=str(scripts_dir / "arraysplitter"),
            source="scripts-dir",
            on_path=False,
            dir_off_path=str(scripts_dir),
        )
        missing = ToolLocation(name="trf")
        text = "\n".join(format_doctor_report(report, [hidden, missing]))

        assert "PROBLEMS (2)" in text
        assert "arraysplitter" in text and "off PATH" in text
        assert "trf" in text and "not installed" in text
        assert "NOT on PATH" in text

    def test_healthy_report_says_so(self, tmp_path, monkeypatch):
        scripts_dir = tmp_path / "bin"
        scripts_dir.mkdir()
        make_script(scripts_dir, "satellome", f"#!{sys.executable}\n")
        monkeypatch.setenv("PATH", str(scripts_dir))
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [scripts_dir])

        report = check_satellome_entrypoint()
        found = ToolLocation(name="fastan", path="/opt/bin/fastan", source="managed", on_path=True)
        text = "\n".join(format_doctor_report(report, [found]))

        assert "No problems found." in text
        assert "/opt/bin/fastan" in text
