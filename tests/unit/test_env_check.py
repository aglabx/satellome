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


class TestEnsurePathEntry:
    def test_appends_marked_block_once(self, tmp_path, monkeypatch):
        rc = tmp_path / ".bashrc"
        rc.write_text("# existing content\n")
        target = tmp_path / "local_bin"
        target.mkdir()
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setenv("SHELL", "/bin/bash")

        first = env_check.ensure_path_entry(target, rc_file=rc)

        assert first.status == "added"
        assert first.changed is True
        body = rc.read_text()
        assert "# existing content" in body, "must append, never rewrite the file"
        assert env_check.PATH_BLOCK_START in body and env_check.PATH_BLOCK_END in body
        assert f'export PATH="{target}:$PATH"' in body

        second = env_check.ensure_path_entry(target, rc_file=rc)

        assert second.status == "already-configured"
        assert second.changed is False
        assert body.count(env_check.PATH_BLOCK_START) == 1
        assert rc.read_text() == body, "second call must not touch the file"

    def test_does_not_duplicate_a_hand_written_entry(self, tmp_path, monkeypatch):
        rc = tmp_path / ".bashrc"
        target = tmp_path / "local_bin"
        target.mkdir()
        rc.write_text(f'export PATH="{target}:$PATH"\n')
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setenv("SHELL", "/bin/bash")

        result = env_check.ensure_path_entry(target, rc_file=rc)

        assert result.status == "already-configured"
        assert env_check.PATH_BLOCK_START not in rc.read_text()

    def test_commented_out_entry_does_not_count_as_configured(self, tmp_path, monkeypatch):
        rc = tmp_path / ".bashrc"
        target = tmp_path / "local_bin"
        target.mkdir()
        rc.write_text(f'# export PATH="{target}:$PATH"\n')
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setenv("SHELL", "/bin/bash")

        assert env_check.ensure_path_entry(target, rc_file=rc).status == "added"

    def test_directory_already_on_path_is_a_no_op(self, tmp_path, monkeypatch):
        rc = tmp_path / ".bashrc"
        rc.write_text("")
        monkeypatch.setenv("PATH", str(tmp_path))
        monkeypatch.setenv("SHELL", "/bin/bash")

        result = env_check.ensure_path_entry(tmp_path, rc_file=rc)

        assert result.status == "already-on-path"
        assert rc.read_text() == "", "must not write when there is nothing to fix"

    def test_unwritable_startup_file_is_reported_not_swallowed(self, tmp_path, monkeypatch):
        rc = tmp_path / ".bashrc"
        rc.write_text("")
        rc.chmod(0o400)
        target = tmp_path / "local_bin"
        target.mkdir()
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setenv("SHELL", "/bin/bash")

        result = env_check.ensure_path_entry(target, rc_file=rc)

        assert result.status == "error"
        assert result.ok is False
        assert str(rc) in result.detail
        rc.chmod(0o600)

    def test_fish_uses_fish_add_path(self, tmp_path, monkeypatch):
        rc = tmp_path / "config.fish"
        rc.write_text("")
        target = tmp_path / "local_bin"
        target.mkdir()
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setenv("SHELL", "/usr/bin/fish")

        env_check.ensure_path_entry(target, rc_file=rc)

        assert f"fish_add_path {target}" in rc.read_text()
        assert "export PATH" not in rc.read_text()

    def test_home_paths_are_written_portably(self, tmp_path, monkeypatch):
        home = tmp_path / "home"
        (home / ".local" / "bin").mkdir(parents=True)
        rc = home / ".bashrc"
        rc.write_text("")
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setenv("SHELL", "/bin/bash")

        env_check.ensure_path_entry(home / ".local" / "bin", rc_file=rc)

        assert 'export PATH="$HOME/.local/bin:$PATH"' in rc.read_text(), (
            "a literal home path breaks if the file is shared between machines"
        )

    def test_activate_command_is_offered_for_the_current_shell(self, tmp_path, monkeypatch):
        rc = tmp_path / ".bashrc"
        rc.write_text("")
        target = tmp_path / "local_bin"
        target.mkdir()
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setenv("SHELL", "/bin/bash")

        result = env_check.ensure_path_entry(target, rc_file=rc)

        assert result.activate_command == f"source {rc}", (
            "a child process cannot change the parent shell's PATH; the user "
            "must be told what to run"
        )


class TestShellRcFiles:
    def test_zsh_respects_zdotdir(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SHELL", "/bin/zsh")
        monkeypatch.setenv("ZDOTDIR", str(tmp_path))
        assert env_check.shell_rc_files()[0] == tmp_path / ".zshrc"

    def test_unknown_shell_falls_back_to_profile(self, monkeypatch):
        monkeypatch.setenv("SHELL", "/bin/some-exotic-shell")
        assert env_check.shell_rc_files()[0].name == ".profile"

    def test_bash_target_is_first(self, monkeypatch):
        monkeypatch.setenv("SHELL", "/bin/bash")
        names = [p.name for p in env_check.shell_rc_files()]
        assert names[0] in (".bashrc", ".bash_profile")
        assert ".profile" in names


class TestAutoFixOnStartup:
    def _broken_install(self, tmp_path, monkeypatch):
        scripts_dir = tmp_path / "local_bin"
        scripts_dir.mkdir()
        make_script(scripts_dir, "satellome")
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("SHELL", "/bin/bash")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.delenv(env_check.ENV_DISABLE, raising=False)
        monkeypatch.delenv(env_check.ENV_NO_FIX, raising=False)
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [scripts_dir])
        # Which startup file bash reads is platform-dependent; follow the same
        # resolution the code uses instead of hardcoding one.
        rc = env_check.shell_rc_files()[0]
        rc.write_text("")
        return scripts_dir, rc

    def test_startup_check_repairs_and_says_so(self, tmp_path, monkeypatch, caplog):
        scripts_dir, rc = self._broken_install(tmp_path, monkeypatch)

        with caplog.at_level("WARNING"):
            report = env_check.warn_if_entrypoint_misconfigured()

        assert report.path_fix is not None and report.path_fix.changed
        assert env_check.PATH_BLOCK_START in rc.read_text()
        # written in $HOME form, see test_home_paths_are_written_portably
        assert f'export PATH="$HOME/{scripts_dir.name}:$PATH"' in rc.read_text()
        assert "Added" in caplog.text
        assert "source" in caplog.text, "must tell the user how to affect this shell"

    def test_no_path_fix_opt_out_warns_without_writing(self, tmp_path, monkeypatch, caplog):
        scripts_dir, rc = self._broken_install(tmp_path, monkeypatch)
        monkeypatch.setenv(env_check.ENV_NO_FIX, "1")

        with caplog.at_level("WARNING"):
            report = env_check.warn_if_entrypoint_misconfigured()

        assert rc.read_text() == "", "opt-out must leave the file untouched"
        assert report.path_fix is None
        assert "not on your PATH" in caplog.text, "the problem is still reported"

    def test_fix_is_not_attempted_when_check_is_disabled(self, tmp_path, monkeypatch):
        _, rc = self._broken_install(tmp_path, monkeypatch)
        monkeypatch.setenv(env_check.ENV_DISABLE, "1")

        env_check.warn_if_entrypoint_misconfigured()

        assert rc.read_text() == ""

    def test_second_run_does_not_append_again(self, tmp_path, monkeypatch):
        _, rc = self._broken_install(tmp_path, monkeypatch)

        env_check.warn_if_entrypoint_misconfigured()
        after_first = rc.read_text()
        env_check.warn_if_entrypoint_misconfigured()

        assert rc.read_text() == after_first
        assert after_first.count(env_check.PATH_BLOCK_START) == 1

    def test_shadowed_launcher_is_not_papered_over_with_a_path_edit(self, tmp_path, monkeypatch, caplog):
        """Two competing installs need one removed; adding a directory cannot help."""
        path_dir = tmp_path / "usr_bin"
        own_dir = tmp_path / "venv_bin"
        path_dir.mkdir()
        own_dir.mkdir()
        make_script(path_dir, "satellome", f"#!{sys.executable}\n")
        make_script(own_dir, "satellome", f"#!{sys.executable}\n")
        rc = tmp_path / ".bashrc"
        rc.write_text("")
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("SHELL", "/bin/bash")
        monkeypatch.setenv("PATH", str(path_dir))
        monkeypatch.delenv(env_check.ENV_DISABLE, raising=False)
        monkeypatch.delenv(env_check.ENV_NO_FIX, raising=False)
        monkeypatch.setattr(env_check, "interpreter_script_dirs", lambda: [own_dir])

        with caplog.at_level("WARNING"):
            report = env_check.warn_if_entrypoint_misconfigured()

        assert report.path_fix is None
        assert rc.read_text() == ""
        assert "competing" in caplog.text


class TestStartupFileSelection:
    def test_prefers_an_existing_file_over_creating_a_new_one(self, tmp_path, monkeypatch):
        """Creating ~/.bash_profile can stop a login shell reading ~/.profile."""
        home = tmp_path / "home"
        home.mkdir()
        bashrc = home / ".bashrc"
        bashrc.write_text("# mine\n")
        target_dir = home / ".local" / "bin"
        target_dir.mkdir(parents=True)
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setenv("SHELL", "/bin/bash")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setattr(
            env_check, "shell_rc_files",
            lambda: [home / ".bash_profile", bashrc, home / ".profile"],
        )

        result = env_check.ensure_path_entry(target_dir)

        assert result.rc_file == str(bashrc)
        assert env_check.PATH_BLOCK_START in bashrc.read_text()
        assert not (home / ".bash_profile").exists(), "must not create a new startup file"

    def test_creates_the_preferred_file_when_none_exist(self, tmp_path, monkeypatch):
        home = tmp_path / "home"
        home.mkdir()
        target_dir = home / ".local" / "bin"
        target_dir.mkdir(parents=True)
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setenv("SHELL", "/bin/bash")
        monkeypatch.setenv("PATH", "/nonexistent-dir-for-test")
        monkeypatch.setattr(env_check, "shell_rc_files", lambda: [home / ".bashrc"])

        result = env_check.ensure_path_entry(target_dir)

        assert result.status == "added"
        assert (home / ".bashrc").exists()
