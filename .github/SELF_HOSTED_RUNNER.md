# Self-hosted GitHub Actions runner (77.234.216.98)

Why: the GitHub-hosted Actions quota for this account is exhausted. Both
workflows (`tests.yml`, `publish.yml`) select the runner from a repository
variable, the same pattern AgentMux uses:

```yaml
runs-on: ${{ (vars.CI_RUNNER == 'self-hosted' && github.event.pull_request.head.repo.fork != true)
             && fromJSON('["self-hosted","vm"]') || 'ubuntu-latest' }}
```

* `CI_RUNNER = self-hosted` → jobs go to a runner labelled `self-hosted` + `vm`.
* variable unset or anything else → back to `ubuntu-latest`, no workflow edit.
* a pull request **from a fork** always goes to a GitHub-hosted runner.

## Fork PRs never run here — on purpose

`aglabx/satellome` is a **public** repository. A self-hosted runner executing a
pull request from a fork runs a stranger's code on our own machine, as the
runner user, with our network. That is why the expression above excludes forks
rather than trusting the runner sandbox — there is none; the runner is a shell
on the box.

Also set, once, in **Settings → Actions → General → Fork pull request workflows
from outside collaborators**: *Require approval for all outside collaborators*.

## 1. Install the runner on 77.234.216.98

Run **on the server**, as a non-root user that does not own anything valuable:

```bash
sudo adduser --disabled-password --gecos "" ghrunner
sudo -iu ghrunner

mkdir -p ~/actions-runner && cd ~/actions-runner
curl -sSLO https://github.com/actions/runner/releases/download/v2.330.0/actions-runner-linux-x64-2.330.0.tar.gz
tar xzf actions-runner-linux-x64-2.330.0.tar.gz
```

Get a registration token (expires in an hour) from
**Settings → Actions → Runners → New self-hosted runner**, then:

```bash
./config.sh --url https://github.com/aglabx/satellome \
            --token <REGISTRATION_TOKEN> \
            --name satellome-vm \
            --labels vm \
            --unattended --replace
```

`self-hosted` is added by GitHub automatically; `vm` is the label the workflows
ask for, so it must be spelled exactly.

Install it as a service so it survives reboots:

```bash
exit                       # back to a sudo-capable user
cd /home/ghrunner/actions-runner
sudo ./svc.sh install ghrunner
sudo ./svc.sh start
sudo ./svc.sh status
```

## 2. What the box needs

| Requirement | Used by | Notes |
|---|---|---|
| `git`, `curl`, `tar` | every job | |
| build toolchain (`build-essential`) | `actions/setup-python` fallback builds | |
| Python 3.9 / 3.10 / 3.11 | the test matrix | `setup-python` downloads them into `~/actions-runner/_work/_tool` on first run; nothing to preinstall |
| Docker | **publish job only** | `pypa/gh-action-pypi-publish` is a container action and will fail without it |

```bash
sudo apt-get update
sudo apt-get install -y git curl tar build-essential docker.io
sudo usermod -aG docker ghrunner   # then restart the runner service
```

Disk: the tool cache plus a pip cache is a few GB. `~/actions-runner/_work` is
reused between runs — that is what makes it fast, and also why a job must never
assume a clean machine beyond its own checkout.

## 3. Turn it on

**Settings → Secrets and variables → Actions → Variables → New variable**

```
CI_RUNNER = self-hosted
```

Delete the variable (or set it to anything else) to go straight back to
GitHub-hosted runners. No commit needed either way.

## 4. Check it

Push any commit, or run the Tests workflow from the Actions tab, and confirm
the job header names `satellome-vm` instead of `ubuntu-latest`.

Trusted publishing to PyPI keeps working: the OIDC token is minted by GitHub,
not by the runner, so `id-token: write` behaves the same here.
