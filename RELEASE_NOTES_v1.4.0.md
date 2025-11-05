# Satellome v1.4.0 - Auto-Install External Tools

## 🎉 Major Release: Simplified Installation

This release significantly improves the installation experience by **automatically installing external tools** during `pip install satellome`.

## ✨ New Features

### Automatic Tool Installation (SAT-37)
- External tools (FasTAN, tanbed, modified TRF) now install automatically with `pip install satellome`
- **Cleaner installation**: Binaries installed to `<site-packages>/satellome/bin/` instead of `~/.satellome/bin/`
- **Graceful failure**: Satellome installs successfully even if tool compilation fails
- **Optional**: Skip with `SATELLOME_SKIP_AUTO_INSTALL=1` environment variable
- **Smart fallback**: Uses `~/.satellome/bin/` if no write permissions to site-packages

### FasTAN and tanbed Installers (SAT-35)
- CLI commands: `--install-fastan`, `--install-tanbed`, `--install-all`
- Automatic compilation from source with dependency checking
- Works on Linux and macOS
- **FasTAN**: Alternative tandem repeat finder by Gene Myers
- **tanbed**: BED format converter from alntools

### Modified TRF Installer (SAT-36)
- CLI command: `--install-trf-large`
- Modified TRF from aglabx for large genomes (chromosomes >2GB)
- Supports large plant and amphibian genomes
- Best support on Linux; manual installation recommended for macOS

## 📦 Installation

### Simple Installation (New!)
```bash
pip install satellome
# That's it! External tools install automatically
```

### Skip Auto-Install
```bash
SATELLOME_SKIP_AUTO_INSTALL=1 pip install satellome
# Install tools later:
satellome --install-all
```

### Development Installation
```bash
git clone https://github.com/aglabx/satellome.git
cd satellome
pip install -e .
```

## 🔧 What Changed

### Installation Process
- **Before**: `pip install satellome` → `satellome --install-all` (manual step)
- **After**: `pip install satellome` (everything included!)

### Binary Location
- **Before**: `~/.satellome/bin/` (pollutes home directory)
- **After**: `<site-packages>/satellome/bin/` (cleaner, removed with pip uninstall)

## 📝 Full Changelog

See [CHANGELOG.md](CHANGELOG.md) for detailed changes.

## 🙏 Credits

- FasTAN by Gene Myers: https://github.com/thegenemyers/FASTAN
- alntools (tanbed) by Richard Durbin: https://github.com/richarddurbin/alntools
- Modified TRF by aglabx: https://github.com/aglabx/trf

## 🐛 Known Issues

- Modified TRF compilation may fail on macOS due to platform-specific code
  - Use manual installation or pre-compiled binaries
  - See README.md for details

## 📊 Installation Requirements

For automatic tool installation:
- git
- make
- C compiler (gcc, clang, or cc)

On macOS: `xcode-select --install`
On Ubuntu/Debian: `sudo apt-get install build-essential git`

---

**Note**: This is a minor version bump (1.3.0 → 1.4.0) as we added significant new functionality while maintaining backward compatibility.
