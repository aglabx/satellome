try:
    from importlib.metadata import version
    __version__ = version("satellome")
except Exception:
    __version__ = "unknown"