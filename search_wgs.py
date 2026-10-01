"""Convenience launcher for the blast-wgs command."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from blast_api.cli import main
if __name__ == "__main__":
    raise SystemExit(main())
