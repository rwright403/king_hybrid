from pathlib import Path
script_path = Path(__file__, '..').resolve()

with open(script_path.joinpath("program/t.txt"))