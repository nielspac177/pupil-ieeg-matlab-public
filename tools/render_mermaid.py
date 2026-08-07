"""
Render the Mermaid diagram in docs/methods_diagram.md to a publication PNG.

The `.md` is the source of truth: it diffs cleanly, renders natively on GitHub,
and is what anyone edits. Journals cannot embed it, so this script extracts the
fenced ```mermaid block and rasterises it at print resolution for the
manuscript. Running it is the only supported way to regenerate that figure —
the PNG is a build product, not something to touch by hand.

Requires Node (mermaid-cli is fetched by npx on first run).

Run:
    python tools/render_mermaid.py
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SOURCE = REPO / "docs" / "methods_diagram.md"
OUTPUT = REPO / "results" / "figures" / "Fig0b_methods_pipeline.png"

# 4x scale on a 1400 px canvas lands near 600 dpi at a two-column width.
SCALE = "4"
WIDTH = "1400"


# Emoji are part of the repository style and render well on GitHub, but a
# journal methods figure carrying them reads as informal. The print render
# strips them; the .md that humans read keeps them.
EMOJI = re.compile(
    r"[\U0001F300-\U0001FAFF\u2190-\u21FF\u2600-\u27BF\uFE0F\u2B00-\u2BFF]+\s*")


def strip_emoji(diagram: str) -> str:
    return EMOJI.sub("", diagram)


def extract_diagram(path: Path) -> str:
    text = path.read_text()
    blocks = re.findall(r"```mermaid\n(.*?)```", text, flags=re.S)
    if not blocks:
        raise ValueError(f"no ```mermaid block in {path}")
    if len(blocks) > 1:
        raise ValueError(
            f"{path} has {len(blocks)} mermaid blocks; this renderer expects the "
            "file to contain exactly one, so that the manuscript figure and the "
            "document cannot disagree about which diagram is the pipeline")
    return blocks[0]


def main() -> int:
    plain = "--with-emoji" not in sys.argv
    if shutil.which("npx") is None:
        print("npx not found: install Node.js to render the diagram",
              file=sys.stderr)
        return 1

    diagram = extract_diagram(SOURCE)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as workdir:
        source = Path(workdir) / "diagram.mmd"
        source.write_text(strip_emoji(diagram) if plain else diagram)
        # A white background rather than the default transparent one: a
        # transparent PNG placed in a Word document renders on whatever the
        # viewer's page colour happens to be.
        command = [
            "npx", "-y", "-p", "@mermaid-js/mermaid-cli", "mmdc",
            "-i", str(source), "-o", str(OUTPUT),
            "-b", "white", "-s", SCALE, "-w", WIDTH,
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0 or not OUTPUT.exists():
            print(result.stdout[-2000:], file=sys.stderr)
            print(result.stderr[-2000:], file=sys.stderr)
            return 1

    print(f"Wrote {OUTPUT} ({OUTPUT.stat().st_size/1e3:.0f} kB)"
          f"{'' if plain else ' with emoji'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
