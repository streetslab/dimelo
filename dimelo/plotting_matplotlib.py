from __future__ import annotations

from pathlib import Path


def save_figure(
    fig,
    path,
    *,
    dpi: int = 300,
    bbox_inches: str = "tight",
    transparent: bool = False,
) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches=bbox_inches,
        transparent=transparent,
    )
    return output_path
