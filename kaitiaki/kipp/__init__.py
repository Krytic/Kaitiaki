"""kipp: shaded Kippenhahn diagrams from Cambridge STARS plot files."""
from .io import load_plot
from .decode import decode_row, decode_all
from .rasterise import rasterise, time_edges
from .render import plot_kippenhahn
from .cli import main

__all__ = [
    "load_plot",
    "decode_row",
    "decode_all",
    "rasterise",
    "time_edges",
    "plot_kippenhahn",
]

__version__ = "0.1.0"
