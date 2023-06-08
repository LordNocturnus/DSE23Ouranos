from astropy import units as u
from matplotlib import pyplot as plt
import numpy as np
from typing import List, Tuple

from poliastro.bodies import Venus, Earth, Mars, Jupiter, Uranus, Saturn
from poliastro.plotting.tisserand import TisserandPlotter, TisserandKind
from poliastro.plotting._base import BODY_COLORS


planets = [Venus, Earth, Mars, Jupiter]
def run_tisserand(planets):
    # Build custom axis
    fig, ax = plt.subplots(1, 1, figsize=(15, 7))
    ax.set_title("Energy Tisserand for Venus, Earth, Mars, Saturn")
    ax.set_xlabel("$R_{p} [AU]$")
    ax.set_ylabel("Heliocentric Energy [km2 / s2]")
    ax.set_xscale("log")
    ax.set_xlim(10**-0.4, 10**0.15)
    ax.set_ylim(-1500, 0)

    tp = TisserandPlotter(axes=ax, kind=TisserandKind.ENERGY)

    v_inf_launch = 3
    for planet in planets:
        ax = tp.plot(planet, (1, 14) * u.km / u.s, num_contours=14
                     )

    # Let us label previous figure
    tp.ax.text(0.70, -650, "Venus", color=BODY_COLORS["Venus"])
    tp.ax.text(0.95, -500, "Earth", color=BODY_COLORS["Earth"])
    tp.ax.text(1.35, -350, "Mars", color=BODY_COLORS["Mars"])
    tp.ax.text(2, -200, "Jupiter", color=BODY_COLORS["Jupiter"])
    tp.ax.text(3, -190, "Uranus", color=BODY_COLORS["Uranus"])

    plt.show()