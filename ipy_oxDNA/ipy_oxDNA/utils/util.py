from pathlib import Path
from typing import Union, Optional
import math


import numpy as np
import matplotlib.pyplot as plt

# general-purpose utility functions

def rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2)
    b, c, d = -axis * math.sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def process_path(p: Union[Path, str], prepend: Union[Path, None] = None) -> Path:
    if isinstance(p, str):
        p = Path(p)
    if not p.is_absolute():
        if p.parts[0] == "~":
            p = p.expanduser()
        elif p.parts[0] not in (".", "..") and prepend is not None:
            p = prepend / p
    return p

# todo: more precision
# source
# oxDNA units
# oxRNA units https://dna.physics.ox.ac.uk/index.php?title=RNA_model_introduction
NEWTONS_PER_UNIT = {
    "rna": 4.93e-11,
    "dna": 4.863e-11
}
METERS_PER_UNIT = {
    "rna": 8.4e-10,
    "dna": 8.518e-10
}
# todo: more units
DEGREES_K_PER_UNIT = {
    "rna": 3000,
    "dna": 3000
}
SECONDS_PER_UNIT = {
    "rna": 3.06e-12,
    "dna": 3.03e-12
}
KG_PER_UNIT = {
    "rna": 5.34e-25,
    "dna": 5.25e-25
}
JOULES_PER_UNIT = {
    "rna": 4.142e-20,
    "dna": 4.142e-20
}

def generate_distinct_colors(n: int) -> np.ndarray:
    """
    Generate distinct colors that are colorblind-friendly by interleaving
    multiple perceptually distinct colormaps.
    This code was written by Claude so there's probably a better way to do it.
    """
    # Use multiple colorblind-safe colormaps with different characteristics
    cmaps = [
        'viridis',  # Purple-blue-green-yellow
        'plasma',  # Purple-pink-orange-yellow
        'cividis',  # Blue-yellow (optimized for colorblindness)
        'coolwarm',  # Blue-red diverging
    ]

    colors = []
    cmap_idx = 0

    # Distribute colors across colormaps, cycling through them
    for i in range(n):
        cmap = plt.get_cmap(cmaps[cmap_idx % len(cmaps)])
        # Sample from the colormap at evenly spaced intervals
        position = (i // len(cmaps)) / max(1, (n // len(cmaps)))
        colors.append(cmap(position))
        cmap_idx += 1

    return np.array(colors)

def si_units(measurement: Union[str, float, np.array], interaction_type: str, measurement_type, to: Optional[str] = None) -> float:
    interaction_type = interaction_type.lower()
    if measurement_type == "T" or measurement_type.lower() == "temperature":
        if to is not None and to == "C":
            return si_units(measurement, interaction_type, "T", "K") - 273.15
        elif to == "K":
            return DEGREES_K_PER_UNIT[interaction_type] * measurement
        else:
            raise ValueError(f"Unsupported temperature conversion to {to}")
    elif measurement_type == "distance" or measurement_type == "d" or measurement_type == "du":
        return METERS_PER_UNIT[interaction_type] * measurement


def ox_units(measurement: Union[str, float], interaction_type: str, units: Optional[str] = None) -> float:
    """
    generic function to convert a measurement to oxDNA or oxRNA units
    for hybrid model, we incorrectly assume units are the same
    :param interaction_type: "dna" or "rna", because the two use slightly different units
    :param units: if measurement is passed as a float or a string with no units, look here
    """
    interaction_type = interaction_type.lower()
    if units is None:
        units = "".join([c for i,c in enumerate(measurement) if not c.isdigit()]).strip()
        if not units:
            raise ValueError(f"Measurement was provided as `{measurement}` but no units were provided!")
        measurement = measurement[:-len(units)]

    measurement = float(measurement)
    # SI Units
    if units == "N":
        return measurement / NEWTONS_PER_UNIT[interaction_type]
    elif units == "s":
        return measurement / SECONDS_PER_UNIT[interaction_type]
    elif units == "K":
        return measurement / DEGREES_K_PER_UNIT[interaction_type]
    elif units == "J":
        return measurement / JOULES_PER_UNIT[interaction_type]
    elif units == "kg":
        return measurement / KG_PER_UNIT[interaction_type]
    elif units == "m":
        return measurement / METERS_PER_UNIT[interaction_type]
    # prefixes (only common ones at scale)
    # force
    elif units == "nN":
        # convert 1000x the measurement from pN
        return ox_units(measurement * 1e-9, interaction_type, "N")
    elif units == "pN":
        # 1 force sim unit s.u  = 48.6 piconewtons units
        return ox_units(measurement * 1e-12, interaction_type, "N")
    # distance
    elif units == "nm":
        return ox_units(measurement * 1e-9, interaction_type, "m")
    elif units == "um" or units == "μm": # allow greek letter mu or latin "u"
        return ox_units(measurement * 1e-6, interaction_type, "m")
    elif units == "pm": # is this even relevant?
        return ox_units(measurement * 1e-12, interaction_type, "m")
    ## temperature
    elif units == "C":
        return ox_units(measurement + 273.15, interaction_type, "K")
    elif units == "F":
        raise ValueError("I said a REAL unit of measurement")
    # time
    elif units == "ns":
        return ox_units(measurement * 1e-9, interaction_type, "s")
    elif units == "ps":
        return ox_units(measurement * 1e-12, interaction_type, "s")
    # energy, good god
    elif units == "nJ":
        return ox_units(measurement * 1e-9, interaction_type, "J")
    elif units == "pJ":
        return ox_units(measurement * 1e-12, interaction_type, "J")
    elif units == "fJ":
        return ox_units(measurement * 1e-15, interaction_type, "J")
    elif units == "aJ":
        return ox_units(measurement * 1e-18, interaction_type, "J")
    elif units == "zJ":
        return ox_units(measurement * 1e-21, interaction_type, "J")

    else:
        raise ValueError(f"Invalid or unsupported unit {units}")