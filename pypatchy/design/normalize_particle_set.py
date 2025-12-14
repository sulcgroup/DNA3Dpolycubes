import copy
from collections import defaultdict
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt

from matplotlib.patches import Patch

from ipy_oxdna.utils.util import process_path

from ..interaction_matrix import InteractionMatrix
from ..patchy.pl.plparticleset import PLParticleSet
from ..util import selectColor, get_output_dir


import numpy as np

import colorsys

def desaturate(rgb_hex: str, amount: float = 0.5):
    r = int(rgb_hex[1:3], 16) / 255
    g = int(rgb_hex[3:5], 16) / 255
    b = int(rgb_hex[5:7], 16) / 255
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    s = max(0.2, s * (1 - amount))
    l = min(0.95, l + 0.3 * amount)
    r2, g2, b2 = colorsys.hls_to_rgb(h, l, s)
    return (r2, g2, b2)

def normalize_particle_set(particle_set: PLParticleSet,
                           target_strength: float = 6.,
                           max_teeth: int = 4,
                           min_teeth: int = 3,
                           weight_particle=1., weight_strand=0.) -> tuple[InteractionMatrix, dict[int, int]]:
    import numpy as np
    from itertools import product
    from scipy.optimize import differential_evolution

    original_colors = sorted(set(c for c in particle_set.patch_colors() if c > 0))
    num_colors = len(original_colors)

    best_ratio = float("inf")
    best_teeth = None
    best_strengths = None

    for teeth_combo in product(range(min_teeth, max_teeth + 1), repeat=num_colors):
        tooth_counts = np.array(teeth_combo)
        tooth_array = np.array(teeth_combo)
        def objective(strengths: np.ndarray) -> float:
            per_tooth = strengths / tooth_array
            color_ratio = max(per_tooth) / min(per_tooth)

            total_strengths = get_total_strengths(particle_set, strengths, tooth_counts, max_teeth=max_teeth)

            particle_ratio = (max(total_strengths) - min(total_strengths)) / np.mean(total_strengths)

            return (color_ratio * weight_strand) + (particle_ratio * weight_particle)

        bounds = [(0.1, 10.0)] * num_colors

        result = differential_evolution(
            objective,
            bounds,
            strategy="best1bin",
            maxiter=1000,
            polish=True,
            seed=42,
            disp=False
        )

        if result.success:
            s_opt = result.x
            ratio = objective(s_opt)
            if ratio < best_ratio:
                best_ratio = ratio
                best_strengths = s_opt.copy()
                best_teeth = tooth_counts.copy()

    if best_strengths is None or best_teeth is None:
        raise RuntimeError("Joint optimization failed")

    scale_factor = target_strength / get_total_strengths(particle_set, best_strengths, best_teeth).mean()

    best_strengths = scale_factor * best_strengths

    intmat = InteractionMatrix()
    for i, c in enumerate(original_colors):
        intmat[(c, -c)] = best_strengths[i]

    return intmat, {i: int(best_teeth[i]) for i in range(num_colors)}



def get_total_strengths(particle_set: PLParticleSet, int_mat: np.ndarray, tooth_counts: dict[int: int], max_teeth: int=4) -> np.ndarray:
    """
    given a particle set and a candidate interaction matrix (continuous and sorted by color number)
    compute the total patch strength of each particle type in the set
    """
    assert (int_mat > 0).all()
    sorted_colors = list(filter(lambda c: c>0, particle_set.patch_colors()))
    sorted_colors = [c - min(sorted_colors) for c in sorted_colors]

    return np.array([
        sum([
            int_mat[sorted_colors.index(abs(patch.color())-21)] * tooth_counts[sorted_colors.index(abs(patch.color())-21)] / max_teeth
            for patch in particle_type.patches()
        ])
        for particle_type in particle_set.particles()
    ])


def plot_interaction_strengths(
    intmats: dict[str, Union[InteractionMatrix, tuple[InteractionMatrix, dict[int, int]]]],
    fp: Union[Path, None] = None
):

    dataset_labels = list(intmats.keys())

    all_colors = set()
    for val in intmats.values():
        mat = val[0] if isinstance(val, tuple) else val
        all_colors.update(k[0] - 20 for k in mat.intmap().keys() if k[0] == -k[1] and k[0] > 0)
    base_colors = sorted(all_colors)
    color_labels = [str(c) for c in base_colors]
    num_datasets = len(dataset_labels)
    bar_width = 0.8 / num_datasets

    # Layout: left column is tall, right column (ax3) is compact
    fig = plt.figure(figsize=(14, 10))
    ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=3)
    ax2 = plt.subplot2grid((2, 4), (1, 0), colspan=3, sharex=ax1)
    ax3 = plt.subplot2grid((2, 4), (0, 3), rowspan=2)

    dataset_handles = []
    color_handles = []
    ratios = []
    top_colors = []

    for idx, label in enumerate(dataset_labels):
        value = intmats[label]
        if isinstance(value, tuple):
            intmat, tooth_counts = value
        else:
            intmat = value
            tooth_counts = {i: 4 for i in range(len(base_colors))}

        x = [i + idx * bar_width for i in range(len(base_colors))]
        saturation_level = idx / (num_datasets + 1)

        per_tooth_strengths = []
        per_tooth_by_color = {}

        for i, c in enumerate(base_colors):
            key = (c + 20, -(c + 20))
            if key not in intmat.intmap():
                continue
            total_strength = intmat.intmap()[key]
            n_teeth = tooth_counts.get(i, 4)
            per_tooth = total_strength / n_teeth
            per_tooth_strengths.append(per_tooth)
            per_tooth_by_color[i] = per_tooth

            base_hex = selectColor(c - 1)
            bar_color = desaturate(base_hex, saturation_level)

            # ax1: stacked per-tooth blocks
            bottom = 0
            for t in range(n_teeth):
                ax1.bar(x[i], per_tooth, bottom=bottom, width=bar_width,
                        color=bar_color, edgecolor='black', linewidth=0.5)
                bottom += per_tooth
            ax1.text(x[i], bottom, f'{total_strength:.2f}', ha='center', va='bottom', fontsize=8, rotation=90)

            # ax2: per-tooth strength
            ax2.bar(x[i], per_tooth, width=bar_width,
                    color=bar_color, edgecolor='black', linewidth=0.5)
            ax2.text(x[i], per_tooth, f'{n_teeth} teeth', ha='center', va='bottom', fontsize=8, rotation=90)

        if per_tooth_strengths:
            ratio = max(per_tooth_strengths) / min(per_tooth_strengths)
            ratios.append(ratio)
            max_idx = max(per_tooth_by_color, key=lambda i: per_tooth_by_color[i])
            top_colors.append(selectColor(base_colors[max_idx] - 1))
        else:
            ratios.append(0)
            top_colors.append("#999999")

        dataset_handles.append(
            Patch(facecolor=desaturate("#888888", saturation_level), label=label)
        )

    # ax3: max/min ratio bars
    ax3.set_title("Max / Min\nPer-Tooth Strength")
    ax3.set_ylabel("Ratio")
    ax3.set_xticks(range(len(dataset_labels)))
    ax3.set_xticklabels(dataset_labels, rotation=45)
    for i, (r, col) in enumerate(zip(ratios, top_colors)):
        ax3.bar(i, r, color=col, width=0.6, edgecolor='black', linewidth=0.5)
        ax3.text(i, r, f'{r:.2f}', ha='center', va='bottom', fontsize=9, rotation=90)

    # X ticks
    mid_xticks = [i + bar_width * (num_datasets - 1) / 2 for i in range(len(base_colors))]
    ax1.set_xticks(mid_xticks)
    ax1.set_xticklabels(color_labels)
    ax2.set_xticks(mid_xticks)
    ax2.set_xticklabels(color_labels)

    ax1.set_ylabel("Interaction Strength")
    ax1.set_title("Stacked Interaction Strength per Tooth")
    ax2.set_ylabel("Per-Tooth Strength")
    ax2.set_xlabel("Color")

    # Legends
    for c in base_colors:
        color_handles.append(Patch(facecolor=selectColor(c - 1), label=f'Color {c}'))

    fig.legend(handles=color_handles,
               labels=[h.get_label() for h in color_handles],
               title="Patch Color",
               ncol=len(color_handles),
               loc='upper left',
               bbox_to_anchor=(0.03, 1.02),
               fontsize=9)

    fig.legend(handles=dataset_handles,
               labels=[h.get_label() for h in dataset_handles],
               title="Interaction Sets",
               ncol=len(dataset_handles),
               loc='upper right',
               bbox_to_anchor=(0.97, 1.02),
               fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.93])  # leave space for top legends
    if fp is None:
        plt.show()
    else:
        fp = process_path(fp, get_output_dir())
        plt.savefig(fp, format='svg')

def plot_patch_strengths_stacked(
    particle_sets: Union[dict[str, PLParticleSet], PLParticleSet],
    fp: Union[Path, str, None] = None):
    if not isinstance(particle_sets, dict):
        particle_sets = {"Dataset 1": particle_sets}

    hatch_patterns = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']
    all_colors = set()
    all_particle_color_strengths = []

    dataset_labels = list(particle_sets.keys())
    for label in dataset_labels:
        particle_set = particle_sets[label]
        particle_color_strengths = []
        for particle_type in particle_set:
            color_strengths = defaultdict(float)
            for patch in particle_type.patches():
                color = abs(patch.color())
                interaction_strength = particle_set.interaction_matrix()[(color, -color)]
                color_strengths[color - 20] += interaction_strength
                all_colors.add(color - 20)
            particle_color_strengths.append(color_strengths)
        all_particle_color_strengths.append(particle_color_strengths)

    all_colors = sorted(all_colors)
    num_sets = len(all_particle_color_strengths)
    bar_width = 0.8 / num_sets

    plt.figure(figsize=(10, 5))

    dataset_handles = []
    color_handles = []

    for set_index, particle_color_strengths in enumerate(all_particle_color_strengths):
        num_particles = len(particle_color_strengths)
        x = [i + set_index * bar_width for i in range(num_particles)]
        bottoms = [0.0] * num_particles
        totals = [0.0] * num_particles
        hatch = hatch_patterns[set_index % len(hatch_patterns)]
        label = dataset_labels[set_index]

        # Dataset handle for legend
        dataset_handles.append(
            Patch(facecolor='white', edgecolor='black', hatch=hatch, label=label)
        )

        for color in all_colors:
            heights = [
                particle_data.get(color, 0.0)
                for particle_data in particle_color_strengths
            ]
            bar_color = selectColor(color - 1)
            plt.bar(
                x, heights, bottom=bottoms, width=bar_width,
                color=bar_color, hatch=hatch, edgecolor='black',
                label=None
            )
            # Update totals and bottoms for stacking
            totals = [t + h for t, h in zip(totals, heights)]
            bottoms = [b + h for b, h in zip(bottoms, heights)]

        # Add total labels above full bar
        for xi, total in zip(x, totals):
            if total > 0:
                plt.text(
                    xi, total,
                    f'{total:.2f}',
                    ha='center', va='bottom',
                    fontsize=8, rotation=90
                )

    # Color legend (no hatching)
    for color in all_colors:
        color_handles.append(
            Patch(facecolor=selectColor(color - 1), label=f'Color {color}', edgecolor='black')
        )

    plt.xlabel('Particle Index')
    plt.ylabel('Total Patch Strength (per color)')
    plt.title('Patch Strengths per Particle (Stacked by Color and Dataset)')

    legend1 = plt.legend(handles=color_handles, title='Patch Color', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.gca().add_artist(legend1)
    plt.legend(handles=dataset_handles, title='Particle Sets', bbox_to_anchor=(1.05, 0.4), loc='upper left')

    plt.tight_layout()
    if fp is None:
        plt.show()
    else:
        fp = process_path(fp, get_output_dir())
        plt.savefig(fp, format='svg')


def draw_gradient_bar(ax, x, y, width, height, color1, color2):
    from matplotlib.colors import to_rgb
    import numpy as np

    grad = np.linspace(to_rgb(color1), to_rgb(color2), 256).reshape(1, 256, 3)
    extent = (x - width/2, x + width/2, 0, height)
    ax.imshow(grad, aspect='auto', extent=extent, origin='lower', interpolation='bicubic', zorder=0)
    rect = plt.Rectangle((x - width/2, 0), width, height,
                         edgecolor='black', facecolor='none', linewidth=0.5, zorder=1)
    ax.add_patch(rect)

def plot_combined_patch_and_interaction_strengths(
    particle_sets: dict[str, Union[PLParticleSet, tuple[PLParticleSet, dict[int, int]]]],
    fp: Union[Path, None] = None,
):
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import to_rgb
    import numpy as np

    def draw_gradient_bar(ax, x, height, width, color1, color2):
        grad = np.linspace(to_rgb(color1), to_rgb(color2), 256).reshape(1, 256, 3)
        extent = (x - width / 2, x + width / 2, 0, height)
        ax.imshow(grad, aspect='auto', extent=extent, origin='lower', interpolation='bicubic', zorder=0)
        rect = plt.Rectangle((x - width / 2, 0), width, height,
                             edgecolor='black', facecolor='none', linewidth=0.5, zorder=1)
        ax.add_patch(rect)

    dataset_labels = list(particle_sets.keys())
    normalized_particle_sets = dict()
    for label, ps_entry in particle_sets.items():
        if isinstance(ps_entry, tuple):
            ps, teeth_map = ps_entry
            normalized_particle_sets[label] = (ps, defaultdict(lambda: 4, teeth_map))
        else:
            normalized_particle_sets[label] = (ps_entry, defaultdict(lambda: 4))
    particle_sets = normalized_particle_sets
    interaction_sets = {
        label: ps.interaction_matrix()
        for label, (ps, _) in particle_sets.items()
    }

    # Collect all color ids used in interaction matrix (e.g., 1,2,3...)
    all_colors = set()
    for mat in interaction_sets.values():
        all_colors.update(k[0] - 20 for k in mat.intmap().keys() if k[0] == -k[1] and k[0] > 0)
    base_colors = sorted(all_colors)
    color_labels = [str(c) for c in base_colors]
    num_datasets = len(dataset_labels)
    bar_width = 0.8 / num_datasets
    ratio_bar_width = 0.35  # For ax3

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 2, figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])  # Stacked interaction strength
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)  # Per-tooth strength
    ax3 = fig.add_subplot(gs[0, 1])  # Variability
    ax4 = fig.add_subplot(gs[1, 1])  # Patch strengths per particle

    dataset_handles = []
    color_handles = []

    per_tooth_ratios = []
    per_particle_ratios = []
    minmax_color_indices = []
    all_per_tooth_strengths = []
    minmax_particle_indices = []

    for idx, label in enumerate(dataset_labels):

        particle_set, tooth_counts = particle_sets[label]
        intmat = particle_set.interaction_matrix()

        x = [i + idx * bar_width for i in range(len(base_colors))]
        saturation_level = idx / (num_datasets + 1)

        per_tooth_strengths = []
        active_colors = []

        for i, c in enumerate(base_colors):
            key = (c + 20, -(c + 20))
            if key not in intmat.intmap():
                continue
            total_strength = intmat.intmap()[key]
            n_teeth = tooth_counts.get(i, 4)
            per_tooth = total_strength / n_teeth
            per_tooth_strengths.append(per_tooth)
            active_colors.append(c)
            base_hex = selectColor(c - 1)
            bar_color = desaturate(base_hex, saturation_level)

            # ax1: stacked per-tooth blocks
            bottom = 0
            for t in range(n_teeth):
                ax1.bar(x[i], per_tooth, bottom=bottom, width=bar_width,
                        color=bar_color, edgecolor='black', linewidth=0.5)
                bottom += per_tooth
            ax1.text(x[i], bottom, f'{total_strength:.2f}', ha='center', va='bottom', fontsize=8, rotation=90)

            # ax2: per-tooth strength
            ax2.bar(x[i], per_tooth, width=bar_width,
                    color=bar_color, edgecolor='black', linewidth=0.5)
            ax2.text(x[i], per_tooth, f'{n_teeth} teeth', ha='center', va='bottom', fontsize=8, rotation=90)

        all_per_tooth_strengths.append(per_tooth_strengths)

        # Compute total patch strength per particle
        total_strengths = []
        for particle in particle_set:
            total = 0.0
            for patch in particle.patches():
                c = abs(patch.color())
                total += particle_set.interaction_matrix()[(c, -c)]
            total_strengths.append(total)

        # Compute ratios for ax3
        if per_tooth_strengths:
            per_tooth_ratios.append(max(per_tooth_strengths) / min(per_tooth_strengths))
            min_idx = np.argmin(per_tooth_strengths)
            max_idx = np.argmax(per_tooth_strengths)
            minmax_color_indices.append((active_colors[min_idx], active_colors[max_idx]))

        else:
            per_tooth_ratios.append(0)
            minmax_color_indices.append((None, None))

        if total_strengths and min(total_strengths) > 0:
            per_particle_ratios.append(max(total_strengths) / min(total_strengths))
            min_particle_idx = int(np.argmin(total_strengths))
            max_particle_idx = int(np.argmax(total_strengths))
            minmax_particle_indices.append((min_particle_idx, max_particle_idx))
        else:
            per_particle_ratios.append(0)
            minmax_particle_indices.append((None, None))

        dataset_handles.append(
            Patch(facecolor=desaturate("#888888", saturation_level), label=label)
        )

    # ax3: grouped bars with gradient + solid fill
    x_vals = np.arange(len(dataset_labels))
    ax3.set_title("Interaction Strength Variability")
    ax3.set_ylabel("Max / Min Ratio")

    for i, (r, (min_c, max_c)) in enumerate(zip(per_tooth_ratios, minmax_color_indices)):
        if min_c is None or max_c is None:
            continue
        color1 = selectColor(min_c - 1)
        color2 = selectColor(max_c - 1)
        draw_gradient_bar(ax3, x_vals[i] - ratio_bar_width / 2, r, ratio_bar_width, color1, color2)
        ax3.text(x_vals[i] - ratio_bar_width * 2/ 3, r + 0.1, f"{r:.2f}",
                 ha="center", va="bottom", fontsize=7, rotation=90)

        ax3.text(x_vals[i] - ratio_bar_width / 3, r, f"{min_c} vs {max_c}",
                 ha="center", va="bottom", fontsize=7, rotation=90)

    ax3.bar(x_vals + ratio_bar_width / 2, per_particle_ratios, width=ratio_bar_width,
            label="Per-Particle", color="#999999", edgecolor='black')
    for i, (r, (min_i, max_i)) in enumerate(zip(per_particle_ratios, minmax_particle_indices)):
        if min_i is None or max_i is None:
            continue
        ax3.text(x_vals[i] + ratio_bar_width / 3, r, f"{min_i} vs {max_i}",
                 ha="center", va="bottom", fontsize=7, rotation=90)
        ax3.text(x_vals[i] + ratio_bar_width * 2 / 3, r + 0.1, f"{r:.2f}",
                 ha="center", va="bottom", fontsize=7, rotation=90)

    ax3.set_xticks(x_vals)
    ax3.set_xticklabels(dataset_labels, rotation=45)
    ax3.legend()

    # Set X ticks
    mid_xticks = [i + bar_width * (num_datasets - 1) / 2 for i in range(len(base_colors))]
    ax1.set_xticks(mid_xticks)
    ax1.set_xticklabels(color_labels)
    ax2.set_xticks(mid_xticks)
    ax2.set_xticklabels(color_labels)

    ax1.set_ylabel("Interaction Strength")
    ax1.set_title("Stacked Interaction Strength per Tooth")
    ax2.set_ylabel("Per-Tooth Strength")
    ax2.set_xlabel("Color")

    # ax4: patch strength per particle
    all_patch_colors = sorted({
        abs(patch.color()) - 20
        for ps_entry in particle_sets.values()
        for p in ps_entry[0]
        for patch in p.patches()
    })

    for set_index, (label, (particle_set, teeth)) in enumerate(particle_sets.items()):
        saturation_level = set_index / (num_datasets + 1)
        patch_strengths = []
        for particle in particle_set:
            color_strengths = defaultdict(float)
            for patch in particle.patches():
                c = abs(patch.color())
                strength = particle_set.interaction_matrix()[(c, -c)]
                color_strengths[c - 20] += strength
            patch_strengths.append(color_strengths)

        x = [i + set_index * bar_width for i in range(len(patch_strengths))]
        bottoms = [0.0] * len(patch_strengths)
        totals = [0.0] * len(patch_strengths)

        for color in all_patch_colors:
            heights = [data.get(color, 0.0) for data in patch_strengths]
            base_hex = selectColor(color - 1)
            bar_color = desaturate(base_hex, saturation_level)

            ax4.bar(x, heights, bottom=bottoms, width=bar_width,
                    color=bar_color, edgecolor='black')
            bottoms = [b + h for b, h in zip(bottoms, heights)]
            totals = [t + h for t, h in zip(totals, heights)]

        for xi, total in zip(x, totals):
            if total > 0:
                ax4.text(xi, total, f'{total:.2f}', ha='center', va='bottom', fontsize=8, rotation=90)

    ax4.set_title('Patch Strength per Particle')
    ax4.set_xlabel('Particle Index')
    ax4.set_ylabel('Total Patch Strength')

    # Legends
    for c in base_colors:
        color_handles.append(Patch(facecolor=selectColor(c - 1), label=f'Color {c}'))

    fig.legend(handles=color_handles,
               title="Patch Color", ncol=len(color_handles),
               loc='upper center', bbox_to_anchor=(0.5, 1.02), fontsize=9)

    fig.legend(handles=dataset_handles,
               title="Datasets", ncol=len(dataset_handles),
               loc='lower center', bbox_to_anchor=(0.5, -0.01), fontsize=9)

    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    if fp is None:
        plt.show()
    else:
        fp = process_path(fp, get_output_dir())
        plt.savefig(fp, format='svg')
