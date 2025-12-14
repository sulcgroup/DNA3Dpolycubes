from __future__ import annotations

import itertools
import json
import os
import pickle
from pathlib import Path
from typing import Union

import libtlm
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from .. import tlm_data
from ..polycubesRule import PolycubesRule
from ...design.solve_utils import toPolycube
from ...util import get_output_dir

NUMBER_WHICH_IS_LAST = -2 # todo fix this

class TLMDataGroup:
    __raws: list[list[libtlm.TLMHistoryRecord]]
    rule: PolycubesRule
    T: float
    totalCubeCount: int
    data: pd.DataFrame

    num_sims = property(lambda self: len(self.__raws))

    def __init__(self, raws: Union[list[list[libtlm.TLMHistoryRecord]], Path, str], T: Union[float, None] = None,
                 rule: Union[PolycubesRule, None] = None, totalCubes: Union[int, None] = None):
        # if raw data is provided as a pickle
        if T is not None:
            self.T = T
        if rule is not None:
            self.rule = rule
        if totalCubes is not None:
            self.totalCubeCount = totalCubes
        if isinstance(raws, (str, Path)):
            if type(raws) == str:
                raws = Path(raws)
            if raws.is_dir():
                # if i gave it a directory with json files, god help me
                self.__raws = []
                json_files = list(raws.glob("sim*.json"))
                assert len(json_files) > 0, f"Directory {str(raws)} does not contain any files with name pattern sim*.json"
                for fp in json_files:
                    with fp.open('r') as file:
                        json_data = json.load(file)
                        if isinstance(json_data, dict):
                            settings = json_data["settings"]

                            if T is None: # and implicitly the rest
                                self.T = settings["T"]
                                self.rule = PolycubesRule(rule_json=settings["cube_types"])
                                self.totalCubeCount = sum([sum(stage["cube_type_counts"]) for stage in settings["stages"]])

                            else:
                                assert self.T - settings["T"] < 1e-5
                                assert self.rule == PolycubesRule(rule_json=settings["cube_types"])
                                assert self.totalCubeCount == sum([sum(stage["cube_type_counts"]) for stage in settings["stages"]])
                            records = [tlm_data.TLMHistoryRecord(**record) for record in json_data["records"]]
                        else:
                            assert totalCubes is not None and T is not None and rule is not None,\
                                "Must provide either records json with tlm " \
                                "settings or provide setting values in constructor"
                            records = [tlm_data.TLMHistoryRecord(**record) for record in json_data]
                        self.__raws.append(records)
            elif raws.is_file():

                # raise Exception(f"Invalid file format for records. file name: {raws}")
                # if i have pickled data
                with raws.open("rb") as f:
                    self.__raws = pickle.load(f)
            else:
                raise Exception(f"No directory/path {str(raws)}")

            # todo: type assertions
        else:
            pass
        self.data = pd.DataFrame(itertools.chain.from_iterable([
            [{
                "sim": i,
                "tidx": tidx,
                "step": datapoint.stepCount(),
                "energy": datapoint.energy(),
                "nPCs": datapoint.numPolycubes(),
                "nCs": datapoint.numCubeInstances(),
            } for tidx, datapoint in enumerate(sim_records) if datapoint.energy() is not None]
            for i, sim_records in enumerate(self.__raws)
        ]))
        self.data.dropna(inplace=True)  # TODO: make not problem
        assert self.rule != PolycubesRule()

    def polycube_size_histograms(self, timepoint: int = NUMBER_WHICH_IS_LAST, log_scale: bool = False, ax: plt.Axes = None,
                                 title: str = None) -> plt.Figure:
        # Extract polycube sizes from the raw data
        polycube_sizes = [[pc.numCubes() for pc in record[timepoint].getPolycubes()] for record in self.__raws]

        # Create a new figure and axes if not provided
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = ax.get_figure()

        # Define a color map for distinct groups
        colors = plt.cm.tab10.colors  # Use a colormap with up to 10 distinct colors

        # Compute the range of data
        all_sizes = [size for sizes in polycube_sizes for size in sizes]
        bins = np.linspace(min(all_sizes), max(all_sizes), 21)  # Create 20 bins

        # Plot histograms for each group separately with slight offsets
        bar_width = (bins[1] - bins[0]) / (len(polycube_sizes) + 1)  # Adjust bar width based on the number of records
        for i, sizes in enumerate(polycube_sizes):
            bin_heights, _ = np.histogram(sizes, bins=bins)
            bin_centers = bins[:-1] + bar_width * (i + 0.5)  # Offset each bar group slightly
            ax.bar(
                bin_centers,
                bin_heights,
                width=bar_width,
                color=colors[i % len(colors)],
                edgecolor='black',
                alpha=0.7,
                label=f'Record {i + 1}'
            )

        # Set log scale if requested
        if log_scale:
            ax.set_yscale('log')

        # Set titles and labels
        default_title = 'Polycube Size Distribution by Record'
        ax.set_title(title if title else default_title)
        ax.set_xlabel('Number of Cubes')
        ax.set_ylabel('Frequency (log scale)' if log_scale else 'Frequency')
        ax.grid(axis='y', linestyle='--', alpha=0.7)

        # Add a legend
        ax.legend(title="Records", loc='upper right')

        return fig

    def raw_data(self) -> list[list[libtlm.TLMHistoryRecord]]:
        return self.__raws

    def sim_raw_data(self, i: int) -> list[libtlm.TLMHistoryRecord]:
        return self.__raws[i]

    def get_slice(self, start: Union[int, None] = None, stop: Union[int, None] = None):
        assert start is not None or stop is not None

    def last_records(self) -> list[libtlm.TLMHistoryRecord]:
        # for reasons i don't undertand and am afraid to fix record[-1] is bad, but this should work too
        return [record[NUMBER_WHICH_IS_LAST] for record in self.raw_data()]

    def __len__(self) -> int:
        assert all(len(x) == len(self.__raws[0]) for x in self.__raws[1:]), "Replicas have inconsistant lengths!"
        return len(self.__raws[0])

    def save_data(self, filename: str):
        """
        save data about the simulation as a csv
        """
        if filename.startswith("~"):
            filename = os.path.expanduser(filename)
        if not filename.startswith("/"):
            filename = str(get_output_dir() / "tlm" / filename)
        self.data.to_csv(filename, sep="\t", index=False)

    def save_raws(self, filename: str):
        """
        save the data themselves as pickle
        """
        if filename.startswith("~"):
            filename = os.path.expanduser(filename)
        if not filename.startswith("/"):
            filename = str(get_output_dir() / "tlm" / filename)
        with open(filename, "wb") as f:
            pickle.dump(self.__raws, f)

    def save_polycubes(self, directory: str, timepoint: int = NUMBER_WHICH_IS_LAST, sim: int = None, istart: int=1):
        """
        save the polycubes that make up this record as json files
        """
        if sim is None:
            i = istart
            for sim in range(self.num_sims):
                self.save_polycubes(directory, timepoint, sim, i)
                i += len(self.sim_raw_data(sim))
        else:
            record = self.sim_raw_data(sim)[timepoint]
            for i, polycube in enumerate(record.getPolycubes()):
                toPolycube(self.rule, polycube).save_polycube(directory + f"/pc{i+istart}.json")