from os.path import join, dirname
import numpy as np
import pandas as pd

path_here = dirname(dirname(__file__))


def load_file(name):
    """ Return a requested data file. """
    return pd.read_csv(join(path_here, "" + name + ".csv"), delimiter=",", comment="#")

data_a = load_file("fig3a")
data_b = load_file("fig3b")