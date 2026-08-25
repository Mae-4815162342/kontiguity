# data managment
from collections import defaultdict
import pandas as pd
import numpy as np

# ignoring warnings for median and mean in numpy (expected for Hi-C data)
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN slice encountered")
warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")

# system
import subprocess
import datetime
import shutil
import sys
import os

# multitasking
from concurrent.futures import ProcessPoolExecutor
from queue import Queue, Empty
import threading

# network
import requests
import json

# mathematical operations
from itertools import permutations
from functools import reduce
import math

# biological data
# import cooltools # if re-integrated, add "cooltools>=0.7.1" to dependencies
# import bioframe # if re-integrated, add "bioframe" to dependencies
import cooler

# display
import matplotlib.patches as patches
import matplotlib.gridspec as grid
import matplotlib.pyplot as plt