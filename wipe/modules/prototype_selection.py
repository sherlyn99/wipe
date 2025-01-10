import subprocess

from skbio.stats.distance import DistanceMatrix

dm = DistanceMatrix.read("all.dm")
k = 10000
prototype = "prototype()"
