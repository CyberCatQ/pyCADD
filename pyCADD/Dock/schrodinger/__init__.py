try:
    import schrodinger
except ImportError:
    raise ValueError("Current environment does not include Schrodinger package. \nPlease activate Schrodinger python environment first.")

from schrodinger import structure as struc
from schrodinger.structure import StructureReader, StructureWriter, Structure
from schrodinger.job import jobcontrol as jc
from schrodinger.application.glide import poseviewconvert as pvc