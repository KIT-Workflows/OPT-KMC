'''required packages from python lib'''
import itertools
from math import exp 
from random import choice

'''                              SITE CLASS                                 '''
###############################################################################
# every site on lattice is an object of this class with 3 attributes for coordinate
# type of species, and species itself
class Site:
    def __init__(self, coordinate, species, bonds, nbr, status):
        self.coordinate = coordinate
        self.species = species
        self.nbr = nbr
        self.bonds = bonds
        self.status = status

    def setN(self, coordinate):
        self.coordinate = coordinate

    def setS(self, species):
        self.species = species

    def setnb(self, nbr):
        self.nbr = nbr

    def setb(self, bonds):
        self.bonds = bonds

    def setss(self, status):
        self.status = status

    def getC(self, coordinate):
        return self.coordinate

    def getS(self, species):
        return self.species

    def getnb(self, nbr):
        return self.nbr

    def getb(self, bonds):
        return self.bonds

    def getss(self, status):
        return self.status
    
'''                                 EVENT CLASS                             '''
###############################################################################
# events/reaction are objects of this class with four attributes
class Event:
    # Generic attributes
    def __init__(self, reactant, product, barrier, rate):
        self.reactant = reactant
        self.rate = rate
        self.product = product
        self.barrier = barrier

    # -------------------------------------------------------------------------------
    # assiging parameters in class

    def set_reactant(self, reactant):
        self.reactant = reactant

    def set_rate(self, rate):
        self.rate = rate

    def set_product(self, product):
        self.product = product

    def set_statt(self, barrier):
        self.barrier = barrier
    # ------------------------------------------------------------------------------
    # Get parameters

    def get_reactant(self, reactant):
        return self.reactant

    def get_rate(self, rate):
        return self.rate

    def get_product(self, product):
        return self.product

    def get_stat(self, stat):
        return self.stat
# -----------------------------NEW EVENT SUBCLASS-------------------------------
"""
this subclass helps in collecting events
"""
class New_event(Event):
    def __init__(self, stat, reactant, product, rate, new_coord):
        super().__init__(stat, reactant, product, rate)
        self.new_coord = new_coord

    def set_newC(self, new_coord):
        self.new_coord = new_coord

    def get_newC(self, new_coord):
        return self.new_coord


'''                             FUNCTIONS                                   '''
###############################################################################
# Boltzman anf Planck constants
K_B = 8.617333262145e-5  # EVK-1
h = 4.135667696e-15  # EV.S
R = 8.314
###############################################################################
def bar_rate(e, T):
    return ((K_B * T) / h) * exp(-e / (K_B * T))
###############################################################################
# prepare float values out of tracking dictionary
def round_float_list(float_list, decimal_points):
    float_list = [round(float(item), decimal_points) for item in float_list]
    return float_list
# ----------------------------------- lattice-----------------------------------
"""
this function provides a NxN lattice with all grid points as site object which 
were defined before.

generating a square lattice:
x : 1 to xdim +1 to cover all xdim
y = 0 to ydim +2 as we need y=0 for the Electrode and Y=ydim+1 for the last layer
"""
def Lattice(xdim, ydim):
    # lattice points as a list of 2-elements sublists as [x,y]
    l = [[i, j] for i, j in itertools.product(range(1, xdim + 1), range(0, ydim + 1))]
    # initialize lattice as a list of site objects// coordinate,species,bonds,nbr,status
    sites = [Site(i, [], [], [[i[0] + 1, i[1]], [i[0] - 1, i[1]], [i[0], i[1] + 1], [i[0], i[1] - 1]], [[], [], '', [0,'']])for i in l]
    # predefining the neighbors
    for i in sites:
        temp = i.nbr
        temp2 = [i for i in sites if i.coordinate in temp]
        i.nbr = temp2
    #PBC from left and right walls
    for i in sites:
        if i.coordinate[0] == 1:
           new = [j for j in sites if j.coordinate == [xdim, i.coordinate[1]]]
           i.nbr = i.nbr + new
        if i.coordinate[0] == xdim:
           new = [j for j in sites if j.coordinate == [1, i.coordinate[1]]]
           i.nbr = i.nbr + new
    return sites
###############################################################################
# ------------------------------selection function-----------------------------
def ratee(i, r, R_index):
    if R_index[i - 1] < r <= R_index[i]:
        picked = R_index[i]
        # this index indicates the ith element in index list
        # e.g if this one is 2 one should look at the index list
        # to find the true index for EVENT which is 3
        picked_index = i - 1
        return [picked, picked_index]
###############################################################################
# ------------------------------ Full species ----------------------------------
def full(x):
    if x.species in ["C", "F"] or len(x.bonds) == 0:
        return x.species
    else:
        if x.bonds:
            temp = x.species
            for i in x.bonds:
                if i != x:
                    temp = temp + i.species
            return temp
###############################################################################
# to species the color in visualization
def colors(s):

    if full(s) == "A":
        return 'grey'
    elif full(s) == "E":
        return 'black'
    elif full(s) == "I":
        return 'pink'
    elif full(s) == "EC-Li":
        return 'olive'
    elif full(s) == "P":
        return 'green'
    elif full(s) == "S":
        return 'w'
    elif full(s) == "O":
        return 'orange'
    elif full(s) == "F":
        return 'red'
    elif full(s) == "EC-Li+":
        return 'brown'
    elif full(s) == "OO":
        return 'blue'
    elif full(s) == "C":
        return 'm'
    
###############################################################################
def delt(x, y):
    temp = [x.coordinate[0] - y.coordinate[0], x.coordinate[1] - y.coordinate[1]]
    return temp
###############################################################################
def list_splitter(list_to_split, ratio):
    elements = len(list_to_split)
    middle = int(elements * ratio)
    return [list_to_split[:middle], list_to_split[middle:]]
###############################################################################
# cluster collector function
def cluster_x(single):
    listx = [single]
    stack = [j for j in single.bonds if j.species == "C"]
    while stack != []:
        picked = choice(stack)
        stack = [j for j in stack if j != picked]
        listx.append(picked)
        temp = [k for k in picked.bonds if k.species == "C" and k not in listx]
        stack.extend(temp)
    # to remove repetitive elements
    listx = list(dict.fromkeys(listx))
    return listx
###############################################################################
# electrode size
def EE_func(lattice):
    return [i for i in lattice if i.species == "E"]
###############################################################################
#empty lattice sites
def empty_func(lattice):
    return [i for i in lattice if i.species == "S"]
###############################################################################
# Collecting all the EC-Li+, E in the lattice and return the list
def top_list_func(lattice):
    return [i for i in lattice if i.species in ["E", "EC-Li+"]]
###############################################################################
# calculate the concentration of EC-Li+ in the lattice
def Li_con_func(lattice,xdim,zdim):
    return len([i for i in lattice if i.species == "EC-Li+"])/(xdim*zdim)
##############################################################################