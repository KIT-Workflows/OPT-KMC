from random import choice

"""
update function
"""
def update(ev, t, lattice, e1,e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12):
    """"
    ["Surf", "EC-Li+"], ["Surf", "EC-Li"]
    """
    if ev.reactant == e1.reactant and ev.product == e1.product:
        ev.new_coord[0].species = "EC-Li"
        return lattice

    """
    ["EC-Li"], ["P"]
    """
    if ev.reactant == e2.reactant and ev.product == e2.product:
        ev.new_coord[0].species = "P"
        return lattice
    """
    ["Surf", "P"], ["Surf", "I"]
    """
    if ev.reactant == e3.reactant and ev.product == e3.product:
        ev.new_coord[0].species = "I"
        return lattice

    """
    ["I", "EC-Li+"], ["F", "S"]
    """
    if ev.reactant == e4.reactant and ev.product == e4.product:
        ev.new_coord[0].species = "S"
        ev.new_coord[1].species = "F"
        return lattice

    """
    [["P"],["P"]],[["O"],["S"]]
    """
    if ev.reactant == e5.reactant:
        # pick one randomly
        pick = [ev.new_coord[0], ev.new_coord[1]]
        p = choice(pick)
        p.species = "O"
        g = [i for i in pick if i != p]
        g[0].species = "S"
        g[0].bonds = []
        p.bonds = []
        p.status = [[p.coordinate, t], [],'', [0,'']]
        return lattice

    """
    [["O"],["S"]],[["S"],["O"]]
    """
    if ev.reactant == e6.reactant:
        s1 = ev.new_coord[0]
        s2 = ev.new_coord[1]
        s1.species = "S"
        s2.species = "O"
        s2.status = s1.status
        s1.status = [[], [], '', [0,'']]
        s1.bonds = []
        s2.bonds = []
        return lattice

    """
    [["P"],["S"]],[["S"],["P"]]
    """
    if ev.reactant == e7.reactant:
        s1 = ev.new_coord[0]
        s2 = ev.new_coord[1]
        s1.species = "S"
        s2.species = "P"
        s1.bonds = []
        s2.bonds = []
        return lattice

    """
    [["EC-Li"],["S"]],[["S"],["EC-Li"]]
    """
    if ev.reactant == e8.reactant :
        s1 = ev.new_coord[0]
        s2 = ev.new_coord[1]
        s1.species = "S"
        s2.species = "EC-Li"
        s1.bonds = []
        s2.bonds = []
        return lattice

    """
     [["EC-Li+"],["S"]],[["S"],["EC-Li+"]]
    """
    if ev.reactant == e9.reactant or ev.reactant == e10.reactant or ev.reactant == e11.reactant \
            or ev.reactant == e12.reactant:
        s1 = ev.new_coord[0]
        s2 = ev.new_coord[1]
        s1.species = "S"
        s2.species = "EC-Li+"
        s1.bonds = []
        s2.bonds = []
        return lattice
