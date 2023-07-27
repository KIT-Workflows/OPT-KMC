from .support import *

"""                             Collection function                         """
################################################################################
def pre_event(top_list, Li_con, events, e1,e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12):
    pre_list = list()
    del pre_list[:]
    for s in top_list:
        # determine the type of sites
        ss = full(s)
        NBR = [i for i in s.nbr if i not in s.bonds]
        for nbr in NBR:
            ns = full(nbr)
            #-------------------------------  EC-Li+  ----------------------------------
            if ss == "EC-Li+":
                #reduction of EC-Li+ to EC-Li
                if ns in ["E","F"]:
                # checking for electron reduction close  SEI inorganic layer, and electrode
                        if nbr.coordinate[1] <= 4 and s.coordinate[1] <= 5:
                            pre_list.append(New_event(e1.reactant, e1.product, e1.barrier, e1.rate, [s, nbr, events.index(e1)]))
                # EC-Li+ diffusion
                elif ns == "S":
                        # which direction
                        deltaX = delt(nbr, s)
                        if deltaX == [1, 0]:
                            '''RIGHT'''
                            ev = e9
                            # modifying the rate with FACTOR
                            pre_list.append(
                                New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))
                        '''LEFT'''
                        if deltaX == [-1, 0]:
                            ev = e10
                            # modifying the rate with FACTOR
                            pre_list.append(
                                New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))
                        '''UP'''
                        if deltaX == [0, 1]:
                            ev = e11
                            # modifying the rate with FACTOR
                            pre_list.append(
                                New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))
                        '''DOWN'''
                        if deltaX == [0, -1]:
                            ev = e12
                            pre_list.append(New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))

                elif ns == "I":
                        ev = e4
                        pre_list.append(
                            New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))
            #-------------------------------  EC-Li  ----------------------------------
            elif ss == "EC-Li":
                # P
                ev = e2
                pre_list.append(New_event(ev.reactant, ev.product, ev.barrier, ev.rate,
                                              [s, s, events.index(ev)]))
                #diffusion
                if ns == "S":
                        ev = e8
                        pre_list.append(New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))

            #-------------------------------  All P reacions ----------------------------------
            elif ss == "P" :
                # Li2EDC production
                if ns == "P":
                    ev = e5
                    pre_list.append(New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))
                #Diffusion
                elif ns == "S":
                    ev = e7
                    pre_list.append(New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))

                elif  ns in ["E", "F"] and nbr.coordinate[1] <= 4 and s.coordinate[1] <= 4 :
                    # I
                    ev = e3
                    pre_list.append(New_event(ev.reactant, ev.product, ev.barrier,ev.rate,
                                                [s, nbr, events.index(ev)]))

            #-------------------------------  All O reacions ----------------------------------
            elif ss == "O":
                # Li2EDC_2 production
                if ss == "S":
                        ev = e6
                        pre_list.append(
                            New_event(ev.reactant, ev.product, ev.barrier, ev.rate, [s, nbr, events.index(ev)]))
    return pre_list
