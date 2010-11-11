def retreive_coordinates(sec):
    ''' Only works with cell which have an xtra mechanism '''
    # Make sure to run h.setpointers() before running
    from neuron import h
    x=sec.sec(.5).x_xtra
    y=sec.sec(.5).y_xtra
    z=sec.sec(.5).z_xtra
    sec.sec.push()
    diam = h.diam3d(0)
    h.pop_section()
    return [x, y, z, diam]

def find_neuron_branches():
    from numpy import array
    from neuron import h
    def cas_index():
        for ii in xrange(len(h.s)):
            if h.s[ii].is_cas():
                return ii
        return -1

    h.s[10].root.push()
    root_index = cas_index()
    h.pop_section()

    branches = []
    def recurse_compartments(index, branches):
        for ii in xrange(int(h.s[index].nchild())):
            h.s[index].child[ii].push()
            child_index = cas_index()
            #print index,',',child_index
            branches.append([index,child_index])
            h.pop_section()
            branches = recurse_compartments(child_index, branches)
        return branches
    return array(recurse_compartments(root_index, branches))

def find_new_point(origPoint, direction, distance):
    ''' determines new point based on the original point
    in the direction of another point with a given length '''
    from numpy import sqrt
    newPoint = [0, 0, 0]
    distance0 = sqrt((direction[0]) ** 2 + (direction[1]) ** 2 + (direction[2]) ** 2)
    distanceRatio = distance / distance0
    for ii in xrange(3):
        newPoint[ii] = origPoint[ii] + (direction[ii]) * distanceRatio
    return newPoint
