#!/bin/env ipython -wthread
from enthought.mayavi import mlab
from neuron import h, gui

def example(choice=1):
    print "-"*80
    filenames = {1:"Examples/1_Spherical_Harmonics.py",
                 2:"Examples/2_Single_Sphere.py",
                 3:"Examples/3_Multiple_Spheres.py",
                 4:"Examples/4_Electrode.py",
                 5:"Examples/5_DBS_Electrode.py",
                 6:"Examples/6_Nuclei.py",
                 7:"Examples/7_Subthalamic_Neuron.py",
                 8:"Examples/8_All_Together.py",
                 9:"Examples/9_Network.py"}
    try:
        print """Executing: %s""" % (filenames[choice])
    except KeyError:
        print "%d is not a valid choice" % (choice)
        return 0
    try:
        __IPYTHON__.magic_run("%s" % filenames[choice])
        from enthought.mayavi import mlab
        mlab.show()
    except:
        print "This interactive prompt must be run from ipython with threads -wthread"
    print "Simulation complete"
    return choice

def interactive(choice=0):
    import sys
    orig_choice = choice
    print "-"*80
    print "# Type Selection.  For help type h."
    print "\ndemo>",
    info = """ 
Select an example:
    1 - Spherical harmonics
    2 - Single sphere
    3 - Multiple spheres
    4 - Simple electrode
    5 - DBS electrode
    6 - Subthalamic NUCLEUS
    7 - Subthalamic NEURON
    8 - DBS + STN + STh Neuron
    9 - Network model
    b - Run previous
    h - Print this help
    x - Exit interactive Prompt
    q - Quit
"""
    text = sys.stdin.readline().rstrip()
    try:
        choice = int(text)
    except:
        if text is "b":
            choice -= 2
        elif text is "":
            pass
        elif text is "q":
            print "Quitting..."
            return -2
        elif text is "x":
            return -1
        elif text is "h":
            print info
            return orig_choice
        else:
            print "Did not understand input"
            return choice
    if 1 <= choice <= 9:
        example(choice)
        return choice+1
    else:
        if choice == 10:
            print "Simulations complete!"
            return -1
        else:
            print "Invalid selection: %d" % choice
            return orig_choice

if __name__=='__main__':
    import sys
    try:
        choice = sys.argv[1]
        choice = int(choice)
        example(choice)
        choice += 1
    except:
        choice = interactive(1)
    while choice >= 0:
        choice = interactive(choice)
    if choice == -2:
        __IPYTHON__.magic_Exit()    
