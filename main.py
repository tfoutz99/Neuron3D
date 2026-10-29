#!/usr/bin/env python
""" Run the Neuron3D examples: python main.py [example number] """
import os
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def example(choice=1):
    print("-"*80)
    filenames = {1:"Examples/1_Spherical_Harmonics.py",
                 2:"Examples/2_Single_Sphere.py",
                 3:"Examples/3_Electrode.py",
                 4:"Examples/4_DBS_Electrode.py",
                 5:"Examples/5_Nuclei.py",
                 6:"Examples/6_Subthalamic_Neuron.py",
                 7:"Examples/7_All_Together.py",
                 8:"Examples/8_Network.py"}
    try:
        print("Executing: %s (close its windows to continue)" % (filenames[choice]))
    except KeyError:
        print("%d is not a valid choice" % (choice))
        return 0
    runpy.run_path(str(ROOT / filenames[choice]), run_name="__main__")
    print("Simulation complete")
    return choice

def interactive(choice=0):
    orig_choice = choice
    print("-"*80)
    print("# Type Selection.  For help type h.")
    print("\ndemo> ", end="", flush=True)
    info = """
Select an example:
    1 - Spherical harmonics
    2 - Single sphere
    3 - Simple electrode
    4 - DBS electrode
    5 - Subthalamic NUCLEUS
    6 - Subthalamic NEURON
    7 - DBS + STN + STh Neuron
    8 - Network model
    b - Run previous
    h - Print this help
    x - Exit interactive Prompt
    q - Quit
"""
    text = sys.stdin.readline().rstrip()
    try:
        choice = int(text)
    except ValueError:
        if text == "b":
            choice -= 2
        elif text == "":
            pass
        elif text == "q":
            print("Quitting...")
            return -2
        elif text == "x":
            return -1
        elif text == "h":
            print(info)
            return orig_choice
        else:
            print("Did not understand input")
            return choice
    if 1 <= choice <= 8:
        example(choice)
        return choice+1
    else:
        if choice == 9:
            print("Simulations complete!")
            return -1
        else:
            print("Invalid selection: %d" % choice)
            return orig_choice

if __name__=='__main__':
    # The examples use paths relative to the project directory
    os.chdir(ROOT)
    choice = 1
    if len(sys.argv) > 1:
        choice = example(int(sys.argv[1])) + 1
    while choice >= 0:
        choice = interactive(choice)
