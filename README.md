# Rocket-Simulator
Simulated the Flight of a Rocket by plotting its Altitude, Vertical Velocity, Vertical Acceleration and Vertical Forces acted upon the Rocket vs Time. <br />
The functions used in the Program for calculating properties like atmospheric density, pressure, temperature, etc are valid upto 10km only, so the program may not work if the rocket goes above 10km

# Requirements
Language Used = Python3<br />
Modules/Packages used:
* math
* sympy
* pickle
* optparse
* datetime
* colorama
* time
* matplotlib

# main.py
It is the python program that does the simulation.
It takes in the following arguments:
* '-m',"--mass" : Mass of the Rocket without the Fuel (Dry Mass)
* '-f',"--fuel" : Mass of theFuel to be loaded in the Rocket
* '-a',"--area" : Area of the Rocket as seen from Top View
* '-o',"--angle" : Vertex Angle of the Right Circular Cone on top of the Rocket (in degrees)
* '-t',"--throttle" : Fuel entering the Engine (in Kg/s) on 100% Throttle (burn rate)
* '-e',"--engine" : Thrust as a function of the fuel entering the engine (x) in Kg/s
* '-r',"--time-resolution" : Time resolution for the simulation (Default = 0.01 seconds)
* '-d',"--displacement-resolution" : Displacement Resolution for the simulation (Default = 0.1 metres)

# read_AVA.py
It reads the pickle dump file that contains data about the Altitude, Vertical Velocity and Vertical Acceleration vs Time created by main.py

# read_focres.py
It reads the pickle dump file that contains data about the Vertical Forces acting upon the Rocket vs Time created by main.py
