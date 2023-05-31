import math, sympy
from pickle import dump
from optparse import OptionParser
from datetime import date
from colorama import Fore, Back, Style
from time import strftime, localtime, time
from matplotlib import pyplot as plot

G = 6.6743e-11
M = 5.97219e24
R = 6.3781e6

status_color = {
    '+': Fore.GREEN,
    '-': Fore.RED,
    '*': Fore.YELLOW,
    ':': Fore.CYAN,
    ' ': Fore.WHITE
}

def display(status, data):
    print(f"{status_color[status]}[{status}] {Fore.BLUE}[{date.today()} {strftime('%H:%M:%S', localtime())}] {status_color[status]}{Style.BRIGHT}{data}{Fore.RESET}{Style.RESET_ALL}")
def display_dictionary(heading, dict):
    display('+', heading)
    spacing = max([len(key) for key in dict.keys()])+5
    for key, value in dict.items():
        display(':', f"\t{key}{' '*(spacing-len(key))}{Back.MAGENTA}{value}{Back.RESET}")

def get_arguments(*args):
    parser = OptionParser()
    for arg in args:
        parser.add_option(arg[0], arg[1], dest=arg[2], help=arg[3])
    return parser.parse_args()[0]

class Planet():
    def __init__(self, mass, radius):
        self.mass = mass
        self.radius = radius
        self.volume = 4*math.pi*math.pow(self.radius, 3)/3
        self.density = self.mass/self.volume
        self.surface_temperature = 288.16
        self.temperature_gradient = 0.0065
        self.surface_pressure = 101325
        self.surface_atmospheric_denstiy = 1.225
        self.adiabatic_exponent = 1.2349
        self.power_factor = self.adiabatic_exponent/(self.adiabatic_exponent-1)
    def feild(self, altitude):
        return G*self.mass/math.pow(altitude+self.radius, 2)
    def atmospheric_density(self, altitude=0):
        return self.surface_atmospheric_denstiy*math.pow(self.temperature(altitude)/self.surface_temperature, self.power_factor-1)
    def pressure(self, altitude=0):
        return self.surface_pressure*math.pow(self.temperature(altitude)/self.surface_temperature, self.power_factor)
    def temperature(self, altitude=0):
        return self.surface_temperature-self.temperature_gradient*altitude
    def speed_of_sound(self, altitude=0):
        return math.pow(self.adiabatic_exponent*self.pressure(altitude)/self.atmospheric_density(altitude), 0.5)
class Rocket():
    def __init__(self, mass, fuel_mass, area, vertex_angle_of_cone, throttle, thrust_function):
        self.mass = mass
        self.fuel_mass = fuel_mass
        self.area = area
        self.throttle = throttle
        self.vertex_angle_of_cone = vertex_angle_of_cone*math.pi/180
        self.thrust_function = sympy.sympify(thrust_function)
        self.thrust = self.thrust_function.subs('x', self.throttle).evalf()
        self.position = 0
        self.velocity = 0
        self.acceleration = 0
        self.drag_coefficient = 1.17
    def total_mass(self):
        return self.mass+self.fuel_mass
    def mach_number(self, Planet):
        return self.velocity/Planet.speed_of_sound(self.position)
    #def drag_coefficient(self, planet):
    #    return 0.5 * (1 + math.pow(math.tan(self.vertex_angle_of_cone), 2)) * (2/(math.pow(self.mach_number(planet), 2) * math.cos(self.vertex_angle_of_cone))) * (1 - (1/((1 + 0.5 * math.pow(self.mach_number(planet), 2) * math.pow(math.sin(self.vertex_angle_of_cone), 2))^(1.5))))
    def air_drag(self, planet):
        return 0.5*planet.atmospheric_density(self.position)*math.pow(self.velocity, 2)*self.drag_coefficient*self.area

def simulate(data):
    display('+', "Rocket Properties")
    max_spacing = max([len(key) for key in vars(data)]) + 5
    display(':', f"Mass{' '*(max_spacing-4)}{Back.MAGENTA}{data.mass}{Back.RESET} Kilogram")
    display(':', f"Fuel Mass{' '*(max_spacing-9)}{Back.MAGENTA}{data.fuel_mass}{Back.RESET} Kilogram")
    display(':', f"Area{' '*(max_spacing-4)}{Back.MAGENTA}{data.area}{Back.RESET} Square Metre")
    display(':', f"Angle{' '*(max_spacing-5)}{Back.MAGENTA}{data.angle}{Back.RESET} Degrees")
    display(':', f"Throttle{' '*(max_spacing-8)}{Back.MAGENTA}{data.throttle}{Back.RESET} Kilogram / second")
    display(':', f"Thrust Function{' '*(max_spacing-15)}{Back.MAGENTA}{data.thrust_function}{Back.RESET}")

    rocket = Rocket(float(data.mass), float(data.fuel_mass), float(data.area), float(data.angle), float(data.throttle), str(data.thrust_function))

    print()

    display('+', "Planet Properties")
    earth = Planet(M, R)
    display(':', f"Mass{' '*28}{Back.MAGENTA}{earth.mass}{Back.RESET} Kilograms")
    display(':', f"Radius{' '*26}{Back.MAGENTA}{earth.radius}{Back.RESET} Metres")
    display(':', f"Volume{' '*26}{Back.MAGENTA}{earth.volume}{Back.RESET} Cubic Metres")
    display(':', f"Density{' '*25}{Back.MAGENTA}{earth.density}{Back.RESET} Kilograms / Cubic Metres")
    display(':', f"Surface Temperature{' '*13}{Back.MAGENTA}{earth.surface_temperature}{Back.RESET} Kelvin")
    display(':', f"Temperature Gradient{' '*12}{Back.MAGENTA}{earth.temperature_gradient}{Back.RESET} Kelvin / Metre")
    display(':', f"Surface Pressure{' '*16}{Back.MAGENTA}{earth.surface_pressure}{Back.RESET} Pascal")
    display(':', f"Surface Atomospheric Density{' '*4}{Back.MAGENTA}{earth.surface_atmospheric_denstiy}{Back.RESET} Kilograms / Cubic Metres")
    display(':', f"Adiabatic Exponent{' '*14}{Back.MAGENTA}{earth.adiabatic_exponent}{Back.RESET}")
    display(':', f"Power Factor{' '*20}{Back.MAGENTA}{earth.power_factor}{Back.RESET}")

    print()

    display('+', "Simulation Settings")
    display(':', f"Time Resolution{' '*13}{Back.MAGENTA}{data.time_resolution}{Back.RESET} Seconds")
    display(':', f"Displacement Resolution{' '*5}{Back.MAGENTA}{data.displacement_resolution}{Back.RESET} Metres")

    print()

    display('+', f"Graphing Earth's Properties")
    t1 = time()
    feild, atmospheric_density, pressure, temperature, speed_of_sound, distances, distance = [], [], [], [], [], [], 0
    while distance < 10000:
        feild.append(earth.feild(distance))
        atmospheric_density.append(earth.atmospheric_density(distance))
        pressure.append(earth.pressure(distance))
        temperature.append(earth.temperature(distance))
        speed_of_sound.append(earth.speed_of_sound(distance))
        distances.append(distance)
        distance += data.displacement_resolution
    t2 = time()
    display('+', f"Done")
    display(':', f"Time Taken = {Back.MAGENTA}{t2-t1}{Back.RESET} seconds")

    print()
    
    plot.title("Gravitational Feild vs Altitude")
    plot.xlabel("Altitude (in metres)")
    plot.ylabel("Gravitational Feild (in N/Kg)")
    plot.plot(distances, feild)
    plot.show()

    plot.title("Atmospheric Density vs Altitude")
    plot.xlabel("Altitude (in metres)")
    plot.ylabel("Atmospheric Density (in Kilograms / Cubic Metres)")
    plot.plot(distances, atmospheric_density)
    plot.show()

    plot.title("Pressure vs Altitude")
    plot.xlabel("Altitude (in metres)")
    plot.ylabel("Pressure (in Pascal)")
    plot.plot(distances, pressure)
    plot.show()

    plot.title("Temperature")
    plot.xlabel("Altitude (in metres)")
    plot.ylabel("Temperature (in Kelvin)")
    plot.plot(distances, temperature)
    plot.show()

    plot.title("Speed of Sound")
    plot.xlabel("Altitude (in metres)")
    plot.ylabel("Speed of Sound (in metres / second)")
    plot.plot(distances, speed_of_sound)
    plot.show()
    
    print()

    t1 = time()
    display('+', f"Starting Simulation")
    position, velocity, acceleration, thrust, gravitational_force, drag_force, net_force, mass, t, total_time = [], [], [], [] ,[] ,[], [], [], [], 0
    while rocket.fuel_mass > 0:
        t.append(total_time)
        position.append(rocket.position)
        velocity.append(rocket.velocity)
        thrust.append(rocket.thrust)
        mass.append(rocket.total_mass())
        gravitational_force.append(mass[-1]*earth.feild(rocket.position))
        drag_force.append(rocket.air_drag(earth))
        net_force.append(thrust[-1]-(gravitational_force[-1]+drag_force[-1]))
        rocket.acceleration = net_force[-1]/mass[-1]
        rocket.velocity += rocket.acceleration * data.time_resolution
        rocket.position += rocket.velocity * data.time_resolution
        if rocket.position < 0:
            rocket.position = 0
            if rocket.velocity < 0:
                rocket.velocity = 0
            if rocket.acceleration < 0:
                rocket.acceleration = 0
        acceleration.append(rocket.acceleration)
        rocket.fuel_mass -= rocket.throttle * data.time_resolution
        total_time += data.time_resolution
    t2 = time()

    rocket.fuel_mass = 0

    acceleration_time = total_time
    simulation_acceleration_time = t2-t1
    display('*', f"Simulation for Acceleration Done")
    display(':', f"Time taken to Simulate = {Back.MAGENTA}{simulation_acceleration_time}{Back.RESET} seconds")
    display(':', f"Time while acceleration = {Back.MAGENTA}{acceleration_time}{Back.RESET} seconds")

    t1 = time()
    while rocket.velocity > 0:
        t.append(total_time)
        thrust.append(0)
        mass.append(rocket.total_mass())
        gravitational_force.append(mass[-1]*earth.feild(rocket.position))
        drag_force.append(rocket.air_drag(earth))
        net_force.append(-(gravitational_force[-1]+drag_force[-1]))
        rocket.acceleration = net_force[-1]/mass[-1]
        rocket.velocity += rocket.acceleration * data.time_resolution
        rocket.position += rocket.velocity * data.time_resolution
        if rocket.position < 0:
            rocket.position = 0
            if rocket.velocity < 0:
                rocket.velocity = 0
            if rocket.acceleration < 0:
                rocket.acceleration = 0
        acceleration.append(rocket.acceleration)
        velocity.append(rocket.velocity)
        position.append(rocket.position)
        total_time += data.time_resolution
    t2 = time()

    apogee_time = total_time
    simulation_apogee_time = t2-t1
    display('*', f"Simulation till Apogee Done")
    display(':', f"Time taken to Simulate = {Back.MAGENTA}{simulation_apogee_time}{Back.RESET} seconds")
    display(':', f"Time taken to reach Apogee = {Back.MAGENTA}{apogee_time}{Back.RESET} seconds")

    t1 = time()
    rocket.drag_coefficient = 0.82
    while rocket.position > 0:
        t.append(total_time)
        thrust.append(0)
        mass.append(rocket.total_mass())
        gravitational_force.append(mass[-1]*earth.feild(rocket.position))
        drag_force.append(rocket.air_drag(earth))
        net_force.append(drag_force[-1]-gravitational_force[-1])
        rocket.acceleration = net_force[-1]/mass[-1]
        rocket.velocity += rocket.acceleration * data.time_resolution
        rocket.position += rocket.velocity * data.time_resolution
        acceleration.append(rocket.acceleration)
        velocity.append(rocket.velocity)
        position.append(rocket.position)
        total_time += data.time_resolution
    t2 = time()

    time_of_flight = total_time
    simulation_declining_time = t2-t1
    display('*', f"Complete Simulation Done")
    display(':', f"Time taken to Simulate = {Back.MAGENTA}{simulation_declining_time}{Back.RESET} seconds")
    display(':', f"Ground Hit Velocity = {Back.MAGENTA}{rocket.velocity}{Back.RESET} metres / second")
    display(':', f"Time of Flight = {Back.MAGENTA}{time_of_flight}{Back.RESET} seconds")

    print()

    display('+', f"Total Simulation Time = {Back.MAGENTA}{simulation_acceleration_time+simulation_apogee_time+simulation_declining_time}{Back.RESET} seconds")

    print()

    max_alt = max(position)
    max_alt_time = t[position.index(max_alt)]
    display(':', f"Maximum Altitude = {Back.MAGENTA}{max_alt}{Back.RESET} metres @ {Back.MAGENTA}{max_alt_time}{Back.RESET} seconds")
    plot.title("Altitude vs Time")
    plot.xlabel("Time (in seconds)")
    plot.ylabel("Altitude (in metres)")
    plot.plot(t, position)
    plot.plot((0, 0), (0, max_alt), 'k-')
    plot.plot((acceleration_time, acceleration_time), (0, max_alt), 'y-')
    plot.plot((apogee_time, apogee_time), (0, max_alt), 'c-')
    plot.plot((time_of_flight, time_of_flight), (0, max_alt), 'm-')
    plot.plot(max_alt_time, max_alt, 'r-o')
    plot.legend(["Altitude (m)", "Launch", "Burnout", "Apogee", "Ground Hit", "Maximum Altitude (m)"])
    plot.show()

    max_vel, min_vel = max(velocity), min(velocity)
    max_vel_time = t[velocity.index(max_vel)]
    display(':', f"Maximum Velocity = {Back.MAGENTA}{max_vel}{Back.RESET} metres / second @ {Back.MAGENTA}{max_vel_time}{Back.RESET} seconds")
    plot.title("Velocity vs Time")
    plot.xlabel("Time (in seconds)")
    plot.ylabel("Velocity (in metres / second)")
    plot.plot(t, velocity)
    plot.plot((0, 0), (min_vel, max_vel), 'k-')
    plot.plot((acceleration_time, acceleration_time), (min_vel, max_vel), 'y-')
    plot.plot((apogee_time, apogee_time), (min_vel, max_vel), 'c-')
    plot.plot((time_of_flight, time_of_flight), (min_vel, max_vel), 'm-')
    plot.plot(max_vel_time, max_vel, 'r-o')
    plot.legend(["Velocity (m/s)", "Launch", "Burnout", "Apogee", "Ground Hit", "Maximum Velocity (m/s)"])
    plot.show()

    max_acc, min_acc = max(acceleration), min(acceleration)
    max_acc_time = t[acceleration.index(max_acc)]
    display(':', f"Maximum Acceleration = {Back.MAGENTA}{max_acc}{Back.RESET} metres / second square @ {Back.MAGENTA}{max_acc_time}{Back.RESET} seconds")
    plot.title("Acceleration vs Time")
    plot.xlabel("Time (in seconds)")
    plot.ylabel("Acceleration (in metres / second sqaure)")
    plot.plot(t, acceleration)
    plot.plot((0, 0), (min_acc, max_acc), 'k-')
    plot.plot((acceleration_time, acceleration_time), (min_acc, max_acc), 'y-')
    plot.plot((apogee_time, apogee_time), (min_acc, max_acc), 'c-')
    plot.plot((time_of_flight, time_of_flight), (min_acc, max_acc), 'm-')
    plot.plot(max_acc_time, max_acc, 'r-o')
    plot.legend(["Acceleration (m/s^2)", "Launch", "Burnout", "Apogee", "Ground Hit", "Maximum Acceleration (m/s^2)"])
    plot.show()

    max_val, min_val = max(max_alt, max_vel, max_acc), min(min_vel, min_acc)
    plot.title("Altitude/Velocity/Acceleration vs Time")
    plot.xlabel("Time (in seconds)")
    plot.ylabel("Altitude/Velocity/Acceleration")
    plot.plot(t, position, 'r-')
    plot.plot(t, velocity, 'b-')
    plot.plot(t, acceleration, 'g-')
    plot.plot((0, 0), (min_val, max_val), 'k-')
    plot.plot((acceleration_time, acceleration_time), (min_val, max_val), 'y-')
    plot.plot((apogee_time, apogee_time), (min_val, max_val), 'c-')
    plot.plot((time_of_flight, time_of_flight), (min_val, max_val), 'm-')
    plot.plot(max_alt_time, max_alt, 'r-o')
    plot.plot(max_vel_time, max_vel, 'b-o')
    plot.plot(max_acc_time, max_acc, 'g-o')
    plot.legend(["Altitude (m)", "Velocity (m/s)", "Acceleration (m/s^2)", "Launch", "Burnout", "Apogee", "Ground Hit",  "Maximum Altitude (m)", "Maximum Velocity (m/s)", "Maximum Acceleration (m/s^2)"])
    plot.show()

    #max_gravitational_force, max_drag_force, max_net_force, max_mass = max(gravitational_force), max(drag_force), max(net_force), mass[0]
    #min_gravitational_force, min_drag_force, min_net_force, min_mass = min(gravitational_force), min(drag_force), min(net_force), mass[0]
    #max_val, min_val = max(thrust[0], max_gravitational_force, max_drag_force, max_net_force, max_mass), min(min_gravitational_force, min_drag_force, min_drag_force, min_mass, min_net_force)
    plot.title("Thrust/Gravitational Force/Air Drag/Mass vs Time")
    plot.xlabel("Time (in seconds)")
    plot.ylabel("Thrust/Gravitational Force/Air Drag/Mass")
    plot.plot(t, thrust, 'r-')
    plot.plot(t, gravitational_force, 'b-')
    plot.plot(t, drag_force, 'g-')
    plot.plot(t, net_force, 'k-')
    plot.legend(["Thrust (in N)", "Gravitational Force (in N)", "Air Drag (in N)", "Net Force (in N)"])
    plot.show()

    with open("AVA.rsiitk", 'wb') as file:
        dump((position, velocity, acceleration, t), file)
    with open("Forces.rsiitk", 'wb') as file:
        dump((thrust, gravitational_force, drag_force, net_force, t), file)

if __name__ == "__main__":
    data = get_arguments(('-m', '--mass', "mass", "Mass of the Rocket without the Fuel (Dry Mass)"),
                         ('-f', '--fuel', "fuel_mass", "Mass of the Fuel to be loaded in the Rocket"),
                         ('-a', '--area', "area", "Area of the Rocket as Seen from Top View"),
                         ('-o', '--angle', "angle", "Vertex Angle of the Right Circular Cone on top of the Rocket (in degrees)"),
                         ('-t', '--throttle', "throttle", "Fuel entering the Engine (in Kg/s) on 100% throttle"),
                         ('-e', '--engine', "thrust_function", "Thrust as a function of the fuel entering the engine (x) in Kg/s"),
                         ('-r', '--time-resolution', "time_resolution", "Time resolution for the simulation (Default = 0.01 seconds)"),
                         ('-d', '--displacement-resolution', "displacement_resolution", "Displacement Resolution for the simulation (Default = 0.1 metres)"))
    if not data.time_resolution:
        data.time_resolution = 0.01
    if not data.displacement_resolution:
        data.displacement_resolution = 0.1
    simulate(data)