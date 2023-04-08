from pickle import load
from matplotlib import pyplot as plot

if __name__ == "__main__":
    with open("Forces.rsiitk", 'rb') as file:
        thrust, gravitational_force, drag_force, net_force, t = load(file)
    
    plot.title("Thrust/Gravitational Force/Air Drag/Mass vs Time")
    plot.xlabel("Time (in seconds)")
    plot.ylabel("Thrust/Gravitational Force/Air Drag/Mass")
    plot.plot(t, thrust, 'r-')
    plot.plot(t, gravitational_force, 'b-')
    plot.plot(t, drag_force, 'g-')
    plot.plot(t, net_force, 'k-')
    plot.legend(["Thrust (in N)", "Gravitational Force (in N)", "Air Drag (in N)", "Net Force (in N)"])
    plot.show()