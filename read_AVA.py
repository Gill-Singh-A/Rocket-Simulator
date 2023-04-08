from pickle import load
from matplotlib import pyplot as plot

if __name__ == "__main__":
    with open("AVA.rsiitk", 'rb') as file:
        position, velocity, acceleration, t = load(file)
    
    time_of_flight = t[-1]

    max_alt = max(position)
    max_alt_time = t[position.index(max_alt)]
    apogee_time = max_alt_time

    max_vel, min_vel = max(velocity), min(velocity)
    max_vel_time = t[velocity.index(max_vel)]
    acceleration_time = max_vel_time

    max_acc, min_acc = max(acceleration), min(acceleration)
    max_acc_time = t[acceleration.index(max_acc)]

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