import matplotlib.pyplot as plt


# Data
x = [1,2,4,8,16,32,64]
y = [140.170206, 140.535529, 141.889092, 146.007486, 160.375434, 161.555665, 180]

y_teo = [y[0] for _ in range(len(x))]


# Calculate efficiency
efficiency = [y[0]/y[i] for i in range(len(y))]
ideal_efficiency = [1 for _ in range(len(x))]


# Create a subplot
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# First subplot
axs[0].plot(x, y, label="Measurements")
axs[0].plot(x, y_teo, label="Ideal")
axs[0].set_xlabel('Processes')
axs[0].set_ylabel('Time (s)')
axs[0].set_title('Weak Scaling of MC Simulation')
axs[0].legend()
axs[0].grid(True)

# Second subplot
axs[1].plot(x, efficiency, label="Efficiency")
axs[1].plot(x, ideal_efficiency, label="Ideal Efficiency")
axs[1].set_xlabel('Processes')
axs[1].set_ylabel('Efficiency')
axs[1].set_title('Efficiency of MC Simulation')
axs[1].legend()
axs[1].grid(True)

# Settings
plt.rcParams.update({'font.size': 14})

# Show the plot
plt.savefig('weak_scaling.png')
plt.show()
