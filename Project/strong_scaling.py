import matplotlib.pyplot as plt


# Data
x = [1,2,4,8,16,32,64]
y = [140.408184, 70.668614, 35.957631, 19.271297, 9.895813, 5.098819, 2.693381]

y_teo = [y[0]/ x[i] for i in range(len(x))]


# Calculate speedup
speedup = [y[0]/y[i] for i in range(len(y))]
ideal_speedup = x


# Create a subplot
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# First subplot
axs[0].plot(x, y, label="Measurements")
axs[0].plot(x, y_teo, label="Ideal")
axs[0].set_xlabel('Processes', fontsize=14)
axs[0].set_ylabel('Time (s)', fontsize=14)
axs[0].set_title('Strong Scaling - Time Measurements', fontsize=18)
axs[0].legend()
axs[0].grid(True)

# Second subplot
axs[1].plot(x, speedup, label="Speedup")
axs[1].plot(x, ideal_speedup, label="Ideal speedup")
axs[1].set_xlabel('Processes', fontsize=14)
axs[1].set_ylabel('Speedup', fontsize=14)
axs[1].set_title('Strong Scaling - Speedup', fontsize=18)
axs[1].legend()
axs[1].grid(True)

# Settings
plt.rcParams.update({'font.size': 14})


# Show the plot
plt.savefig('strong_scaling.png')
plt.show()
