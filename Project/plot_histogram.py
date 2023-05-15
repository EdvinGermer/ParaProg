import matplotlib.pyplot as plt
import numpy as np

# Read file
with open('histogram.txt', 'r') as f:
    lines = f.read().split('\n')

# Convert each line to a list of integers
boundaries = list(map(int, lines[0].split(', ')))
counts = list(map(int, lines[1].split(', ')))

# Calculate mid-points of each bin for plotting histogram
mid_points = [(boundaries[i]+boundaries[i+1])/2 for i in range(len(boundaries)-1)]

# Plot histogram
plt.bar(mid_points, counts, width=np.diff(boundaries), align='center', edgecolor='black')

# Add labels and title
plt.xlabel('Boundaries')
plt.ylabel('Counts')
plt.title('MC Simulation of Malara Epidemic')

# Show plot
plt.savefig("histogram.png")
