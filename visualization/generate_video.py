import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio

# Replace 'your_data_file.csv' with the path to your data file
data = pd.read_csv('../output_files/positions_data.txt', sep=' ', names=['particle', 'time', 'x', 'y'])

# Prepare the plot
fig, ax = plt.subplots()
scat = ax.scatter(data['x'], data['y'])

# Update function for animation
def update(frame):
    current_data = data[data['time'] == frame]
    scat.set_offsets(current_data[['x', 'y']])
    return scat,

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=data['time'].unique(), interval=10, blit=True)

# Save the animation
ani.save('../output_files/particle_motion.mp4', writer='ffmpeg')

#plt.show()