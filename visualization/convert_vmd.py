import pandas as pd

# Replace 'your_data_file.csv' with the path to your data file
data = pd.read_csv('../output_files/positions_data.txt', sep=' ', names=['particle', 'time', 'x', 'y', 'z'])
# Add a Z-coordinate
#data['z'] = 0.0


min_time_step = data['time'].unique()[1] - data['time'].unique()[0]
data['time'] = (data['time']/min_time_step).astype(int)

# Define the atom type (assuming all particles are of the same type)
atom_type = "C"  # You can change this to the appropriate atom type

# Write to XYZ file
with open('../output_files/simulation_data.xyz', 'w') as file:
    for time_step in sorted(data['time'].unique()):
        frame_data = data[data['time'] == time_step]
        file.write(f"{len(frame_data)}\n")
        file.write(f"Time step: {time_step}\n")
        for _, row in frame_data.iterrows():
            file.write(f"{atom_type} {row['x']} {row['y']} {row['z']}\n")

print("XYZ file created: simulation_data.xyz")
