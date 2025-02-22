import csv
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

def convert_csv_to_float_arrays(filename):
    try:
        # Open the CSV file
        with open(filename, 'r') as csvfile:
            csv_reader = csv.reader(csvfile)
            
            # Read the header row
            column_names = next(csv_reader)
            
            # Initialize empty lists for each column
            columns = [[] for _ in range(len(column_names))]
            
            # Read rows and populate columns, converting to float
            for row in csv_reader:
                for i, value in enumerate(row):
                    try:
                        columns[i].append(float(value))
                    except ValueError:
                        print(f"Warning: Non-float value '{value}' encountered in column {column_names[i]}. Skipping this value.")
            
            # Convert lists to NumPy arrays with float dtype
            arrays = [np.array(col, dtype=np.float64) for col in columns]
            
            return column_names, arrays
    
    except FileNotFoundError:
        print(f"The file '{filename}' was not found.")
        return None, None
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None, None

def plot_data(arrays, labels, filename, interpolated=False):
    x = arrays[0]
    y = arrays[1]
    z = arrays[2]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if interpolated:
        # Define a grid of points where you want to interpolate the data
        grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]

        # Interpolate the data onto the grid
        grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')

        # Create a surface plot
        ax.plot_surface(grid_x, grid_y, grid_z, cmap='viridis')
    else:
        # Create a  3D scatter plot
        ax.scatter(x, y, z)

    # Set labels for the axes
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])

    # Show the plot
    plt.savefig(filename)

# Usage
column_names, arrays = convert_csv_to_float_arrays('data.csv')


time = arrays[4]
print("Number of points:", len(time))
if column_names and arrays:
    print("Column names:", column_names)
    for i, array in enumerate(arrays):
        if column_names[i] != 'time':
            plt.plot(time,array, label=column_names[i])
plt.legend()
plt.savefig('plot.png')
plt.close()
column_names, arrays = convert_csv_to_float_arrays('dets0.csv')
time = arrays[1]
if column_names and arrays:
    print("Column names:", column_names)
    for i, array in enumerate(arrays):
        if column_names[i] != 'time':
            plt.plot(time,array, label=column_names[i])
plt.legend()
plt.savefig('plotdets0.png')
plt.close()
column_names, arrays = convert_csv_to_float_arrays('dets01.csv')
time = arrays[1]
if column_names and arrays:
    print("Column names:", column_names)
    for i, array in enumerate(arrays):
        if column_names[i] != 'time':
            plt.plot(time,array, label=column_names[i])
plt.legend()
plt.savefig('plotdets01.png')
plt.close()
column_names, arrays = convert_csv_to_float_arrays('dets102.csv')
time = arrays[1]
if column_names and arrays:
    print("Column names:", column_names)
    for i, array in enumerate(arrays):
        if column_names[i] != 'time':
            plt.plot(time,array, label=column_names[i])
plt.legend()
plt.savefig('plotdets102.png')
plt.close()
column_names, arrays = convert_csv_to_float_arrays('dets1001.csv')
time = arrays[1]
if column_names and arrays:
    print("Column names:", column_names)
    for i, array in enumerate(arrays):
        if column_names[i] != 'time':
            plt.plot(time,array, label=column_names[i])
plt.legend()
plt.savefig('plotdets1001.png')
plt.close()


column_names, arrays = convert_csv_to_float_arrays('globalError.csv')
time = arrays[1]
plt.plot(time,arrays[0], label=column_names[0])
plt.legend()
plt.savefig('error.png')
plt.close()
