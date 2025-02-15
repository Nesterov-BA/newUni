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
    
    
    # def max_difference(array):
    #     return np.diff(array).max()
names, arrays = convert_csv_to_float_arrays("plot1001.csv")
names, arrays = convert_csv_to_float_arrays("plot0.csv")
time = arrays[4]
plt.plot(time,arrays[2], label="0")
names, arrays = convert_csv_to_float_arrays("plot01.csv")
time = arrays[4]
plt.plot(time,arrays[2], label="0.01")
names, arrays = convert_csv_to_float_arrays("plot102.csv")
time = arrays[4]
plt.plot(time,arrays[2], label="1.02")
names, arrays = convert_csv_to_float_arrays("plot1001.csv")
time = arrays[4]
plt.plot(time,arrays[2], label="10.01")
plt.legend()
plt.savefig("plotx1.png")
plt.close()

names, arrays = convert_csv_to_float_arrays("plot0.csv")
time = arrays[4]
plt.plot(time,arrays[1], label="0.0")
names, arrays = convert_csv_to_float_arrays("plot01.csv")
time = arrays[4]
plt.plot(time,arrays[1], label="0.01")
names, arrays = convert_csv_to_float_arrays("plot102.csv")
time = arrays[4]
plt.plot(time,arrays[1], label="1.02")
names, arrays = convert_csv_to_float_arrays("plot1001.csv")
time = arrays[4]
plt.plot(time,arrays[1], label="10.01")
plt.legend()
plt.savefig("plotp2.png")
