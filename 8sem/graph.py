import csv
import numpy as np
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

# Usage
column_names, arrays = convert_csv_to_float_arrays('data.csv')

if column_names and arrays:
    print("Column names:", column_names)
    
    # Print first 5 elements of each array
    for i, array in enumerate(arrays):
        print(f"\nArray {i+1} ({column_names[i]}):")
        print(array[:5])  # Print first 5 elements
        
        # Print every 100th element
        print("\nEvery 100th element:")
        hundredth_elements = array[99::100]  # Start from 99th index (100th element) and take every 100th element
        print(hundredth_elements[:10])  # Print first 10 elements (or less if available)
        
        # Print total number of elements
        print(f"\nTotal elements in array: {len(array)}")
        
        # Print data type of the array
        print(f"Data type of array: {array.dtype}")
        
        # Print summary statistics
        print(f"Mean: {np.mean(array):.4f}, Median: {np.median(array):.4f}, Standard Deviation: {np.std(array):.4f}")

# If no arrays were returned, print an error message
else:
    print("No data was loaded successfully.")

# Usage
column_names, arrays = convert_csv_to_float_arrays('data.csv')
time = arrays[4]
if column_names and arrays:
    print("Column names:", column_names)
    for i, array in enumerate(arrays):
        print(len(array))
        print(f"\nArray {i+1} ({column_names[i]}):")
        print(array[:5])  # Print first 5 elements of each array
        plt.plot(time,array)
        plt.show()

