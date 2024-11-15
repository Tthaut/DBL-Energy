import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def read_csv(path):
    input_file = path  # Original CSV file path
    output_file = path.replace('.csv','_fixed.csv')  # Fixed CSV file path

    # Process the file to fix delimiters and remove commas in numbers
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for i, line in enumerate(infile):
            # For the first line, replace the first comma with a semicolon
            if i == 0:
                line = line.replace(",", ";", 1)
                line = line.replace(' ', '',1)
            else:
                # For other lines, replace the second comma with a semicolon
                parts = line.split(",")
                if len(parts) > 2:
                    line = ",".join(parts[:2]) + ";" + ",".join(parts[2:])
                    line = line.replace(' ', '',1)
                    line = line.replace(',','.',2)
            outfile.write(line)

    # Read the fixed CSV file using numpy
    data = np.genfromtxt(output_file, delimiter=";", dtype=None, names=True, encoding='ISO-8859-1')
    
    return data

data_0 = read_csv("Data\\13-11-24\\azobenzene 50uM.Sample.Raw.csv")
data_5min = read_csv('Data\\13-11-24\\aB 50uM 5min 365nm.Sample.Raw.csv')
data_10min = read_csv('Data\\13-11-24\\aB 50uM 10min 365nm.Sample.Raw.csv')
data_15min = read_csv('Data\\13-11-24\\aB 50uM 15min 365nm.Sample.Raw.csv')
data_45min = read_csv('Data\\13-11-24\\aB 50uM 45min 365nm.Sample.Raw.csv')
data_60min = read_csv('Data\\13-11-24\\aB 50um 60min 365nm.Sample.Raw.csv')

plt.plot(data_0['nm'], data_0['A'], label='0 min', color='blue')  # Plot data_0
plt.plot(data_5min['nm'], data_5min['A'], label='5 min', color='green')  # Plot data_5min
plt.plot(data_10min['nm'], data_10min['A'], label='10 min', color='red')  # Plot data_10min
plt.plot(data_15min['nm'], data_15min['A'], label='15 min', color='orange')  # Plot data_15min
plt.plot(data_45min['nm'], data_45min['A'], label='45 min', color='purple')  # Plot data_45min
plt.plot(data_60min['nm'], data_60min['A'], label='60 min', color='brown')  # Plot data_60min


max_y_value = np.max(data_0['A'])
max_y_index = np.argmax(data_0['A'])
max_x_value = data_0['nm'][max_y_index]  # The corresponding x-value for max y

plt.axvline(max_x_value, color='black', linestyle='--', linewidth=1, label=f'Max at nm = {max_x_value:.2f}')  # Vertical line at max y
plt.axvline(365, color='black', linestyle='-.', linewidth=0.5, alpha=0.5, label='Compounds exposed with 365nm light')
# Add labels and title
plt.axhline(0, color='black', linestyle='--', linewidth=1)  # Horizontal line at y=0
plt.xlim(275, 800)  # Set x-axis limits
plt.ylim(-0.5, 1.1)  # Set y-axis limits
plt.xlabel('Wavelength (nm)')  # Label x-axis
plt.ylabel('Absorbance (A)')  # Label y-axis
plt.suptitle('Absorbance vs. Wavelength for Azobenzene Samples', weight='bold')  # Title of the plot
plt.title('Samples Exposed to 365nm Light')
# Add legend to identify each line
plt.legend()

# Show the plot
plt.show()

## Plotting time vs abs_max

time = [0,5,10,15,45,60]
# Collect maximum absorbance values at each time point
abs_max = [np.max(data_0["A"]), np.max(data_5min["A"]), np.max(data_10min["A"]),np.max(data_15min["A"]), np.max(data_45min["A"]), np.max(data_60min["A"])]

def exponential_decay(t, A, B, C):
    return A * np.exp(-B * t) + C

# Fit the exponential decay curve to the time vs abs_max data
params, covariance = curve_fit(exponential_decay, time, abs_max, p0=[1, 0.1, 0])

# Extract the fitted parameters
A_fit, B_fit, C_fit = params

# Generate fitted data points
time_fine = np.linspace(min(time), max(time), 1000)  # Fine time points for a smooth curve
abs_max_fit = exponential_decay(time_fine, A_fit, B_fit, C_fit)

# Plot Time vs Maximum Absorbance
plt.figure(figsize=(8, 6))
plt.plot(time, abs_max, marker='o', linestyle='-', color='b', label='Max Absorbance')
plt.plot(time_fine, abs_max_fit, linestyle='--', color='r', label='Exponential Fit')

# Labels and Title
plt.xlabel('Time (min)')
plt.ylabel('Maximum Absorbance (A)')
plt.title('Maximum Absorbance with Exponential Decay Fit')
plt.grid(True)
plt.legend()

# Show the plot
plt.show()

# Print fitted parameters
print(f"Fitted parameters: A = {A_fit:.4f}, B = {B_fit:.4f}, C = {C_fit:.4f}")