import csv
import math

# Initialize variables to store max, min, average, and standard deviation of CHL values
max_chl = float("-inf")
min_chl = float("inf")
total_chl = 0
count = 0

# Read the CSV file and find max, min, average, and standard deviation of CHL values
with open('/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/traits/orders/b74733107ab240e64f230e7bb4bfce48/SHIFT_Leaf_Traits_Chl_SB_CA/data/SHIFT_Leaf_Traits_LMA_LWC_Chl.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    chl_values = []
    for row in reader:
        chl_value = float(row['CHL'])
        if chl_value != -9999:  # Exclude invalid CHL values
            max_chl = max(max_chl, chl_value)
            min_chl = min(min_chl, chl_value)
            total_chl += chl_value
            count += 1
            chl_values.append(chl_value)

    if count > 0:
        average_chl = total_chl / count
        standard_deviation_chl = math.sqrt(sum((x - average_chl) ** 2 for x in chl_values) / count)
        print(f'Maximum CHL value: {max_chl} mg.m-2')
        print(f'Minimum CHL value: {min_chl} mg.m-2')
        print(f'Average CHL value: {average_chl} mg.m-2')
        print(f'Standard Deviation of CHL values: {standard_deviation_chl} mg.m-2')
    else:
        print('No valid CHL values found in the CSV file.')

