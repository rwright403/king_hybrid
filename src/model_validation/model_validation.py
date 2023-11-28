import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from src import constants

def read_csv(file_path):   
    np.x = []
    np.y = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:

            if not row:
                continue

            try:
                np.x.append(float(row[0]))
                np.y.append(float(row[1]))
            except ValueError:
                continue

    return np.x, np.y

###ADD FILEPATH HERE!!!!!
model_file_path = constants.model_file_path
data_file_path = constants.data_file_path

xm, ym = read_csv(model_file_path)
xd, yd = read_csv(data_file_path)

plt.plot(xd, yd, label='data')
plt.plot(xm, ym, label='model output')



plt.xlabel('Time (s)')
plt.ylabel('Thrust (N)')
plt.title('Thrust Curve Validation')

for index, fruit in enumerate(yd):
    print(fruit)

plt.legend()
plt.show()