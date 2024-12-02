import csv
import numpy as np
import matplotlib.pyplot as plt

def extract_data(filename):
    x = []
    y = []

    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            x.append(float(row[0]))
            y.append(float(row[1]))
    return x, y

def trapezoidal_integration(x, y):
    integral = 0.0
    for i in range(1, len(x)):
        dx = x[i] - x[i-1]
        integral += (y[i] + y[i-1]) * dx / 2.0
    return integral

def velocity_integration(filename):
    x, y = extract_data(filename)
    total_velocity = trapezoidal_integration(x, y)
    return total_velocity


def plot_csv(csv_files,labels):
    x_data = []
    y_data = []
    for i in range(len(csv_files)):
        x, u = extract_data(csv_files[i])
        x_data.append(x)
        y_data.append(u)
    for x, y, label in zip(x_data, y_data, labels):
        plt.plot(x, y, label=label)

    plt.xlabel('X [m]')
    plt.ylabel('Velocity [m/s]')
    plt.legend()
    plt.savefig('aaaa.png')


path='/home/stelio/thixotropicSimulation/PreProcessing/CoatingHangerMeng2/'
fileName = 'CoatingHangerMeng2'

Q1_cases = ['SMD4',
            'Thixotropic5a',
            'Thixotropic5b',
            'Thixotropic5c'
            ]
Q1_labels = ['GNM',
             r'$\Lambda = 5x10^{-2}$',
             r'$\Lambda = 1$',
             r'$\Lambda = 5$'
            ]

Q2_cases = ['SMD3',
            'Thixotropic6a',
            'Thixotropic6b',
            'Thixotropic6c'
            ]
Q2_labels = ['GNM',
             r'$\Lambda = 5x10^{-2}$',
             r'$\Lambda = 1$',
             r'$\Lambda = 5$'
             ]


x_data = []
y_data = []

for i in range(len(Q1_cases)):
    x, u = extract_data(f'{path}{fileName}{Q1_cases[i]}.csv')

    x_data.append(x/np.mean([168e-3]))
    y_data.append(u/np.mean(u))
    
plt.figure()
for x, y, label in zip(x_data, y_data, Q1_labels):
        plt.plot(x, y, label=label, linewidth=2)

plt.axhline(y=1, color='black', linestyle='--', linewidth=1, alpha=0.7)

plt.xlabel(r'$x/L$', fontsize=16)
plt.ylabel(r'$u/u_m$', fontsize=16)
plt.legend(loc='upper left', fontsize=14)
plt.ylim(0.6, 1.6)
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)
plt.savefig(f'{path}Q_1.png')
plt.close()

x_data = []
y_data = []

for i in range(len(Q2_cases)):
    x, u = extract_data(f'{path}{fileName}{Q2_cases[i]}.csv')

    x_data.append(x/np.mean([168e-3]))
    y_data.append(u/np.mean(u))
plt.figure()
for x, y, label in zip(x_data, y_data, Q2_labels):
        plt.plot(x, y, label=label, linewidth=2)
plt.axhline(y=1, color='black', linestyle='--', linewidth=1, alpha=0.7)
      
plt.xlabel(r'$x/L$', fontsize=16)
plt.ylabel(r'$u/u_m$', fontsize=16)
plt.legend(loc='upper left', fontsize=14)
plt.ylim(0.6, 1.6)
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)
plt.savefig(f'{path}Q_2.png')
plt.close()