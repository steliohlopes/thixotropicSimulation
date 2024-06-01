import csv
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


path='/home/stelio/thixotropicSimulation/PreProcessing/CoatingHangerMeng/'
plot_csv([f'{path}CoatingHangerMengSMD.csv',
          f'{path}CoatingHangerMengThixotropic1.csv',
          f'{path}CoatingHangerMengThixotropic10.csv',
          'paperMengExp.csv'],['SMD','1','10','Experimental'])
