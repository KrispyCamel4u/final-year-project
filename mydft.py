import math
import time

N = 100  # number of samples per cycle


def dft(wave):
    x1 = (2/N)*(sum((wave[i]*math.cos(2*math.pi*i/N) for i in range(N))))
    x2 = (2/N)*(sum((wave[i]*math.sin(2*math.pi*i/N) for i in range(N))))
    mag = math.sqrt(math.pow(x1, 2)+math.pow(x2, 2))
    phaseAngle = math.atan(x2/x1)
    return mag, phaseAngle


# data gathering from csv file
phase = []
with open('C:/Users/camel/Documents/MATLAB/Iabs.csv', 'r') as ifile:
    for s in ifile.readlines():
        phase.append(float(s.split(',')[0]))

# cycle by cycle dft
mCycle = int(len(phase)/100)
for cycle in range(mCycle):
    print(cycle, ' ', dft(phase[N*cycle:N+N*cycle]))

# sliding window
for i in range(len(phase)-N+1):
    print(i+N-1, ' ', (i+N-1)/5000, ' ', dft(phase[i:i+N]))
