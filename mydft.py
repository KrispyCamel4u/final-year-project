import sys
import math
import time
import csv

N = 100  # number of samples per cycle


def dft(wave):
    x1 = (2/N)*(sum((wave[i]*math.cos(2*math.pi*i/N) for i in range(0, N))))
    x2 = -(2/N)*(sum((wave[i]*math.sin(2*math.pi*i/N) for i in range(0, N))))
    mag = math.sqrt(math.pow(x1, 2)+math.pow(x2, 2))
    phaseAngle = math.atan(x2/x1)
    return mag, math.degrees(phaseAngle)#-(360/N)*(i%N)


# data gathering from csv file
phase = []
with open(sys.argv[1], 'r') as ifile:
    ifile.readline()
    for s in ifile.readlines():
        phase.append(float(s.split(',')[-1]))

# cycle by cycle dft
mCycle = int(len(phase)/N)
output = []
for cycle in range(mCycle):
    params = dft(phase[N*cycle:N+N*cycle])
    output.append([cycle, params[0], params[1]])
    #print(cycle+1, params[0], params[1])

# for i in range(len(phase)-N+1):
#     params=dft(phase[i:i+N],i)
#     output.append([(i+N)*0.0002, params[0], params[1]])
    # print(i+N-1, ' ', (i+N-1)/5000, ' ', dft(phase[i:i+N]))

with open('Output/output.csv', 'w+', newline='') as file: 
    writer = csv.writer(file)
    writer.writerows(output)
'''
# sliding window
for i in range(len(phase)-N+1):
    print(i+N-1, ' ', (i+N-1)/5000, ' ', dft(phase[i:i+N]))
'''
