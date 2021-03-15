import sys
import math
import time
import csv

N = 100  # number of samples per cycle
f0=50
del_t=1/(f0*N)

def dft(wave):
    x1={}
    x2={}
    for k in range(3):
        x1[k] = (2/N)*(sum((wave[i-1]*math.cos(2*math.pi*i/N) for i in range(k+1, N+1+k))))
        x2[k] = (2/N)*(sum((wave[i-1]*math.sin(2*math.pi*i/N) for i in range(k+1, N+1+k))))
    # print(x1)

    eq=((x1[2]-x2[1])*math.cos(2*math.pi/N))/((x1[1]-x2[0])*math.cos(4*math.pi/N))
    A= (0.5*N*(x1[0]-x1[1]))/ (math.cos(2*math.pi/N) * eq * (pow(eq,N)-1))
    print(A,eq,(0.5*N*(x1[0]-x1[1])),(math.cos(2*math.pi/N) * eq * (pow(eq,N)-1)))

    for i in range(N):
        wave[i]-=A*pow(eq,i*del_t)

    x11 = (2/N)*(sum((wave[i-1]*math.cos(2*math.pi*i/N) for i in range(0+1, 1+N+k))))
    x21= (2/N)*(sum((wave[i-1]*math.sin(2*math.pi*i/N) for i in range(0+1, 1+N+k))))
    
    mag = math.sqrt(math.pow(x11, 2)+math.pow(x21, 2))
    phaseAngle = math.atan(x11/x21)

    # mag = math.sqrt(math.pow(x1[0], 2)+math.pow(x2[0], 2))
    # phaseAngle = math.atan(x1[0]/x2[0])

    return mag, math.degrees(phaseAngle)


# data gathering from csv file
phase = []
with open(sys.argv[1], 'r') as ifile:
    ifile.readline()
    for s in ifile.readlines():
        phase.append(float(s.split(',')[1]))
# print(phase)
# cycle by cycle dft
mCycle = int(len(phase)/N)
output = []
for cycle in range(mCycle-1):
    params = dft(phase[N*cycle:N+N*cycle+2])
    output.append([cycle, params[0], params[1]])
    #print(cycle+1, params[0], params[1])

# # sliding window
# for i in range(len(phase)-N+1):
#     params=dft(phase[i:i+N])
#     output.append([(i+N)*0.0002, params[0], params[1]])
    # print(i+N-1, ' ', (i+N-1)/5000, ' ', dft(phase[i:i+N]))

with open('DC_offset/Mag_phasor.csv', 'w', newline='') as file: 
    writer = csv.writer(file)
    writer.writerows(output)
'''
# sliding window
for i in range(len(phase)-N+1):
    print(i+N-1, ' ', (i+N-1)/5000, ' ', dft(phase[i:i+N]))
'''
