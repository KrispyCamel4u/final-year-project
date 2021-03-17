import sys
import math
import time
import csv

N = 100  # number of samples per cycle
f0=50
del_t=1/(f0*N)

one_time=True
A=0
eq=0

def dft(wave,mcycle):
    x1={}
    x2={}
    global one_time
    global A
    global eq
    for k in range(3):
        x1[k] = (2/N)*(sum((wave[i]*math.cos(2*math.pi*i/N) for i in range(k, N+k))))
        x2[k] = -(2/N)*(sum((wave[i]*math.sin(2*math.pi*i/N) for i in range(k, N+k))))
    # print(x1)
    if one_time:
        eq=((x1[2]-x2[1])*math.cos(2*math.pi/N))/((x1[1]-x2[0])*math.cos(4*math.pi/N))
        A= (0.5*N*(x1[1]-x1[0]))/ (math.cos(2*math.pi/N) * eq * (pow(eq,N)-1))
        print(A,eq,(0.5*N*(x1[1]-x1[0])),(math.cos(2*math.pi/N) * eq * (pow(eq,N)-1)))
        # one_time=False

    for i in range(N):
        print("DDC",A*pow(eq,(i+N*mcycle)*del_t))
        print('org',wave[i])
        wave[i]-=A*pow(eq,(i+N*mcycle)*del_t)
        print('mod',wave[i])

    x11 = (2/N)*(sum((wave[i]*math.cos(2*math.pi*i/N) for i in range(0, N))))
    x21= -(2/N)*(sum((wave[i]*math.sin(2*math.pi*i/N) for i in range(0, N))))
    
    mag = math.sqrt(math.pow(x11, 2)+math.pow(x21, 2))
    phaseAngle = math.atan(x21/x11)

    # mag = math.sqrt(math.pow(x1[0], 2)+math.pow(x2[0], 2))
    # phaseAngle = math.atan(x2[0]/x1[0])

    return mag, math.degrees(phaseAngle)


# data gathering from csv file
phase = []
with open(sys.argv[1], 'r') as ifile:
    ifile.readline()
    for s in ifile.readlines():
        phase.append(float(s.split(',')[-1]))
# print(phase)
# cycle by cycle dft
mCycle = int(len(phase)/N)
output = []
for cycle in range(mCycle-1):
    params = dft(phase[N*cycle:N+N*cycle+2],cycle)
    output.append([cycle, params[0], params[1]])
    #print(cycle+1, params[0], params[1])

# # sliding window
# for i in range(len(phase)-N+1):
#     params=dft(phase[i:i+N])
#     output.append([(i+N)*0.0002, params[0], params[1]])
    # print(i+N-1, ' ', (i+N-1)/5000, ' ', dft(phase[i:i+N]))

with open('Output/output.csv', 'w', newline='') as file: 
    writer = csv.writer(file)
    writer.writerows(output)
'''
# sliding window
for i in range(len(phase)-N+1):
    print(i+N-1, ' ', (i+N-1)/5000, ' ', dft(phase[i:i+N]))
'''
