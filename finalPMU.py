import math
import cmath
import sys
import csv
import os

import matplotlib.pyplot as plt
import filter
from SVA import shift_samples
from TVE import TVE

root = os.path.join(os.path.abspath(os.path.dirname(__file__)), './')

N = 100  # number of samples per cycle
f0=50
# phase = 60
dt = 1/(50*N)
cutoff=50
fs=1/dt
order=8
offset=5


Bz,Az=filter.butter(N=order,Wn=cutoff,btype='low',fs=fs)

def parabola_aprox(wave,cycle,sign=-1):
    wave=list(wave)
    wave_sort=sorted(wave)
    y1=wave_sort[sign*1]
    inm1=wave.index(y1)
    of=1
    if sign==-1:
        of=-1
    if inm1==0:
        y1=wave_sort[sign*1+of]
        inm1=wave.index(y1)
        # inm2=1
        # inm3=2
    elif inm1==N-1+offset:
        y1=wave_sort[sign*1+of]
        inm1=wave.index(y1)
        # inm2=inm1-1
        # inm3=inm2-1
    # else:
    inm2=inm1-1
    inm3=inm1+1

    y2=wave[inm2]
    y3=wave[inm3]

    x1=cycle*N*dt+inm1*dt
    x2=cycle*N*dt+inm2*dt
    x3=cycle*N*dt+inm3*dt
    k=( math.pow(x2,2)*(y1-y3) + math.pow(x1,2)*(y3-y2) + math.pow(x3,2)*(y2-y1))/( x2*(y1-y3) + x1*(y3-y2) + x3*(y2-y1) )
    return k/2

def pricise_freq(wave,cycle):
    
    tm1=parabola_aprox(wave,cycle,-1)
    tm2=parabola_aprox(wave,cycle,0)
    T=2*abs(tm1-tm2)
    return 1/T

# filtered = butter_lowpass_filter(noise, cutoff, fs, order)

def jexp(angle):
    return complex(math.cos(angle), math.sin(angle))

def dft(wave):
    x1 = (2/N)*(sum((wave[i]*math.cos(2*math.pi*i/N) for i in range(0, N))))
    x2 = -(2/N)*(sum((wave[i]*math.sin(2*math.pi*i/N) for i in range(0, N))))
    #mag = math.sqrt(math.pow(x1, 2)+math.pow(x2, 2))
    phaseAngle = math.atan(x2/x1)
    return math.degrees(phaseAngle)#-(360/N)*(i%N)

def dft_modified(wave):

    I_dft = (2/N)*(sum([wave[k]*jexp(-2*math.pi*k/N) for k in range(N)]))

    I_even = (2/N)*(sum([wave[2*k]*jexp(-2*math.pi*2*k/N)
                         for k in range(int(N/2))]))

    I_odd = (2/N)*(sum([wave[2*k+1]*jexp(-2*math.pi*(2*k+1)/N)
                        for k in range(int(N/2))]))

    Kre = (I_even-I_odd).real
    Kimg = (I_even-I_odd).imag

    # E=e^(-del_t/tau)
    E = (Kimg)/(Kre*math.sin(2*math.pi/N)-Kimg*math.cos(2*math.pi/N))

    I_dc_dft = (I_even-I_odd)*(1+E*jexp(-2*math.pi/N))/(1-E*jexp(-2*math.pi/N))

    I_fundamental_dft = I_dft-I_dc_dft

    return abs(I_fundamental_dft), math.degrees(cmath.phase(I_fundamental_dft))

INsyncdata=[0]*order
OUTsyncdata=[0]*order
mx=0
prev_frequency=f0
def my_pmu(wave,cycle):
    # freq estimation
    global mx,INsyncdata,OUTsyncdata,shifter,prev_frequency
    filtered_wave=filter.lfilter(Bz,Az, wave,True,cycle,INsyncdata,OUTsyncdata)
    est_freq=pricise_freq(filtered_wave,cycle)
    INsyncdata=wave[-order-offset:-offset]
    OUTsyncdata=filtered_wave[-order-offset:-offset]
    if cycle>10 and cycle<mCycle-5:
        if mx<abs(est_freq-f0):
            mx=abs(est_freq-f0)
    # /freq estimation
    # ROCOF
    est_rocof=(est_freq-prev_frequency)*f0
    # rocof end
    prev_frequency = est_freq
    if cycle>=5:
        if abs(est_freq-f0)>0.01:

            shifted_samples = shifter.shift(wave[:-offset+1], est_freq, N)
            est_mag, est_phase = dft_modified(shifted_samples)

        else:
            est_mag, est_phase= dft_modified(wave[:-offset])

        return round(est_mag, 4), round(est_phase, 4), round(est_freq, 4), round(est_rocof, 4)
    return 0,0,0,0

# data gathering from csv file
input_data = []
with open("Dataset/normal.csv", 'r') as ifile:
    ifile.readline()
    for s in ifile.readlines():
        input_data.append(float(s.split(',')[1]))

with open(f"{root}/{sys.argv[1]}", 'r') as ifile:
    ifile.readline()
    for s in ifile.readlines():
        input_data.append(float(s.split(',')[1]))
        input_data2.append(float(s.split(',')[2]))
        input_data3.append(float(s.split(',')[3]))
        input_data4.append(float(s.split(',')[4]))

off_angle_51 = []
off_angle_52 = []
off_angle_50_25 = []
with open("Dataset/off_freq_angle.csv", 'r') as ifile:
    for s in ifile.readlines():
        off_angle_50_25.append(float(s.split(',')[9]))
        off_angle_52.append(float(s.split(',')[16]))

# cycle by cycle dft
mCycle = int(len(input_data)/N)
output = []

shifter = shift_samples() #shift samples
tve = TVE()
for cycle in range(mCycle-1):
    params=my_pmu(input_data[N*cycle:N+N*cycle+offset],cycle)
    if (cycle >19):
        calc_tve = round(tve.compute_TVE(150, 32, params[0], params[1]), 4)
    else :
        calc_tve = round(tve.compute_TVE(150, 32, params[0], params[1]), 4)
    output.append([cycle, params[0], params[1],params[2],params[3], calc_tve])

with open(f'{root}/Output/output.csv', 'w+', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(output)
