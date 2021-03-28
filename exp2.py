import math
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz,filtfilt
from scipy import signal

###### using scipy filter

# wave,wave_abs = [],[]
# x=[]
# with open(sys.argv[1], 'r') as ifile:
#     ifile.readline()
#     for s in ifile.readlines():
#         wave_abs.append(abs(float(s.split(',')[1])))
#         x.append(float(s.split(',')[0]))
#         wave.append(float(s.split(',')[1]))

f0, phase, N = 48, 60, 50
dt = 1/(50*N)
f = 50.3
# dtt=1/(f*N)
dtt = dt
cutoff=50
fs=1/dt
order=8
offset=5
Ah=0.06

t = np.arange(0, 2, dt)
ideal = np.sin(2*np.pi*f0*t+math.radians(phase))
# ideal=np.zeros(len(t))
noise = ideal + Ah*np.sin(2*np.pi*2*f0*t+math.radians(phase/2))+Ah*\
     np.sin(2*np.pi*3*f0*t + math.radians(phase/3)) + Ah*np.sin(2*np.pi*4*f0*t + math.radians(phase/4)) + \
    Ah*np.sin(2*np.pi*5*f0*t + math.radians(phase/5))+Ah * \
    np.sin(2*np.pi*6*f0*t + math.radians(phase/6)) + Ah*np.sin(2*np.pi*7*f0*t + math.radians(phase/7)) +\
        Ah*np.sin(2*np.pi*8*f0*t + math.radians(phase/8))


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    # y = filtfilt(b, a, data)
    return y

# Get the filter coefficients so we can check its frequency response.
b, a = butter_lowpass(cutoff, fs, order)


def parabola_aprox(wave,cycle,sign=-1):
    wave=list(wave)
    # print(wave)
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
    # print(x1,x2,x3,y1,y2,y3)
    k=( math.pow(x2,2)*(y1-y3) + math.pow(x1,2)*(y3-y2) + math.pow(x3,2)*(y2-y1))/( x2*(y1-y3) + x1*(y3-y2) + x3*(y2-y1) )
    return k/2


def pricise_freq(wave,cycle):
    
    tm1=parabola_aprox(wave,cycle,-1)
    tm2=parabola_aprox(wave,cycle,0)
    T=2*abs(tm1-tm2)
    # print(1/T,'Hz')
    # print("its T",tm1,tm2)
    return 1/T


# Filter the data, and plot both the original and filtered signals.
# bb, aa =  signal.cheby1(5, 1, 0.05)
filtered = butter_lowpass_filter(noise, cutoff, fs, order)
# filtered = lfilter(bb, aa, noise)

w, h = freqz(b=bb, a=aa)
x = w * fs * 1.0 / (2 * np.pi)
y = 20 * np.log10(abs(h))
plt.figure(figsize=(10,5))
plt.semilogx(x, y)
plt.ylabel('Amplitude [dB]')
plt.xlabel('Frequency [Hz]')
plt.title('Frequency response')
plt.grid(which='both', linestyle='-', color='grey')
plt.xticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000], ["20", "50", "100", "200", "500", "1K", "2K", "5K", "10K", "20K"])
plt.show()


ncycle = len(ideal)/N
freq_est2 = []
freq_est1 = []
error_list=[]
maxFE=0
for cycle in range(int(ncycle)):
    freq_est2.append(pricise_freq(filtered[cycle*N:cycle*N+N +offset],cycle))
    error=abs(freq_est2[cycle]-f0)
    print(cycle,"error",error)
    error_list.append(error)
    # if error>0.01 and cycle>6:
    #     plt.plot(t[cycle*N:cycle*N+N +offset],filtered[cycle*N:cycle*N+N +offset])
    #     # plt.show()
    if cycle>15 and cycle<96 and error>maxFE:
        maxFE=error

# print(freq_est2)
# print(freq_est1)
# x=range(N)
print("max FE",maxFE)
plt.plot(range(int(ncycle)),freq_est2)
plt.plot(range(int(ncycle)),error_list)
plt.show()

# print('helo')
# # print(sum((fx*fxx)*dt))
# fig, axs = plt.subplots(3, 1)
# plt.sca(axs[0])
# plt.plot(t, ideal)
# plt.plot(t, filtered2)
# plt.plot(t, filtered)
# plt.sca(axs[1])
# plt.plot(t, noise)
# plt.sca(axs[2])

# # plt.sca(axs[3])
# plt.plot(t, filtered)

# plt.show()
