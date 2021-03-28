import math
import numpy as np
import sys
import matplotlib.pyplot as plt

############## my filter same as butter worth

# wave,wave_abs = [],[]
# x=[]
# with open(sys.argv[1], 'r') as ifile:
#     ifile.readline()
#     for s in ifile.readlines():
#         wave_abs.append(abs(float(s.split(',')[1])))
#         x.append(float(s.split(',')[0]))
#         wave.append(float(s.split(',')[1]))

f0, phase, N = 48, 30, 50
dt = 1/(50*N)
f = 50.3
# dtt=1/(f*N)
dtt = dt
fs=1/dt

t = np.arange(0, 2, dt)
ideal = np.sin(2*np.pi*f0*t+math.radians(phase))
# ideal=np.zeros(len(t))
noise = ideal + 0.06*np.sin(2*np.pi*2*f0*t+math.radians(phase/2))+0.06                                                                                                                                                                                                                                                                                                                                                                               * \
    np.sin(2*np.pi*3*f0*t + math.radians(phase/3)) + 0.06*np.sin(2*np.pi*4*f0*t + math.radians(phase/4)) + \
    0.06*np.sin(2*np.pi*5*f0*t + math.radians(phase/5))+0.06 * \
    np.sin(2*np.pi*6*f0*t + math.radians(phase/6))


def bw_filter_init(fc, fs,order=2):
    wd = 2*math.pi*fc
    Ts = 1/fs
    temp = math.tan((wd*Ts) / 2)
    if order==1:
        M = 1/temp
        B=1-M
        C=M+1
        return B,C
    elif order==2:
        M = 1/(temp**2)
        A = M-math.sqrt(2)+1
        B = 2-2*M
        C = M+1+math.sqrt(2)
        return A, B, C

B1, C1 = bw_filter_init(50, 1/dt,1)
A2, B2, C2 = bw_filter_init(50, 1/dt)

def first_order_butter(noise):
    out = np.zeros(len(noise))
    # out[0] = noise[0]
    for i in range(1, len(noise)):
        out[i] = (1/C1)*(noise[i] + noise[i-1] - B1*out[i-1])
        # out[i] = noise[i]+out[i-1]
    return out


def second_order_butter(noise):
    out = np.zeros(len(noise))
    # out[0:2] = noise[0:2]
    for i in range(2, len(noise)):
        # out[i] = 0.0078*(noise[i]+noise[i-2])+0.0156 * \
        #     noise[i-1]-0.7656*out[i-2]+1.734*out[i-1]
        out[i] = (1/C2)*((1/10)*(noise[i] + noise[i-2] + 2 *
                        noise[i-1]) - A2*out[i-2] - B2*out[i-1])
    return out


def find_frequency(wave):
    wave_abs = list(np.abs(wave))
    wave = list(wave)
    mxIndex = wave.index(max(wave))
    mnIndex = wave.index(min(wave))
    T = (abs(mxIndex-mnIndex))*dt*2
    # print(wave_abs[:12])
    # mn1 = wave_abs.index(min(wave_abs[:50]))
    # mx1 = wave_abs.index(max(wave_abs[:50]))
    # mn2 = wave_abs.index(min(wave_abs[50:]))
    # mx2 = wave_abs.index(max(wave_abs[50:]))
    # print('first', 1/(abs(mx1-mn1)*dt*4))
    # print('second', 1/(abs(mx2-mn2)*dt*4))
    return 1/T

def parabola_aprox(wave,cycle,sign=-1):
    wave=list(wave)
    # print(wave)
    wave_sort=sorted(wave)
    y1=wave_sort[sign*1]
    inm1=wave.index(y1)
    if inm1==0:
        inm2=1
        inm3=2
    elif inm1==N-1:
        inm2=inm1-1
        inm3=inm2-1
    else:
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
    tm2=parabola_aprox(wave,cycle,1)
    T=2*abs(tm1-tm2)
    # print(1/T,'Hz')
    # print("its tm1",tm1)
    return 1/T



filtered = first_order_butter(noise)
filtered2 = second_order_butter(noise)

ncycle = len(ideal)/N
freq_est2 = []
freq_est1 = []
for cycle in range(int(ncycle)):
    freq_est2.append(pricise_freq(filtered2[cycle*N:cycle*N+N],cycle))
    print(cycle,"error",abs(freq_est2[cycle]-f0))

# print(freq_est2)
# print(freq_est1)
# x=range(N)

plt.plot(range(int(ncycle)),freq_est2)
plt.show()

# print('helo')
# # print(sum((fx*fxx)*dt))
# fig, axs = plt.subplots(3, 1)
# plt.sca(axs[0])
# plt.plot(t, ideal)
# plt.plot(t, filtered2)
# plt.sca(axs[1])
# plt.plot(t, noise)
# plt.sca(axs[2])

# # plt.sca(axs[3])
# plt.plot(t, filtered)

# plt.show()
