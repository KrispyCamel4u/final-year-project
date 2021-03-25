import math
import numpy as np
import sys
import matplotlib.pyplot as plt

# wave,wave_abs = [],[]
# x=[]
# with open(sys.argv[1], 'r') as ifile:
#     ifile.readline()
#     for s in ifile.readlines():
#         wave_abs.append(abs(float(s.split(',')[1])))
#         x.append(float(s.split(',')[0]))
#         wave.append(float(s.split(',')[1]))

f0, phase, N = 49.9, 60, 50
dt = 1/(f0*N)
f = 50.3
# dtt=1/(f*N)
dtt = dt

t = np.arange(0, 0.5, dt)
ideal = np.sin(2*np.pi*f0*t+math.radians(phase))
# ideal=np.zeros(len(t))
noise = ideal + 1*np.sin(2*np.pi*2*f0*t+math.radians(phase/2))+0.1 * \
    np.sin(2*np.pi*3*f0*t + math.radians(phase/3)) + 0.1*np.sin(2*np.pi*4*f0*t + math.radians(phase/4)) + \
    0.06*np.sin(2*np.pi*5*f0*t + math.radians(phase/5))+0.06 * \
    np.sin(2*np.pi*6*f0*t + math.radians(phase/6))


def bw_filter_init(fc, fs):
    wd = 2*math.pi*fc
    Ts = 1/fs
    temp = math.tan((wd*Ts) / 2)
    M = 1/(temp**2)
    A = M-math.sqrt(2)+1
    B = 2-2*M
    C = M+1+math.sqrt(2)
    return A, B, C


A, B, C = bw_filter_init(50, 2500)


def filter(noise):
    out = np.zeros(len(noise))
    # out[0] = noise[0]
    for i in range(1, len(noise)):
        out[i] = 0.0000046*(noise[i]+noise[i-1])+0.999999*out[i-1]
        # out[i] = noise[i]+out[i-1]
    return out


def second_order_butter(noise):
    out = np.zeros(len(noise))
    # out[0:2] = noise[0:2]
    for i in range(2, len(noise)):
        # out[i] = 0.0078*(noise[i]+noise[i-2])+0.0156 * \
        #     noise[i-1]-0.7656*out[i-2]+1.734*out[i-1]
        out[i] = (1/C)*(noise[i] + noise[i-2] + 2 *
                        noise[i-1] - A*out[i-2] - B*out[i-1])
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


filtered = filter(noise)
filtered2 = second_order_butter(noise)

ncycle = len(ideal)/N
freq_est = []
for cycle in range(int(ncycle)):
    freq_est.append(find_frequency(filtered2[cycle*N:cycle*N+N]))

print(freq_est)
# x=range(N)
print('helo')
# print(sum((fx*fxx)*dt))
fig, axs = plt.subplots(3, 1)
plt.sca(axs[0])
plt.plot(t, ideal)
plt.plot(t, filtered2)
plt.sca(axs[1])
plt.plot(t, noise)
plt.sca(axs[2])

# plt.sca(axs[3])
plt.plot(t, filtered)

plt.show()
