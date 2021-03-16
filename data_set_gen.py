import math
import csv
import numpy as np

### parameters ###
N = 100
Rm = 150
f0 = 50
theta = 30

w = 2*math.pi*f0
del_T = 1/(N*f0)

num_cycle = 10
# decaying DC
# gamma 30% to 90%
# tau 10 to 100ms

# Harmonics
# R2 1.33 to 6.667%
# Rk+1 = 0.2*Rk
# theta_k=theta/k

# modulation
ka = 0.1   # magnitude modulation index
kx = 0.1   # phase modulation index

##################


def compute_ddc(n, gamma, tou):
    actual = Rm * math.cos(w*del_T*n + math.radians(theta))

    with_off = actual + gamma * math.exp(-1*(del_T*n)/tou)

    return [round(actual, 4), round(with_off, 4)]


## normal
def normal():
    output_mag_err = []
    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        temp = compute_ddc(i, 0, 1)
        temp2.append(temp[0])
        output_mag_err.append(temp2)

    with open("Dataset/normal.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['time(s)', 'signal'])
        writer.writerows(output_mag_err)

def ddc_offset():
    ### Fixed tou with varying mag ###
    output_mag_err = []
    tou = 0.03
    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in range(3, 10, 2):
            temp = compute_ddc(i, Rm*j/10, tou)
            temp2.append(temp[1])
        output_mag_err.append(temp2)

    with open("Dataset/DC_offset/ddc_offset_mag.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['time(s)', '30%', '50%', '70%', '90%'])
        writer.writerows(output_mag_err)

    ### End of making dataset for varying mag ###

    ### Fixed mag with varying tou ###
    output_tou_err = []
    gamma = 150*0.9
    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in range(1, 11, 3):
            temp = compute_ddc(i, gamma, j/100)
            temp2.append(temp[1])
        output_tou_err.append(temp2)

    with open("Dataset/DC_offset/ddc_offset_tau.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['time(s)', '10ms', '40ms', '70ms', '100ms'])
        writer.writerows(output_tou_err)

    ### End of making dataset for varying tou ###

# harmonics


def compute_harmonics(n, R2):

    actual = Rm * math.cos(w*del_T*n + math.radians(theta))
    with_har = actual
    mag_k = R2
    theta_k = theta
    for i in range(2, 11):
        with_har += mag_k * math.cos(i*w*del_T*n + math.radians(theta_k))
        mag_k *= 0.2
        theta_k /= i

    return [round(actual, 4), round(with_har, 4)]


def harmonics():
    output_mag_err = []
    header = ['time(s)']

    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in np.arange(1.33, 6.66, 1):
            temp = compute_harmonics(i, Rm*j/100)
            temp2.append(temp[1])
            if i == 1:
                header.append(str(j)+'%')
        output_mag_err.append(temp2)

    with open("Dataset/with_harmonics.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(header)
        writer.writerows(output_mag_err)

    ### End of making dataset for varying harmonics magnitude ###

# off nominal frequency


def compute_off_nominal(n, freq):
    return round(Rm * math.cos(2*math.pi*freq*del_T*n + math.radians(theta)), 4)


off_nominal_step = 0.25


def off_nominal_freq():
    output_mag_err = []
    header = ['time(s)']

    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in np.arange(f0-2, f0+2+off_nominal_step, off_nominal_step):
            temp = compute_off_nominal(i, j)
            temp2.append(temp)
            if i == 1:
                header.append(str(j)+'Hz')
        output_mag_err.append(temp2)

    with open("Dataset/off_nominal_freq.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(header)
        writer.writerows(output_mag_err)

######### for magnitude modulation #####


def compute_modulation(n, fm, kx, ka):
    return round(Rm*(1+kx*math.cos(2*math.pi*fm*del_T*n)) * math.cos(w*del_T*n + ka*math.cos(2*math.pi*fm*del_T*n) + math.radians(theta)), 4)


def modulation():
    output_mag_err = []
    header = ['time(s)']

    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in np.arange(f0-2, f0+2+off_nominal_step, off_nominal_step):   # modulation freq
            temp = compute_modulation(i, j, kx, 0)
            temp2.append(temp)
            if i == 1:
                header.append(str(j)+'Hz')
        output_mag_err.append(temp2)

    with open("Dataset/magnitude_modulation.csv", "w", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(header)
        writer.writerows(output_mag_err)

    output_mag_err = []
    header = ['time(s)']

    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in np.arange(f0-2, f0+2+off_nominal_step, off_nominal_step):   # modulation freq
            temp = compute_modulation(i, j, 0, ka)
            temp2.append(temp)
            if i == 1:
                header.append(str(j)+'Hz')
        output_mag_err.append(temp2)

    with open("Dataset/phase_modulation.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(header)
        writer.writerows(output_mag_err)

#### step changes   ####


def compute_step(n, step_start_n, kx, ka, sign=1):
    if n >= step_start_n:
        step_on = sign
    else:
        step_on = 0
    return round(Rm*(1+kx*step_on) * math.cos(w*del_T*n + ka*step_on + math.radians(theta)), 4)

# by default produces +ve step but using sign var in compute step can be changed


def step():
    output_mag_err = []
    header = ['time(s)']

    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in range(300, 420, 20):  # diffent time instants
            temp = compute_step(i, j, kx, 0, 1)
            temp2.append(temp)
            if i == 1:
                header.append(str(j)+'sample')
        output_mag_err.append(temp2)

    with open("Dataset/magnitude_step.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(header)
        writer.writerows(output_mag_err)

    output_mag_err = []
    header = ['time(s)']

    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in range(300, 420, 20):  # diffent time instants
            temp = compute_step(i, j, 0, math.pi/18)
            temp2.append(temp)
            if i == 1:
                header.append(str(j)+'sample')
        output_mag_err.append(temp2)

    with open("Dataset/phase_step.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(header)
        writer.writerows(output_mag_err)

# freq ramp


def compute_freq_ramp(n, ramp_start_n, ramp_rate, sign=1):
    if n >= ramp_start_n:
        ramp_on = sign
    else:
        ramp_on = 0
    return round(Rm * math.cos(w*del_T*n + ramp_on*math.pi*ramp_rate*(del_T*n)**2 + math.radians(theta)), 4)


def freq_ramp():
    output_mag_err = []
    header = ['time(s)']

    for i in range(0, num_cycle*N+1):
        temp2 = [round(i*del_T, 4)]
        for j in range(300, 420, 20):  # diffent time instants
            temp = compute_freq_ramp(i, j, 1)
            temp2.append(temp)
            if i == 1:
                header.append(str(j)+'sample')
        output_mag_err.append(temp2)

    with open("Dataset/freq_ramp.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(header)
        writer.writerows(output_mag_err)

normal()
ddc_offset()
harmonics()
off_nominal_freq()
modulation()
step()
freq_ramp()
