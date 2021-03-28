import math,cmath,sys,csv

### success

N = 50  # number of samples per cycle

prev_phase=0
one_time=True
add_360=False

def jexp(angle):
    return complex(math.cos(angle),math.sin(angle))

def dft_modified(wave):
    I_dft=(2/N)*(sum( [wave[k]*jexp(-2*math.pi*k/N) for k in range(N)] ))

    I_even=(2/N)*(sum( [wave[2*k]*jexp(-2*math.pi*2*k/N) for k in range(int(N/2))] ))

    I_odd=(2/N)*(sum( [wave[2*k+1]*jexp(-2*math.pi*(2*k+1)/N) for k in range(int(N/2))] ))

    Kre=(I_even-I_odd).real
    Kimg=(I_even-I_odd).imag

    ## E=e^(-del_t/tau)
    E=(Kimg)/(Kre*math.sin(2*math.pi/N)-Kimg*math.cos(2*math.pi/N))  

    I_dc_dft=(I_even-I_odd)*(1+E*jexp(-2*math.pi/N))/(1-E*jexp(-2*math.pi/N))

    I_fundamental_dft=I_dft-I_dc_dft

    return abs(I_fundamental_dft),math.degrees(cmath.phase(I_fundamental_dft))

def dft_modified1(wave,sample):
    global one_time
    global prev_phase
    global add_360

    I_dft=(2/N)*(sum( [wave[k]*jexp(-2*math.pi*k/N) for k in range(N)] ))

    I_even=(2/N)*(sum( [wave[2*k]*jexp(-2*math.pi*2*k/N) for k in range(int(N/2))] ))

    I_odd=(2/N)*(sum( [wave[2*k+1]*jexp(-2*math.pi*(2*k+1)/N) for k in range(int(N/2))] ))

    Kre=(I_even-I_odd).real
    Kimg=(I_even-I_odd).imag

    ## E=e^(-del_t/tau)
    E=(Kimg)/(Kre*math.sin(2*math.pi/N)-Kimg*math.cos(2*math.pi/N))  

    I_dc_dft=(I_even-I_odd)*(1+E*jexp(-2*math.pi/N))/(1-E*jexp(-2*math.pi/N))

    I_fundamental_dft=I_dft-I_dc_dft

    phase_angle=math.degrees(cmath.phase(I_fundamental_dft))
    # factor=(360/N)*(sample%N)
    # if one_time:
    #     prev_phase=phase_angle
    #     one_time=False
    
    # # print("   ",phase_angle,prev_phase)
    # if abs(phase_angle-prev_phase)>90:
    #     add_360=True
    #     print('add360 True')
    # prev_phase=phase_angle
    
    # if add_360:
    #     phase_angle+=360
    # # print(sample,phase_angle,factor,phase_angle-factor)
    # if (sample%N)==99:
    #     add_360=False
    #     print('add360 false')
    # return abs(I_fundamental_dft),phase_angle-factor
    return abs(I_fundamental_dft),phase_angle

# data gathering from csv file
phase = []
with open(sys.argv[1], 'r') as ifile:
    ifile.readline()
    for s in ifile.readlines():
        phase.append(float(s.split(',')[1]))


# cycle by cycle dft
mCycle = int(len(phase)/N)
output = []
for cycle in range(mCycle):
    params = dft_modified(phase[N*cycle:N+N*cycle])
    output.append([cycle, params[0], params[1]])

# for i in range(len(phase)-N+1):
#     params=dft_modified1(phase[i:i+N],i)
#     output.append([(i+N)*0.0002, params[0], params[1]])

with open('Output/output.csv', 'w+', newline='') as file: 
    writer = csv.writer(file)
    writer.writerows(output)