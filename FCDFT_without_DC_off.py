import sys
import math
import csv
import time

##### Parameters #####
N = 100     ## No. of samples

class remove_DC_off:
    def __init__(self):
        self.circular_buffer = [0 for i in range(N)]

    def input_sample(self, sample):
        self.circular_buffer = sample
        return remove_DC_off.compute(self)

    def compute(self):
        sum_of_samples = sum(self.circular_buffer)
        print ("The sum of the samples are", sum_of_samples)
        if (sum_of_samples <  0.05):
            return [0, 0, 0]
        else :
            tou1 = 0.95
            tou2 = 0.99

            iterations = 0
            while(abs(tou1 - tou2) > 0.0004):
                B1 = (sum_of_samples*(tou1 - 1))/(pow(tou1, N+1) - tou1)     ## Calculating B part in the DC offset.
                B2 = (sum_of_samples*(tou2 - 1))/(pow(tou2, N+1) - tou2)     ## Calculating B part in the DC offset.

                err_samples1 = []
                err_samples2 = []
                for k in range(1, N+1):
                    err_samples2.append(self.circular_buffer[k-1] - B2*pow(tou2, k))
                    err_samples1.append(self.circular_buffer[k-1] - B1*pow(tou1, k))
                #print(dft(err_samples1), dft(err_samples2))
                f_tou1 = sum(err_samples1)     ## sum of error samples with tou1        
                f_tou2 = sum(err_samples2)     ## sum of error samples with tou2
                #print("f_tou1:", f_tou1, "f_tou2:", f_tou2)
                '''
                DC_off1 = 0
                DC_off2 = 0
                for k in range(1, N+1):
                    DC_off1 += B1*pow(tou1,k) 
                    DC_off2 += B2*pow(tou2,k) 
                
                f_tou1 = sum_of_samples - B1*(pow(tou1, N+1) - 1)/(tou1 - 1)    ## sum of error samples with tou1
                f_tou2 = sum_of_samples - B2*(pow(tou2, N+1) - 1)/(tou2 - 1)    ## sum of error samples with tou2
                '''
                print(f_tou1, f_tou2)
                print(B1, tou1, B2, tou2)
                try:
                    tou_new = tou2 - f_tou2*((tou1 - tou2)/(f_tou1 - f_tou2))
                    tou1 = tou2
                    tou2 = tou_new
                except ZeroDivisionError:
                    tou1 = tou2
                    print("Float division error encountered")
                    print(tou1, tou2)

                print("tou1:", tou1, "tou2:", tou2)# "B1:", B1, "B2:", B2)

                iterations += 1

            print()
            B = (sum_of_samples*(tou2 - 1))/(pow(tou2, N+1) - 1)     ## Calculating B part in the DC offset.
            #print("B: ", B, " tou1: ", tou1, "tou2:", tou2, " iters: ", iterations)
            return [B, tou2, iterations]


def dft(wave):
    x1 = (2/N)*(sum((wave[i-1]*math.cos(2*math.pi*i/N) for i in range(1, N+1))))
    x2 = (2/N)*(sum((wave[i-1]*math.sin(2*math.pi*i/N) for i in range(1, N+1))))
    #print(x1, x2)
    mag = math.sqrt(pow(x1, 2)+pow(x2, 2))
    phaseAngle = math.atan(x1/x2)
    return mag, math.degrees(phaseAngle)

if __name__ == "__main__":
    phase = []
    with open(sys.argv[1], 'r') as ifile:
        ifile.readline()    ## Heading first line read.
        for s in ifile.readlines():
            phase.append(float(s.split(',')[1]))

    mCycle = int(len(phase)/100)

    pmu = remove_DC_off()
    
    DC_off_removed = []
    output = []
    for cycle in range(int(len(phase)/N)):
        params = pmu.input_sample(phase[N*cycle:N+N*cycle])
        for k in range(N):
            #DC_off_removed.append(phase[N*cycle+k] - params[0]*math.exp((N*cycle+(k+1))*(-0.0002)/params[1]))
            DC_off_removed.append(phase[N*cycle+k] - params[0]*pow(params[1], k+1))
            
        output.append([cycle+1, dft(DC_off_removed[N*cycle:N+N*cycle]), params[2]])

    with open('output/DFT_output.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(output)

