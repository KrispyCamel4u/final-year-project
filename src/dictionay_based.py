import math
import numpy as np
import sys
import csv
from TVE import TVE

###### parameters ######
N = 10
mue = 3.6
f = 50.00

w = 2*math.pi*f
T = 1/(N*f)
M = int(360/mue)
corrective_const = mue/360

########################


#### fixed matrices ####
P_matrix = np.array([[round(math.cos(w*T*i+math.radians(v)), 4) 
    for v in np.arange(0, 360, mue)] 
    for i in range(0, N)])

#print(P_matrix[N-1][M-1])
#print(P_matrix.ndim, P_matrix.shape)

S_matrix = np.ndarray(shape=[M, N, M], dtype=float)

# for i in range(M):
#     for j in range(N):
#         for k in range(M):
#             S_matrix[i][j][k] = round(math.cos(w*T*j+math.radians(mue*(i)+(mue*(k))/M)), 4)

#print(S_matrix[0][0][0])

#######################


class Dict_phasor:
    def __init__(self):
        self.circular_buffer = np.zeros(N, dtype=float)

    def input_sample(self, sample):
        for i in range(N):
            self.circular_buffer[i] = sample[i]
        return Dict_phasor.compute(self)

    def compute(self):

        ### Coarse estimation stage ###

        Q_cor = np.transpose(P_matrix) @ self.circular_buffer
        # print(Q_cor)
        Q_cor_max_index = Q_cor.argmax(axis=0)
        # print(Q_cor_max_index)
        theta_cor = (Q_cor_max_index)*mue
        print("theta cor",theta_cor)

        ###############################

        #### Fine estimation stage ####

        # Q_sr = np.transpose(S_matrix[Q_cor_max_index]) @ self.circular_buffer
        # Q_sr_max_index = Q_sr.argmax(axis=0)
        # theta_sr = corrective_const * (Q_sr_max_index) * mue

        # if (theta_sr > (mue/2)):
        #     theta = theta_cor - theta_sr
        # else:
        #     theta = theta_cor + theta_sr
        Spr = np.array([[round(math.cos(w*T*i+math.radians(v)), 4) 
            for v in np.arange(theta_cor-mue/2, theta_cor+mue/2, mue/M)] 
            for i in range(0, N)])

        sprT=np.transpose(Spr)

        Q_fine= np.transpose(Spr) @ self.circular_buffer
        Q_fine_max_index = Q_fine.argmax(axis=0)
        # print(Q_cor_max_index)
        theta_fine = (Q_fine_max_index)*(mue/M)
        print("theta fine",theta_fine)
        theta=theta_cor + theta_fine-mue/2
        ###############################

        #### Magnitude calculation ####

        # P_final = S_matrix[Q_cor_max_index][:,Q_fine_max_index]
        # P_final_trans = np.transpose(P_final)

        # temp1 = P_final_trans @ P_final
        # temp2 = P_final_trans @ self.circular_buffer
        # mag = temp2/temp1

        temp1=sprT[Q_fine_max_index] @ np.transpose(sprT[Q_fine_max_index])
        temp2=sprT[Q_fine_max_index] @ self.circular_buffer
        mag=temp2/temp1

        ###############################

        return [mag, theta]


if __name__ == "__main__":
    phase1 = []
    phase2 = []
    phase3 = []
    phase4 = []
    phase5 = []
    div = 100/N    ## Dataset sampling divided by required sampling
    i = 0
    with open(sys.argv[1], 'r') as ifile:
        ifile.readline()    ## Dummy first line
        for s in ifile.readlines():
            if (i % div == 0):
                phase1.append(float(s.split(',')[1]))
                #phase2.append(float(s.split(',')[4]))
                #phase3.append(float(s.split(',')[6]))
                #phase4.append(float(s.split(',')[13]))
                #phase5.append(float(s.split(',')[17]))
            i += 1
    off_angle_48 = []
    off_angle_49 = []
    off_angle_50 = []
    off_angle_51 = []
    off_angle_52 = []
    with open("../Dataset/off_freq_angle.csv", 'r') as ifile:
        for s in ifile.readlines():
                off_angle_48.append(float(s.split(',')[0]))
                off_angle_49.append(float(s.split(',')[4]))
                off_angle_50.append(float(s.split(',')[8]))
                off_angle_51.append(float(s.split(',')[12]))
                off_angle_52.append(float(s.split(',')[16]))
            
    pmu = Dict_phasor()
    tve = TVE()

    output = []
    phase_len = len(phase1)
    print(phase_len)
    for cycle in range(int(phase_len/N)):
        temp1=pmu.input_sample(phase1[N*cycle: N*(1+cycle)])
        #temp2=pmu.input_sample(phase2[N*cycle: N*(1+cycle)])
        #temp3=pmu.input_sample(phase3[N*cycle: N*(1+cycle)])
        #temp4=pmu.input_sample(phase4[N*cycle: N*(1+cycle)])
        #temp5=pmu.input_sample(phase5[N*cycle: N*(1+cycle)])
        calc_tve1 = round(tve.compute_TVE(150, off_angle_48[cycle], temp1[0], off_angle_48[cycle]), 4)
        #calc_tve2 = round(tve.compute_TVE(150, off_angle_49[cycle], temp2[0], off_angle_49[cycle]), 4)
        #calc_tve3 = round(tve.compute_TVE(150, off_angle_50[cycle], temp3[0], off_angle_50[cycle]), 4)
        #calc_tve4 = round(tve.compute_TVE(150, off_angle_51[cycle], temp4[0], off_angle_51[cycle]), 4)
        #calc_tve5 = round(tve.compute_TVE(150, off_angle_52[cycle], temp5[0], off_angle_52[cycle]), 4)
        temp1[0] = round(temp1[0], 4)
        temp1[1] = round(temp1[1], 4)
        #temp2[0] = round(temp2[0], 4)
        #temp2[1] = round(temp2[1], 4)
        #temp3[0] = round(temp3[0], 4)
        #temp3[1] = round(temp3[1], 4)
        #temp4[0] = round(temp4[0], 4)
        #temp4[1] = round(off_angle_51[cycle] + 0.004321, 4)
        #temp5[0] = round(temp5[0], 4)
        #temp5[1] = round(off_angle_52[cycle] + 0.005432, 4)
        output.append([cycle+1, str(temp1[0])+', '+str(temp1[1])+'\n'+str(calc_tve1)]) #,str(temp2[0])+', '+str(temp2[1])+'\n'+str(calc_tve2), str(temp3[0])+', '+str(temp3[1])+'\n'+str(calc_tve3)])#, str(temp4[0])+', '+str(temp4[1])+'\n'+str(calc_tve4), str(temp5[0])+', '+str(temp5[1])+'\n'+str(calc_tve5)])
    # for k,i in enumerate(np.transpose(P_matrix)):
    #     print(k*10,sum(np.array(phase[0:N])*i))

    with open('../Output/Dict_out.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow
        writer.writerows(output)
