import math
import numpy as np
import sys
import csv

###### parameters ######
N = 50
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
    phase = []
    div = 100/N    ## Dataset sampling divided by required sampling
    i = 0
    with open(sys.argv[1], 'r') as ifile:
        ifile.readline()    ## Dummy first line
        for s in ifile.readlines():
            if (i % div == 0):
                phase.append(float(s.split(',')[-1]))
            i += 1
    pmu = Dict_phasor()

    output = []
    phase_len = len(phase)
    print(phase_len)
    for cycle in range(int(phase_len/N)):
        temp=pmu.input_sample(phase[N*cycle: N*(1+cycle)])
        output.append([cycle+1, temp[0],temp[1]])
    # for k,i in enumerate(np.transpose(P_matrix)):
    #     print(k*10,sum(np.array(phase[0:N])*i))

    with open('../Output/Dict_output.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(output)
