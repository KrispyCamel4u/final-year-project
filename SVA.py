import math, csv, sys

## Beta_max = 4. So the input buffer should atleast be of length 100+4+1 = 105.

class shift_samples:
    print('yes')
    def __init__(self):
        self.nominal_freq = 50

    def shift(self, wave, est_freq, N):
        out = [0 for i in range(N)]
        out[0] = wave[0]
        del_f=(est_freq - self.nominal_freq)/est_freq
        for k in range(1,N):
            alpha = k*del_f
            # alpha, p = shift_samples.temp(abs(alpha), alpha, k, est_freq<=self.nominal_freq)

            b0=1-alpha**2
            b1=((alpha-1)*alpha)/2
            b2=((alpha+1)*alpha)/2
            
            out[k] = b0*wave[k] + b1*wave[k+1] + b2*wave[k-1]
        return out
    
    # def temp(beta, alpha, k, freq_bool):
    #     if(beta < 0.5):
    #         p = k
    #     elif(beta <= 1.5 and beta > 0.5):
    #         if (freq_bool):
    #             alpha = 1 + alpha
    #             p = k + 1
    #         else:
    #             alpha = alpha -1
    #             p = k -1
    #     elif(beta <= 2.5 and beta > 1.5):
    #         if (freq_bool):
    #             alpha = 2 + alpha
    #             p = k + 2
    #         else:
    #             alpha =  alpha -2
    #             p = k - 2
    #     elif(beta <= 3.5 and beta > 2.5):
    #         if (freq_bool):
    #             alpha = 3 + alpha
    #             p = k + 3
    #         else:
    #             alpha = alpha - 3
    #             p = k - 3
    #     else: 
    #         if (freq_bool):
    #             alpha = 4 + alpha
    #             p = k + 4
    #         else:
    #             alpha =alpha-4
    #             p = k -4

    #     return alpha, p

#if __name__=="__main__":
#    #data gathering from csv file
#    phase = []
#    with open(sys.argv[1], 'r') as ifile:
#        ifile.readline()
#        for s in ifile.readlines():
#            phase.append(float(s.split(',')[1]))
#    ideal = []
#    with open("Dataset/normal.csv", 'r') as ifile:
#        ifile.readline()
#        for s in ifile.readlines():
#            ideal.append(float(s.split(',')[1]))
#
#    N=48
#    shifter=shift_samples()
#    shifted_samples=shifter.shift(phase[:N+1],48,N)+shifter.shift(phase[N:2*N+1],48,N)
#    import matplotlib.pyplot as plt
#    plt.plot(phase[:2*N])
#    plt.plot(shifted_samples)
#    plt.plot(ideal[:2*N])
#
#    plt.show()



