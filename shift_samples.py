import math, csv, sys

## Beta_max = 4. So the input buffer should atleast be of length 100+4+1 = 105.

class shift_samples:
    def __init__(self):
        self.nominal_freq = 50


    def shift(self, wave, est_freq, N):
        out = [0 for i in range(N)]
        wave_len = len(wave)
        out[N-1] = wave[wave_len-1]

        for k in range(wave_len-2, wave_len-N-3, -1):
            alpha = (k-wave_len+1)*(est_freq - self.nominal_freq)/est_freq
            alpha, p = shift_samples.temp(abs(alpha), alpha, k, est_freq<=self.nominal_freq)

            alpha2 = alpha*alpha
            b0 = 1 - alpha2
            b1 = 0.5*(1 + alpha2)
            b2 = 0.5*(-1 + alpha2)

            out[k-5] = b0*wave[p] + b1*wave[p+1] + b2*wave[p-1]
        return out
    
    def temp(beta, alpha, k, freq_bool):
        if(beta < 0.5):
            p = k
        elif(beta <= 1.5 and beta > 0.5):
            if (freq_bool):
                alpha = -1 + alpha
                p = k - 1
            else:
                alpha = 1 - alpha
                p = k + 1
        elif(beta <= 2.5 and beta > 1.5):
            if (freq_bool):
                alpha = -2 + alpha
                p = k - 2
            else:
                alpha = 2 - alpha
                p = k + 2
        elif(beta <= 3.5 and beta > 2.5):
            if (freq_bool):
                alpha = -3 + alpha
                p = k - 3
            else:
                alpha = 3 - alpha
                p = k + 3
        else: 
            if (freq_bool):
                alpha = -4 + alpha
                p = k - 4
            else:
                alpha = 4 - alpha
                p = k + 4

        return alpha, p
