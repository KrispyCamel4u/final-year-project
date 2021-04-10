import math, csv, sys

## Beta_max = 4. So the input buffer should atleast be of length 100+4+1 = 105.

class shift_samples:
    def __init__(self):
        self.nominal_freq = 50


    def shift(self, wave, est_freq, N):
        out = [0 for i in range(N)]
        wave_len = len(wave)
        out[N-1] = wave[wave_len-1]

        for k in range(wave_len-2, wave_len-N-6, -1):
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
        else:
            rounded_beta = int(round(beta, 0))
            if(freq_bool):
                alpha = -rounded_beta + alpha
                p = k - rounded_beta
            else:
                alpha = rounded_beta - alpha
                p = k + rounded_beta
        #elif(beta <= 1.5 and beta > 0.5):
        #    if (freq_bool):
        #        alpha = -1 + alpha
        #        p = k - 1
        #    else:
        #        alpha = 1 - alpha
        #        p = k + 1
        #elif(beta <= 2.5 and beta > 1.5):
        #    if (freq_bool):
        #        alpha = -2 + alpha
        #        p = k - 2
        #    else:
        #        alpha = 2 - alpha
        #        p = k + 2
        #elif(beta <= 3.5 and beta > 2.5):
        #    if (freq_bool):
        #        alpha = -3 + alpha
        #        p = k - 3
        #    else:
        #        alpha = 3 - alpha
        #        p = k + 3
        #else: 
        #    if (freq_bool):
        #        alpha = -4 + alpha
        #        p = k - 4
        #    else:
        #        alpha = 4 - alpha
        #        p = k + 4

        return alpha, p

    def shift2(self, wave, est_freq, N):
        out = [0 for i in range(N)]
        wave_len = len(wave)
        out[N-1] = wave[wave_len-1]
        freq_bool = est_freq > self.nominal_freq

        for k in range(wave_len-2, 2, -1):
            alpha = (k-wave_len+1)*(est_freq - self.nominal_freq)/est_freq
            print(alpha)
            if(freq_bool):
                alpha = -1-alpha

            b0 = 1 + 1.5*alpha + alpha*alpha
            b1 = -alpha*(2 + alpha)
            b2 = 0.5*alpha*(1 + alpha)

            if(freq_bool):
                out[k-2] = b1*wave[k+1] + b0*wave[k] + b2*wave[k-1]
            else:
                out[k-2] = b0*wave[k] + b1*wave[k-1] + b2*wave[k-2]
                
        return out

    def shift3(self, wave, est_freq, N):
        out = [0 for i in range(N)]
        wave_len = len(wave)
        out[0] = wave[0]
        freq_diff = (est_freq - self.nominal_freq)/est_freq

        for k in range(1, N-1):
            alpha = k*(freq_diff)
            alpha, p = shift_samples.temp(abs(alpha), alpha, k, est_freq<=self.nominal_freq)
            alpha2 = alpha*alpha
            b0 = 1 - alpha2
            b1 = 0.5*(1 + alpha2)
            b2 = 0.5*(-1 + alpha2)
            
            out[k] = b0*wave[p] + b1*wave[p+1] + b2*wave[p-1]

        return out
