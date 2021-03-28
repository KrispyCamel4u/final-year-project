import math, csv, sys

## Beta_max = 4. So the input buffer should atleast be of length 100+4+1 = 105.

nominal_freq = 50

class shift_samples:
    def shfit(wave, est_freq, N):
        out = [0 for i in range(N)]
        out[0] = wave[0]

        for k in range(1, N):
            alpha = k*(nominal_freq - est_freq)/nominal_freq
            alpha, p = shift_samples.temp(abs(alpha), alpha, k)

            alpha2 = alpha*alpha
            b0 = 1 - alpha2
            b1 = 0.5*(1 + alpha2)
            b2 = 0.5*(-1 + alpha2)

            out[k] = b0*wave[p] + b1*wave[p+1] + b2*wave[p-1]
        return out
    


    def temp(beta, alpha, k):
        if(beta < 0.5):
            p = k
        else if(beta < 1.5 && beta > 0.5):
            if (est_freq <= nominal_freq):
                alpha = -1 + alpha
                p = k - 1
            else:
                alpha = 1 - alpha
                p = k + 1
        else if(beta < 2.5 && beta > 1.5):
            if (est_freq <= nominal_freq):
                alpha = -2 + alpha
                p = k - 2
            else:
                alpha = 2 - alpha
                p = k + 2
        else if(beta < 3.5 && beta > 2.5):
            if (est_freq <= nominal_freq):
                alpha = -3 + alpha
                p = k - 3
            else:
                alpha = 3 - alpha
                p = k + 3
        else: 
            if (est_freq <= nominal_freq):
                alpha = -4 + alpha
                p = k - 4
            else:
                alpha = 4 - alpha
                p = k + 4

        return alpha, p
