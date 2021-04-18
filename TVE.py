import math
import sys

mag=150
phase=32

class TVE:
    def compute_TVE(self, a_M, a_p, e_M, e_p):
        De_real=e_M*math.cos(math.radians(e_p))
        De_img=e_M*math.sin(math.radians(e_p))
        Da_real=a_M*math.cos(math.radians(a_p))
        Da_img=a_M*math.sin(math.radians(a_p))

        temp1 = (De_real - Da_real)**2 
        temp2 = (De_img - Da_img)**2

        temp3 = (temp1 + temp2)/(Da_real*Da_real + Da_img*Da_img)

        return math.sqrt(temp3) * 100

#error_cal=TVE()
#
#with open(sys.argv[1], 'r') as ifile:
#    for s in ifile.readlines():
#        if (int(s.split(',')[0]) > 9):
#            print(error_cal.compute_TVE(mag,float(s.split(',')[2]),float(s.split(',')[1]),float(s.split(',')[2])))
