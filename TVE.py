import math
import sys

class TVE:
    def compute_TVE(self,a_M, a_p, e_M, e_p):
        De_real=e_M*math.cos(math.radians(e_p))
        De_img=e_M*math.sin(math.radians(e_p))
        Da_real=a_M*math.cos(math.radians(a_p))
        Da_img=a_M*math.sin(math.radians(a_p))

        temp1 = (De_real - Da_real)*(De_real - Da_real) 
        temp2 = (De_img - Da_img)*(De_img - Da_img)

        temp3 = (temp1 + temp2)/(Da_real*Da_real + Da_img*Da_img)

        return math.sqrt(temp3) * 100

error_cal=TVE()

with open(sys.argv[1], 'r') as ifile:
    for s in ifile.readlines():
        print(error_cal.compute_TVE(150,10,float(s.split(',')[1]),float(s.split(',')[2])))