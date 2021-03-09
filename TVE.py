import math

class TVE:
    def compute_TVE(De_real, De_img, Da_real, Da_img):
        temp1 = (De_real - Da_real)*(De_real - Da_real) 
        temp2 = (De_img - Da_img)*(De_img - Da_img)

        temp3 = (temp1 + temp2)/(Da_real*Da_real + Da_img*Da_img)

        return math.sqrt(temp3) * 100
