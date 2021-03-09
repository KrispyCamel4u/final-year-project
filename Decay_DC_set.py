import math
import csv


### parameters ###
N = 100
Rm = 150
f = 50
theta = 10

w = 2*math.pi*f
T = 1/(N*f)
##################


def compute(t, gamma, tou):
    actual = Rm * math.sin(w*T*t + math.radians(theta))

    with_off = actual + gamma * math.exp(-1*(T*t)/tou)

    return [round(actual, 4), round(with_off, 4)]


if __name__ == "__main__":

    ### Fixed tou with varying mag ###
    output_mag_err = []
    output_mag_act = []
    tou = 0.03
    for i in range(1, 15*N+1):
        temp1 = [round(i*T, 4)]
        temp2 = [round(i*T, 4)]
        for j in range(3, 10, 2):
            temp = compute(i, Rm*j/10, tou)
            if (j == 3):
                temp1.append(temp[0])
            temp2.append(temp[1])
        output_mag_err.append(temp2)
        output_mag_act.append(temp1)

    with open("Dataset/DC_offset/My_mag_err.csv", "w", newline='') as ofile:
        writer = csv.writer(ofile)
        # writer.writerow(['time(s)', '30%', '50%', '70%', '90%'])
        writer.writerows(output_mag_err)

    with open("Dataset/DC_offset/My_mag_act.csv", "w", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['time(s)', 'correct_value'])
        writer.writerows(output_mag_act)

    ### End of making dataset for varying mag ###

    ### Fixed mag with varying tou ###
    output_tou_err = []
    output_tou_act = []
    gamma = 150
    for i in range(1, 15*N+1):
        temp1 = [round(i*T, 4)]
        temp2 = [round(i*T, 4)]
        for j in range(1, 11, 3):
            temp = compute(i, gamma, j/100)
            if (j == 1):
                temp1.append(temp[0])
            temp2.append(temp[1])
        output_tou_err.append(temp2)
        output_tou_act.append(temp1)

    with open("Dataset/DC_offset/My_tou_err.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['time(s)', '10ms', '40ms', '70ms', '100ms'])
        writer.writerows(output_tou_err)

    with open("Dataset/DC_offset/My_tou_act.csv", "w+", newline='') as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['time(s)', 'correct_value'])
        writer.writerows(output_tou_act)
    ### End of making dataset for varying tou ###
