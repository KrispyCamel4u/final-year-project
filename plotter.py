import matplotlib.pyplot as plt
import csv

filee = 'Dataset/freq_ramp.csv'
x = []
ys = {}
bypass = 1
bpy = 1
with open(filee) as ifile:
    csv_reader = csv.reader(ifile, delimiter=',')
    for row in csv_reader:
        if bpy:
            bpy = 0
            continue
        x.append(float(row[0]))
        for i, ele in enumerate(row[1:]):
            if bypass:
                ys[i] = []
            ys[i].append(float(ele))
        bypass = 0
plt.figure()
plt.subplot(121)
plt.plot(x, ys[0])
plt.subplot(122)
plt.plot(x, ys[len(ys)-1])
plt.show()
