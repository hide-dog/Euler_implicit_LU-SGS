import pandas
import numpy 
from matplotlib import pyplot

s_in = "time200d-3"
date_1, date_2 ,date_3, date_4 = numpy.loadtxt(s_in, delimiter=' ', skiprows=1, unpack=True)

lst = pandas.read_csv("sod,t=0.2.csv").values.tolist()  # attension please

x=[]
y=[]
for i in range(len(lst)):
    x.append(lst[i][0])
    y.append(lst[i][1])

# axis
x_axis = "length[m]" 
y_axis = "Density[kg/m3]"
pyplot.xlabel(x_axis)
pyplot.ylabel(y_axis)

""" glid line """
pyplot.grid(color="gray", linestyle="dotted", which="both")

pyplot.plot(x, y, 'bo-')
pyplot.plot(date_1, date_3, 'ro-')

""" format adjustment """
pyplot.tight_layout()

pyplot.show()



