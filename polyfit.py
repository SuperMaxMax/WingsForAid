c = [-4.64970233e-08,  1.33281512e-06, -1.39844253e-05,  5.02004763e-05,  #No clue if this is correct
                             1.34340260e-04, -1.05075580e-03, -2.95189891e-03,  1.89082949e-02,
                             8.49297399e-02, -2.79515926e-01, -2.24670348e+00,  2.99288268e+00,
                             5.30704946e+01, -2.00229786e+01, -1.19842392e+03,  6.12640703e+02,
                             2.67539560e+04, -6.27887289e+04, -4.26912912e+05,  3.37435924e+06,
                             -1.17389949e+07,  2.57579109e+07, -3.88663744e+07,  4.13412545e+07,
                             -3.08618888e+07,  1.57188640e+07, -5.14751005e+06,  9.62206912e+05,
                             -7.63817658e+04,  1.02621415e+03]


c = c[::-1]
import numpy as np
x = np.arange(0, 4.5, 0.1)
def polynomial(coefficients, x):
    y = 0
    for i in range(len(coefficients)):
        y+= coefficients[i]*x**i

    return y
        

y = polynomial(c, x)
 
import matplotlib.pyplot as plt

plt.plot(x, y)
plt.show()
