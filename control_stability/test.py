import sympy as sp
import numpy as np
import matplotlib.pyplot as plt 
# Definieer de symbolen x en y
x, y = sp.symbols('x y')

# Bepaal de coördinaten van punt B
B = (1, 0)

# Definieer de vergelijking van de cirkel
circle_eq = (x - 4)**2 + (y)**2 - 13

# Bepaal de afgeleide van de cirkelvergelijking
dy_dx = sp.diff(circle_eq, y) / sp.diff(circle_eq, x)

# Bepaal de helling van de raaklijnen in het punt B
m = dy_dx.subs([(x, B[0]), (y, B[1])])
print(m)
# Bepaal de vergelijkingen van de raaklijnen n1 en n2
y1 = - B[1] - m * (x - B[0])
n2_eq = (y - B[1]) + m * (x - B[0])

# Vereenvoudig de vergelijkingen
#1_eq = sp.simplify(n1_eq)
n2_eq = sp.simplify(n2_eq)
y_1_range = []
y_2_range = []

x_range = np.arange(0, 5, 0.1)
for i in x_range:
    y_1 = B[1] - m * (i - B[0])
    y_2 = B[1] + m * (i - B[0])
    y_1_range.append(y_1)
    y_2_range.append(y_2)

plt.plot(x_range, y_1_range)
plt.show()
