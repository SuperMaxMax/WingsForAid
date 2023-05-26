
import csv

### options of airfoil
### [2412, 4415, 23012, 23015, 63215, 64415, "clarky", "usa35b"]
airfoil = "usa35b"

### CL/ALPHA STUFF ###
file_name = "cl-alpha-" + str(airfoil)

file_path = "/Users/janvonken/Desktop/TU/3rd Year/DSE/" + file_name + ".csv"

cl_max = 0

with open(file_path) as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    linecount = 0
    for row in reader:
        if row:    
            if linecount == 0:
                linecount += 1
            else:
                alpha = float(row[0])
                cl = float(row[1])

                if airfoil != "usa35b":
                    if alpha == 0:
                        cl_0 = cl
                    elif alpha == -4: 
                        cl_min5 = cl
                    elif alpha == 4:
                        cl_plus5 = cl
                else: 
                    cl_alpha = 0.10777
                    cl_0 = 0.62

                if airfoil == "63215":
                    cl_max = 1.175
                    alpha_stall = 9.5
                elif airfoil == "64415":
                    cl_max = 1.225
                    alpha_stall = 11.5
                else:
                    if cl > cl_max: 
                        cl_max = cl
                        alpha_stall = alpha

if airfoil != "usa35b":
    cl_alpha = (cl_plus5 - cl_min5)/8


### CM/ALPHA STUFF ###
file_name = "cm-alpha-" + str(airfoil)

file_path = "/Users/janvonken/Desktop/TU/3rd Year/DSE/" + file_name + ".csv"

with open(file_path) as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    linecount = 0
    for row in reader:
        if row:    
            if linecount == 0:
                linecount += 1
            else:
                alpha = float(row[0])
                cm = float(row[1])
                if airfoil == "usa35b":
                    cm_0 = -0.0875
                    cm_alpha = 0.00731
                else:
                    if alpha == 0:
                        cm_0 = cm

                if airfoil ==  2412:
                    cm_alpha = 0.00578
                elif airfoil == 4415: 
                    cm_alpha = 0.00748
                elif airfoil == 23012: 
                    cm_alpha = 0.00490
                elif airfoil == 23015: 
                    cm_alpha = 0.00645
                elif airfoil == 63215: 
                    cm_alpha = 0.00886
                elif airfoil == 64415: 
                    cm_alpha = 0.0105
                elif airfoil == "clarky": 
                    cm_alpha = 0.00627
                

### CL/CD STUFF ###
file_name = "cl-cd-" + str(airfoil)

file_path = "/Users/janvonken/Desktop/TU/3rd Year/DSE/" + file_name + ".csv"

cl_cd = []
cl32_cd = []
cl12_cd = []

with open(file_path) as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    linecount = 0
    for row in reader:
        if row:    
            if linecount == 0:
                linecount += 1
            else:
                cd = float(row[0])
                cl = float(row[1])
                cl_cd.append(cl/cd)
                if cl > 0:
                    cl32_cd.append(cl**1.5 / cd)
                    cl12_cd.append(cl**0.5 / cd)
                
                if airfoil == "usa35b":
                    cd0 = 0.0085055
                else:
                    if cl == cl_0:
                        cd0 = cd
                


print("cl_max = ", cl_max)
print("alpha_stall = ", alpha_stall)
print("cl_alpha = ", cl_alpha)
print("cl_0 = ", cl_0)
print("Cl/Cd max = ", max(cl_cd))
print("Cl^1.5/Cd max = ", max(cl32_cd))
print("Cl^0.5/Cd max = ", max(cl12_cd))
print("c_m0 = ", cm_0)
print("c_m_alpha = ", cm_alpha)
print("c_d0 = ", cd0)
