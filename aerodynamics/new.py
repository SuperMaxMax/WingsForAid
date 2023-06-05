import old as old

#Open file, use proper path
def rep_var(filepath, variable):
    f = open(filepath, "r")

    txt = f.read()
    txtlist = txt.split("\n")

    for lines in txtlist:
        print(lines)
        print("\n")


    txt = txt.replace("x=14","x=15")
    #f = open(filepath, "w")

    #f.write(txt)
    #print(txt)

filepath = r"C:\Users\jarno\OneDrive\Documenten\WingsForAid\aerodynamics\old.py"

rep_var(filepath, 123)
