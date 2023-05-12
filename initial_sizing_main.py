
import Class_I_weight_estimation as c1
import geometry_determination as geo

class initial_sizing:
    def __init__(self):
        pass

    def iteration(self):
        run = True
        while run:
            c1.weight.iteration()



if __name__ == "__main__":
    weight = initial_sizing()
    weight.iteration()

