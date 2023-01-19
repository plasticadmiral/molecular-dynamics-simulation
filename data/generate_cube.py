import os
import sys
import glob
import numpy as np


# na1 = 3 
# na2 = 3 
# na3 = 3 

#distances
a1 = np.array([1.0, 0.0, 0.0])
a2 = np.array([0.0, 1.0, 0.0])
a3 = np.array([0.0, 0.0, 1.0])

def single_cube(size, file = "_custom_cube.xyz"):
    
    na1 = 3 + size
    na2 = 3 + size
    na3 = 3 + size

    path = "cubes/" + file

    with open(path, 'w') as f:
        f.write(str(na1*na2*na3) + "\n")
        f.write(" " + "\n")
        # print(na1*na2*na3)
        # print(" ")
        for n1 in range(na1):
            for n2 in range(na2):
                for n3 in range(na3):
                    p = n1*a1 + n2*a2 +n3*a3
                    f.write("py" + "    " +str(p[0]) + "    " + str(p[1]) + 
                    "   " + str(p[2]) + "   " + "0.0" + "   " + "0.0" + "   " 
                    + "0.0" + "\n")
                    # print("py" + "   " +str(p[0]) + "  " + str(p[1]) + "    " 
                    # + str(p[2]))

def multiple_cubes(initial, final, step):

    files = glob.glob("cubes/cube_*")
    for f in files:
        os.remove(f)

    for i in range(initial, final, step):
        single_cube(i, "cubes/cube_" + str(i) + ".xyz")


#python generate_cube.py 1 10 
#generates cube with size 40 and writes to cubes/_cube_custom.xyz

#python generate_cube.py 2 10 30 3                              
#generates cubes with sizes 10 13 16 19 ... and stores them in cubes folder

if __name__ == "__main__":
    mode = int(sys.argv[1])
        

    if(mode == 1):
        print("generating single cube to cubes/ folder")
        size = int(sys.argv[2])
        single_cube(size)

    if(mode == 2):
        print("generating multiple cubes to cubes/ folder")
        intial = int(sys.argv[2])   
        final = int(sys.argv[3])
        steps = int(sys.argv[4])
        multiple_cubes(intial, final, steps)