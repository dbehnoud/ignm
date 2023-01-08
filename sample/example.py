from ignm.gasPhaseSolver import solve
import matplotlib.pyplot as plt

input_filepath ='../data/input.yaml'
ts_filepath = '../data/ts.csv'
yin_filepath = '../data/yin.csv'
j0_filepath = '../data/j0.csv'

solution = solve(input_filepath,ts_filepath,yin_filepath,j0_filepath)

plt.plot(solution['time'], solution[solution.columns[-1]])
plt.xlabel('time (s)')
plt.ylabel('Temperature (K)')
plt.show()


    
    