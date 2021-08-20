from ignm.source.GasPhaseSolver import solver
import matplotlib.pyplot as plt


solution = solver('../data/example_input.yaml','../data/ts.csv','../data/yin.csv','../data/j0.csv' )

plt.plot(solution['time'], solution[solution.columns[-1]])
plt.xlabel('time (s)')
plt.ylabel('Temperature (K)')
plt.show()


    
    