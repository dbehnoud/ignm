import ignm.source.GasPhaseSolver
import matplotlib.pyplot as plt

def simulation():
    y, t = ignm.source.GasPhaseSolver('ignm/data/example_input.yaml','ignm/data/ts.csv',\
            'ignm/data/yin.csv','ignm/data/j0.csv' )
    
    
    