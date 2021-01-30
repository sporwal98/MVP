Simulation of Ising Model using Monte Carlo Method:
Glauber Dynamics or Kawasaki Dynamics
By Sahaj Porwal s1705173

To run the simulation
Open command prompt/terminal
Navigate to folder containing the mcising.py and genTPlots files
Open a python terminal/command prompt using:
$ python 
In the python terminal enter:
$ from mcising import *

For a single temperature, enter:
$ main(model = <model>, T = <temperature>, size = <size>)
Replacing <model>, <temperature> and <size> by the mechanism to use, 
the temperature to simulate and the size of the spin matrix.

The data for energies and magnetisation for each time-step will be saved in a file:
energiesT<T>_model_<model>.dat
And the energy plot can be found in:
plot_T<T>_energies_model_<model>.png

To run for multiple temperatures:
$ run(model = <model>, start = <start T>, end = <end T>, step = <step T>, size = <size>)
Where again, input the parameters like above, alongside the new parameters: <start T>, <end T>,
and <step T> corresponding to starting Temperature, Ending Temperature and Temperature step-size

After the simulation is run for multiple temperatures, the final observables can be generated
To calculate the final observables and plot the avg energy, avg magnetisation, susceptibility, 
and specific heat vs temperature, start a new python session in terminal by typing:
$ quit()
If the previous session is still running, and then typing 
$ python
Followed by
$ from genTPlots import *
Then run:
$ main(<model>)

The plots will be saved in files:
EvT_<model>.png, CvT_<model>.png, MvT_<model>.png, and SvT_<model>.png
with the data found in files:
observables_<model>.dat

Additional parameters and their definitions can be found in the in-code documentation