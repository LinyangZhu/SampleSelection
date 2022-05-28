# SampleSelection
The MATLAB version is R2016b

The data in 'SourceData.dat' is from the flowfield of NACA0012 airfoil with the freestream condition Ma = 0.15,Re= 6*10^6 with the angle of attack AOA=10. Please see https://turbmodels.larc.nasa.gov/naca0012_val.html.

Make sure that the MATLAB script and 'SourceData.dat' are in the same file path before running the code.

Run the 'Physics_assited_RecursiveMethod.m' code and the command-line window will output:
The number of original data: 57344
The number of sub-region is: 35 
Time consumption: 22.453125s 
The sample number afer sample selection is: 7998

Run the 'RecursiveMethod.m' code and the command-line window will output:
The number of original data: 57344
Time consumption: 426.890625s 
The sample number afer sample selection is: 7993

It should be note that 'Time consumption' can be different for different computers!

For more details, please see our manuscript 'Physics-assisted recursive method for sample selection from turbulence data'.

Contributor:Linyang Zhu
Email:zhulinyang@mail.nwpu.edu.cn
