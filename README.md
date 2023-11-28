# Introduction:
This project is about the paper *"Steady Mixing State of Black Carbon Aerosols from a Particle-Resolved Model"*. Here, you can find the nine scenario setups (parameter settings) for simulation in the 'PartMC-MOSAIC' model mentioned in the paper. However, this project does not include the results after running the model (you need to run 'PartMC-MOSAIC' with the provided scenarios on your own). Additionally, the project includes data preprocessing scripts, optical calculation scripts, and plotting scripts based on the 'PartMC-MOSAIC' simulation results. In summary, by reading this README.md file and the files in this project, you can reproduce the paper (with minor numerical deviations due to the random nature of the model simulations). It is important to note that 'PartMC-MOSAIC' in this project is run using a 'Docker' image, and you can choose to use the Docker image or install PartMC-MOSAIC on your own server according to your needs. Installation instructions for PartMC can be found at http://lagrange.mechse.illinois.edu/partmc/, and the code for MOSAIC can be obtained from Rahul A. Zaveri.

***
## 1. Project Structure
This section presents the project structure as an unordered list and provides a brief explanation for each part. For detailed explanations of each file, please refer to the corresponding sections. For example, if you want to understand the description and usage of the ```Figure1_a.py``` script, you can refer to section 4 **The plot files of this project**.

- **Project**: *the whole project*
    - **sce11**
    - **sce12**
    - **...**    
    - **sce22**: *the basic scenario*
        - **.dat**: *PartMC input file format*
        - **.spec**: *PartMC setup File*
        - **python**: *directory of Python scripts*
        - **Figure**: *directory of figures*
        - **Data**: *directory of pre-process data*
        - **out**: *directory of PartMC output*
        - **matlab**: *directory of matlab scripts*
    - **...**
    - **sce33**

***
## 2. How to run PartMC-MOSAIC
In this section, you can find instructions on how to install and run 'PartMC-MOSAIC'. Additionally, using Scenario 22 as an example, the ```.dat``` and ```.spec``` files required for running 'PartMC-MOSAIC' are explained.

### 2.1 Installation of PartMC and MOSAIC
You can refer to http://lagrange.mechse.illinois.edu/partmc/ for the installation of 'PartMC'. The link provides various versions of the 'PartMC' program, installation tutorials, and documentation. The code for 'MOSAIC' can be obtained from Rahul A. Zaveri.

### 2.2 Introduction to .spec and .dat files
The ```urban_plume.spec``` file specifies the running mode, total particle number, simulation duration, meteorological parameters, background field, emission sources, and atmospheric reaction processes for PartMC simulations. (For detailed information and explanations, please refer to http://lagrange.mechse.illinois.edu/partmc/partmc-2.6.1/doc/html/index.html )

- **Temperature data throughout the entire simulation process (input):** 
```temp.dat```

- **Pressure data throughout the entire simulation process (input):** 
```pres.dat```

- **Boundary layer height data throughout the entire simulation process (input):** 
```height.dat``` 

- **Emission and initial fields of aerosol species (input):** 
```aero_data.dat, aero_back.dat, aero_back_comp.dat, aero_back_dist.dat, aero_init_comp.dat, aero_init_dist.dat, aero_emit.dat, aero_emit_dist.dat, aero_emit_comp_BC.dat, aero_emit_comp_BCfree.dat, aero_emit_comp_mixBC.dat```:

- **Emission and initial fields of gas species (input):** 
```gas_data.dat, gas_init.dat, gas_back.dat, gas_emit.dat```:


### 2.3 Running PartMC-MOSAIC
In this project, if you are using the Docker image to run 'PartMC-MOSAIC', you only need to enter the command ```./run.sh``` to run it. If you have installed 'PartMC-MOSAIC' locally, you need to manually modify the PartMC-MOSAIC running command in the ```./run.sh``` file.
***
## 3. The pre-process data of this project
In this section, you can learn about the preprocessing data, how to obtain it, and how to use it. Before running the Python script, make sure that the required libraries for the script are installed, otherwise errors may occur.

### 3.1 Introduction of pre-process scripts
The preprocessing scripts are used to process the output data of 'PartMC-MOSAIC' in nc format to obtain the required information for the study. In this project, the main preprocessing script is the set of functions in ```NC_X.py```, which can output information such as individual particle Dc and Dp, mass of each component, and number concentration of black carbon particles. Since the individual particle Dc and Dp information is repeatedly used in this project and both 'Python' and 'MATLAB' are used for data processing and plotting, some data information is stored in CSV files to improve overall efficiency.

### 3.2 How to get the pre-process data
The ```Dp_Dc_Conc.csv``` file is obtained by running the Python script ```saveDcDpconc.py```.
The ```MAC_MACE_Rabs_MACEX_ABC_time.csv``` file is generated from the ```Dp_Dc_Conc.csv``` file using the 'MATLAB' script ```Eabs_perparticle.m```, which processes the optical parameters of black carbon aerosols at each time point.
The ```VF_massBC.csv``` file is the result of running the 'Python' script ```Figure4_a.py```, representing the volume fraction of black carbon components in black carbon aerosols for different core particle sizes.

### 3.3 How to use the pre-process data
The ```Dp_Dc_Conc.csv``` file is used for studying the distribution of coating thickness for black carbon particles (Figure2 and Figure4_a) and calculating optical properties (Figure3).
The ```MAC_MACE_Rabs_MACEX_ABC_time.csv``` file is used to compare the differences in evaluating the optical properties of black carbon aerosols using per-particle calculation and k-value methods (Figure3 and Table S7).
The ```VF_massBC.csv``` file represents the volume fraction of black carbon components in black carbon aerosols for different core particle sizes (Figure 4_b).

***
## 4. The plot and output files of this project
You can learn about the instructions and usage of drawing and screen output scripts used in this project in this section. Please note that the drawing scripts in the python folder all use relative paths, so the python command needs to be run in their respective directories. Before running the Python script, make sure that the required libraries for the script are installed, otherwise errors may occur.

```NC_X.py```&```move.py```: Models needed to import 

```Figure1_a.py```: The plot script of Figure 1a and Table S3 in the manuscript.

```Figure1_b.py```: The plot script of Figure 1b in the manuscript.

```Figure2_ab.py```: The plot script of Figure 2 in the manuscript, need the file ```'../Data/Dp_Dc_Conc.csv'```.

```Figure3.py```: The plot script of Figure 3 in the manuscript; the output script of Table S7 in the SI, need the file ```'../Data/MAC_MACE_Rabs_MACEX_ABC_time.csv'```.

```Figure4_a.py```: The plot script of Figure 4a in the manuscript.

```Figure4_d.py```: The plot script of Figure 4d in the manuscript.

```Valiadation.py```: The output script of Table S2 in the SI.

```Eabs_Kvalue.m```: The output script of the value of MAC and E_{abs} in the Figure 3 and Table S7, need to use the value of k calculated by ```Figure2_ab.py```.

``` Figure4_bc.m```: The plot script of Figure 4b and Figure 4c in the manuscript. This file need to input the value k according to the Figure 2a.

```Figure_S1.py```: The plot script of Figure S1 in the manuscript.

```print_50to90nm.py```: Print the number concentration fraction of BC particle (range from 50 to 90 nm).

```print_meanCT.py```: Print the meanCT calculated by per-particle method, used to make a comparision with 1/k.




