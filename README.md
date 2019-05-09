Welcome to the DeNoGAP_v2 Tutorial

## **1. INSTALLATION**

**Perquisites required**

- Installation of DeNoGAP2 pipeline and programs used within it requires Miniconda3 & Python3.7 or above to be installed on the system.

- Download & install Minicona3 for Python3.7 or above on your system from here [Miniconda](https://conda.io/miniconda.html)

- Miniconda installer for Mac users: [Miniconda (MacOS)](https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh)

- Miniconda installer for Linux users: [Miniconda (Linux)](https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh)

- Run following commands to install miniconda3:

`mkdir miniconda3`<br>
`cd miniconda3`<br>
`wget https://repo.continuum.io/miniconda/Miniconda3-latest-<OS>.sh`<br>
`sh Miniconda3-latest-<OS>.sh`

**Instructions to setup & install DeNoGAP2**

- Create a conda environment for DeNOGAP2

`conda create --name denogap2`<br>

- Add channels for installing required programs & packages

`conda config --add channels defaults`<br>
`conda config --add channels conda-forge`<br>
`conda config --add channels bioconda`<br>

**Install required programs and supporting packages in DeNoGAP2 conda environment**

`conda install --name denogap2 python=3.7`<br>
`conda install --name denogap2 biopython hmmer blast mcl scipy sqlite mafft psutil pandas`<br>

**2. How to run DeNoGAP pipeline**

- Activate DeNoGAP2 Conda environment:

`source activate denogap2`<br>

- Run Commands:

To see help for running DeNoGAP type following command:

`python denogap_hmm.py -h`
