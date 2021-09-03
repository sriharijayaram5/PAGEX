**Proton, Alpha, Gamma, Electron and X-radiation interaction parameters**

Description
===========

A complete package for computation of various radiation interaction
parameters useful in diverse fields. The radiations considered include
X-/gamma rays, electrons, protons and alpha particles.

Authors
-------

| **Original release** : October 2019
| Srilakshmi Prabhu :sup:`1` *srilakshmi.prabhu@res.christuniversity.in*
| Sreehari Jayaram :sup:`2` *s.jayaram@pi3.uni-stuttgart.de*
| Bubbly S G :sup:`1` *bubbly.sg@christuniversity.in*
| S B Gudennavar :sup:`1` *shivappa.b.gudennavar@christuniversity.in*

| :sup:`1` Department of Physics and Electronics, CHRIST (Deemed to be
  University), Bengaluru-560029, Karnataka, India
| :sup:`2` 3. Physikalisches Institut, Universität Stuttgart, Stuttgart-70569, Stuttgart, Germany


Various parameters and their relevant atomic number and energy range 
computed by PAGEX are as follows:

Photon interaction parameters
-----------------------------

.. raw:: html

   <table>

.. raw:: html

   <tbody>

.. raw:: html

   <tr>

   <th>Description</th>

   <th>Parameters</th>

   <th>Atomic number range</th>

   <th>Energy range</th>

   </tr>

.. raw:: html

   <tr>

   <td>Partial/total cross section per atom (b/atom) (with and without coherent scattering)</td>

   <td>σ<sub>coh</sub>, σ<sub>incoh</sub>, σ<sub>e</sub>, σ<sub>pair</sub>, σ<sub>trip</sub>, σ<sub>w/coh</sub>, σ<sub>wo/coh</sub> (for elements)</td>

   <td>1-99</td>

   <td>1 keV - 100 GeV</td>

   </tr>

.. raw:: html

   <tr>

   <td>Mass attenuation coefficients (cm<sup>2</sup>/g) (with and without coherent scattering)</td>

   <td>μ/ρ<sub>coh</sub>, μ/ρ<sub>incoh</sub>, μ/ρ<sub>e</sub>, μ/ρ<sub>pair</sub>, μ/ρ<sub>trip</sub>, μ/ρ<sub>w/coh</sub>, μ/ρ<sub>wo/coh</sub> (for elements/compounds)</td>

   <td>1-99</td>

   <td>1 keV - 100 GeV</td>

   </tr>

.. raw:: html

   <tr>

   <td>Average cross section per atom (b/atom) and average cross section per electron (b/electron)<sup>a</sup></td>

   <td>σ<sub>a</sub> and σ<sub>e</sub></td>

   <td>1-99</td>

   <td>1 keV - 100 GeV</td>

   </tr>

.. raw:: html

   <tr>

   <td>Effective atomic number (Z<sub>eff</sub>) and electron density (N<sub>eff</sub>, electrons/g)</td>

   <td>Z<sub>eff</sub> and N<sub>eff</sub></td>

   <td>1-99</td>

   <td>1 keV - 100 GeV</td>

   </tr>

.. raw:: html

   <tr>

.. raw:: html

   <td>

.. raw:: html

   <tr>

   <td>Mass-energy absorption coefficients (μ<sub>en</sub>/ρ, cm<sup>2</sup>/g) and relative KERMA<sup>b</sup></td>

   <td>μ<sub>en/ρ</sub> and K<sub>R</sub></td>

   <td>1-92</td>

   <td>1 keV - 20 MeV</td>

   </tr>

.. raw:: html

   <tr>

   <td>Photon energy absorption effective atomic number and electron density (electrons/g)</td>

   <td>Z<sub>PEAeff</sub>, N<sub>PEAeff</sub> (electrons/g)</td>

   <td>1-92</td>

   <td>1 keV - 20 MeV</td>

   </tr>

.. raw:: html

   <tr>

   <td>Equivalent atomic number</td>

   <td>Z<sub>eq</sub></td>

   <td>1-99</td>

   <td>1 keV - 100 GeV</td>

   </tr>

.. raw:: html

   <tr>

   <td>G-P fitting parameters and buildup factors</td>

   <td>a, b, c, X<sub>k</sub> and d, EABF and EBF</td>

   <td>Z<sub>eq</sub> ~ 4-82</td>

   <td>15 keV - 15 MeV</td>

   </tr>

.. raw:: html

   </tbody>

.. raw:: html

   </table>

a - Note: 1 barn = 10-24cm2 b - Relative to air and 11 human tissues
(adipose tissue, blood, bone, brain, breast, eye, muscles, ovary,
testis, skin and soft tissue) as reported in ICRU Report 44 (1989)

Charged particle interaction parameters
---------------------------------------

.. raw:: html

   <table id="chtb">

   <tbody>

   <tr>

   <th>Description</th>

   <th>Parameters</th>

   <th>Atomic number range</th>

   <th>Energy range</th>

   </tr>

   <tr>

   <td>Mass stopping power (MeV cm<sup>2</sup>/g)</td>

   <td>S(E)/ρ</td>

   <td rowspan="3">1-98 (electrons) 1-92 (alpha particles) 1-92 (protons)</td>

   <td rowspan="3">10 keV - 1 GeV (electrons) 1 keV - 1 GeV (alpha particles) 1 keV - 10 GeV (protons)</td>

   </tr>

   <tr>

   <td>Mass stopping cross section (MeV cm<sup>2</sup>/atom)</td>

   <td>S<sub>c</sub></td>

   </tr>

   <tr>

   <td>Effective atomic number and electron density (electrons/g)</td>

   <td>Z<sub>eff</sub> and N<sub>eff</sub></td>

   </tr>

   </tbody>

   </table>
|
Data files used in PAGEX
-------------------------

.. raw:: html

   <table>

   <tbody>

   <tr>

   <th>File</th>

   <th>Data</th>

   </tr>

   <tr>

   <td>NIST/MDATX3n.xxx</td>

   <td>Atomic cross section data for each element</td>

   </tr>

   <tr>

   <td>AStar_data/DATAxxx</td>

   <td>Stopping-power data of helium ions in each element</td>

   </tr>

   <tr>

   <td>PStar_data/DATAxxx</td>

   <td>Stopping-power data of protons in each element</td>

   </tr>

   <tr>

   <td>EStar_data/DATAxxx</td>

   <td>Stopping-power data of electrons in each element</td>

   </tr>

   <tr>

   <td>XRay_data1/DATAnxx</td>

   <td>Element photon mass attenuation coefficient (μ/ρ) and mass energy-absorption coefficient (μ<sub>en</sub>/ρ) data</td>

   </tr>

   <tr>

   <td>XRay_Comp1/DATAn_material</td>

   <td>Material photon mass attenuation coefficient (μ/ρ) and mass energy-absorption coefficient (μ<sub>en</sub>/ρ) data</td>

   </tr>

   <tr>

   <td>ANSI_data/ANSI_(A/B)/DATAn_x</td>

   <td>G-P fitting parameters a, b, c, d, X<sub>k</sub></td>

   </tr>

   </tbody>

   </table>

where 'xxx' and 'xx' represent the atomic number.

PAGEX I/O
---------

User Input
..........

The user can input the composition of material, either in terms of molecular formula or as fractional (by mass) elemental composition, 
which is then stored in temporary memory. All the data regarding the basic properties of individual elements is extracted. 
The formula input is case sensitive. The composition data entered will be used to load the corresponding standard data files stored in the program directory. 
The user can choose either only standard energy range or the standard energy range with user defined energies for computation of desired parameters. 
However, the chosen energy should be within the available range.
An input log ('InputLog.log') is maintained that logs all input made to the program
with a time-stamp when PAGEX GUI is used.

PAGEX output
............

An attractive feature of PAGEX is that it allows the user to save the output data in an easily extractable CSV file or plot the desired parameter as 
functions of energy. The data in CSV file contains all the quantities related to the calculation in addition to the specific output parameter. 
For instance, mass
attenuation coefficients are written in the file named "Photon mass
attenuation and interaction cross section parameters” in the folder "Save_Folder" which is automatically created if not present.
Graphs are also plotted using the GUI as well as Matplotlib (Hunter, 2007) when using the API from comman line. These graphs can be easily saved in various publication ready formats.

Installation
============

-  Clone the github repository
   https://github.com/sriharijayaram5/PAGEX.git on your local machine.

-  Install python3.x from https://www.python.org/downloads/.

-  Run the following commands in your command prompt in the directory
   you wish to install a virtual python environment.

-  ::

        pip install virtualenv
        virtualenv my_env_name
        my_env_name\Scripts\activate.bat

-  Navigate to the PAGEX directory.

-  ::

      pip install -r requirements.txt

-  | You can now run the PAGEX program GUI by the following command:

-  ::

      python pagex_worker.py

-  The program can also be run from terminal by importing the
   ``Compound`` class:

-  ::

         from pagex_worker import Compound
         my_comp = Compound(comp='C 6 H 12 O 6')
         my_comp.myu()
         my_comp.plot_parameter()
         my_comp.write_to_csv()
         #To run gui
         import pagex_worker
         pagex_worker.run_gui()

-  An iPython REPL is recommended for terminal use. See docstring for
   help and arguments. Users are encouraged to read :ref:`api-desc` for more details.

-  ::

         help(my_comp.myu)
         
Disclaimer
===========

The developers of PAGEX have put in their sincere efforts to deliver a
program that computes accurate data and relies on other standard
databases. However, the developers make no warranties to that effect,
and the developers or affiliates shall not be liable for any damage that
may result from errors or omissions in the calculations or database
used.

License
========

Copyright 2019 Bubbly S G Permission is hereby granted, free of charge,
to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to
do so, subject to the following conditions: The above copyright notice
and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.