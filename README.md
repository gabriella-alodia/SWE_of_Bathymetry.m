# SWE_of_Bathymetry.m

[![View SWE_of_Bathymetry on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/108139-swe_of_bathymetry)

(c) Gabriella Alodia, Chris M. Green, Andrew M. McCaig &#8212; 2022

<tt>SWE_of_Bathymetry.m</tt> is a Matlab-based geomorphometric algorithm to obtain the numerical description of both magmatic and tectonic crust in a slow-spreading ridge through a series of calculation based on the distribution of the azimuth and plunge observed in the seafloor morphology.

<b>Input:</b>
A gridded shipborne multibeam bathymetry (depths in metres) in *.xyz format (here: <tt>'Input_Bathymetry_15s.xyz'</tt>)

![Input](https://github.com/gabriella-alodia/SWE_of_Bathymetry/blob/main/Figure_Input.png)

<b>Output</b>
1. Terrain eccentricity (here: <tt>'Output_eccentricity.xyz'</tt>)
2. Weight matrix: 1-sin(slope) (here: <tt>'Output_weight.xyz'</tt>)
3. SWE: Slope-weighted eccentricity (here: <tt>'Output_SWE.xyz'</tt>)
4. Masked SWE (here: <tt>'Output_SWE_masked.xyz'</tt>)

![Output](https://github.com/gabriella-alodia/SWE_of_Bathymetry/blob/main/Figure_Output.png)

Each output is exported in *.xyz format. The resulting *.xyz data can be converted into *.grd using the <tt>xyz2grd</tt> function in GMT (http://gmt.soest.hawaii.edu/doc/5.3.2/xyz2grd.html)

The shipborne multibeam bathymetry data sample is downloaded from the GMRT MapTool (https://www.gmrt.org/GMRTMapTool/) with the extent <tt>xmin/xmax/ymin/ymax</tt> of <tt>-46/-44/12.5/13.15</tt>

<b>Full paper:</b>

Alodia, G., Green, C. M., & McCaig, A. M. (2022). SWE_of_Bathymetry.m: A geomorphometric tool to automate discrimination between detachment and magmatic seafloor at slow-spreading ridges from shipboard multibeam bathymetry, Computers & Geosciences, Volume 166, 105177. [https://doi.org/10.1016/j.cageo.2022.105177]

<b>Preprint & EGU Presentation:</b>

Alodia, G., Green, C., McCaig, A., & Paton, D. (2020). A novel approach for oceanic spreading terrain classification at the Mid-Atlantic Ridge using Eigenvalues of high-resolution bathymetry. In EGU General Assembly Conference Abstracts (p. 337). [https://presentations.copernicus.org/EGU2020/EGU2020-337_presentation.pdf]

Alodia, G., Green, C., McCaig, A., & Paton, D. (2020). Slope-Weighted Eccentricity: Automatic Terrain Classification of Atlantic Ocean Crust. ESSOAr: Earth and Space Science Open Archive. [https://doi.org/10.1002/essoar.10502634.1]
