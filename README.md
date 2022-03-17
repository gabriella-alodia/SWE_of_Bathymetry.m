# SWE_of_Bathymetry.m

(c) Gabriella Alodia &#8212; 2021

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

<b>Preprint Citation:</b>

Alodia, G., Green, C., McCaig, A., & Paton, D. (2020). Slope-Weighted Eccentricity: Automatic Terrain Classification of Atlantic Ocean Crust. ESSOAr: Earth and Space Science Open Archive. [https://doi.org/10.1002/essoar.10502634.1]

<b>EGU Presentation:</b>

Alodia, G., Green, C., McCaig, A., & Paton, D. (2020). A novel approach for oceanic spreading terrain classification at the Mid-Atlantic Ridge using Eigenvalues of high-resolution bathymetry. In EGU General Assembly Conference Abstracts (p. 337). [https://presentations.copernicus.org/EGU2020/EGU2020-337_presentation.pdf]

Full paper currently under review.
