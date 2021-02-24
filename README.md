# A systems biology analysis of lipolysis and fatty acid release from adipocytes *in vitro* and from adipose tissue *in vivo*

This repository corresponds to the results and analysis used in the paper: "A systems biology analysis of lipolysis and fatty acid release from adipocytes *in vitro* and from adipose tissue *in vivo*". 
A [preprint version of the paper is available at bioRxiv.](https://www.biorxiv.org/content/10.1101/2020.12.18.423229v1)

### Abstract
Lipolysis and the release of fatty acids to supply energy to other organs, such as between meals, during exercise, and starvation, are fundamental functions of the adipose tissue. The intracellular lipolytic pathway in adipocytes is activated by adrenaline and noradrenaline, and inhibited by insulin. Circulating fatty acids are elevated in type 2 diabetic individuals. The mechanisms behind this elevation are not fully known, and to increase the knowledge a link between the systemic circulation and intracellular lipolysis is key. However, data on lipolysis and knowledge from *in vitro* systems have not been linked to corresponding *in vivo* data and knowledge *in vivo*. Here, we use mathematical modelling to provide such a link. We examine mechanisms of insulin action by combining *in vivo* and *in vitro* data into an integrated mathematical model that can explain all data. Furthermore, the model can describe independent data not used for training the model. We show the usefulness of the model by simulating new and more challenging experimental setups *in silico*, e.g. the extracellular concentration of fatty acids during an insulin clamp, and the difference in such simulations between individuals with and without type 2 diabetes. Our work provides a new platform for model-based analysis of adipose tissue lipolysis, under both non-diabetic and type 2 diabetic conditions. 

## Required software
These scripts are implemented in MATLAB (tested on R2020a), and requires a MATLAB compatible C-compiler. The easiest way to check if a valid compiler is available is to run `mex -setup` in MATLAB and check if it finds a compiler. For more information about compatible compilers, refer to [MathWorks support page](https://se.mathworks.com/support/requirements/supported-compilers.html).

## How to recreate the plots from the paper:
1. Download the files (as .zip or through git)
2. If downloaded as .zip, extract the files
3. Open MATLAB (tested on R2020a)
4. Change folder in MATLAB to the folder with the files
5. To plot specific figures:
    * Run either of the script `PlotAllFigures` or `PlotAllFigures_fast` to plot all figures. Note that the fast version runs with a lower resolution than the one used in the paper to save time.
    * Run the script `PlotSpecificFigures`, and follow the instructions in the MATLAB command window to select specific plots to recreate. 
    

Note that generating high resolution plots will take some time. A low resolution faster option is available to save time.

