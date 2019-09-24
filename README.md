# GibbsIT
Gibbs energy from isotope tracing

Code used in:

Near-equilibrium glycolysis supports metabolic homeostasis and energy yield.

Junyoung O. Park, Lukas B. Tanner, Monica H. Wei, Daven B. Khana, Tyler B. Jacobson, Zheyun Zhang, Sara A. Rubin, Sophia Hsin-Jung Li, Meytal B. Higgins, David M. Stevenson, Daniel Amador-Noguez, and Joshua D. Rabinowitz

University of California, Los Angeles, Princeton University, University of Wisconsin - Madison

Nature Chemical Biology (2019)

https://doi.org/10.1038/s41589-019-0364-9

In each organism folder,
- *_mea.xlsx contains measured labeling fractions and net/exchange flux constraints that went into the model as input
- .xml shows the metabolites, reactions, and carbon and hydrogen mapping
- .m functions contain the stoichiometric matrix, its kernel, and the EMU (elementary metabolite unite) model
- .mat file contains stoichiometric matrix, its kernel, and measured metabolite labeling

The following files contain confidence intervals of dG
- nupshift.mat, pupshift.mat, oupshift.mat, oligomycin.mat, and clostridia.mat; these can be used with ./src/area_overlap_script.m script to find the overlapping area between two dG values

To run MFA, obtain Gibbs energy of reaction, and integrate metabolite concentrations and Gibbs free energies,
1) Open Matlab 2013b or newer
2) Set Path -> Add Folder -> Choose ./src -> Save
3) Open (*_)script.m and execute
4) To run it parallel, enter 'matlabpool local 4' on the command line. '4' can be a different number if a different number of cores is desired
5) Running on clusters requires an additional setup. Please contact system administrator