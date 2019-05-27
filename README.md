# GibbsIT
Gibbs energy from isotope tracing

Code used in:

Spare glycolytic capacity supports metabolic homeostasis and energy efficiency.

Junyoung O. Park, Lukas B. Tanner, Monica H. Wei, Daven Khana, Tyler Jacobson, Zheyun Zhang, Sara A. Rubin, Sophia Hsin-Jung Li, Meytal B. Higgins, David M. Stevenson, Daniel Amador-Noguez, and Joshua D. Rabinowitz

Princeton University, University of Wisconsin - Madison, UCLA

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