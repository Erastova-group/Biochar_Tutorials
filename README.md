# Biochar Molecular Modelling Tutorials

These tutorials support the biochar modelling methodology introduced in the 2024 papers by Wood, Masek & Erastova and Ngambia, Masek & Erastova. They show how to build biochar molecular models from experimental descriptors, and how to use those models in GROMACS adsorption studies.

They are intended as follow-on material from the MD Research Techniques workshop:

[github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques](https://github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques)

That workshop provides the GROMACS, HPC, topology, `.mdp`, visualisation, and analysis background assumed here. These tutorials do not repeat the full introductory setup.

## Tutorials

1. [`biochar_model_building_tutorial.md`](biochar_model_building_tutorial.md)

   Build a biochar molecular model from experimental targets using the iterative approach developed by Wood, Masek & Erastova and extended by Ngambia, Masek & Erastova. The tutorial covers target selection, building-block design, model assembly, condensation, validation against experimental properties, and preparation of a surface-exposed model.

2. [`biochar_solvated_gromacs_tutorial.md`](biochar_solvated_gromacs_tutorial.md)

   Set up and analyse adsorption simulations using published biochar surface models in GROMACS. The tutorial covers selecting a biochar model, solvating the surface, adding solutes and ions, running the simulation, checking convergence, and analysing adsorption at the biochar-water interface.

## Supporting Repositories

The tutorials refer to files, examples, and workflows from the following repositories:

- MD Research Techniques workshop: [github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques](https://github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques)
- Published biochar molecular models: [github.com/Erastova-group/Biochar_MolecularModels](https://github.com/Erastova-group/Biochar_MolecularModels)
- Example 2,4-D adsorption simulations: [github.com/Erastova-group/24D_biochar](https://github.com/Erastova-group/24D_biochar)
- Mn(II)-biochar simulation materials: [github.com/Erastova-group/Mn_Biochar](https://github.com/Erastova-group/Mn_Biochar)
- DynDen convergence analysis: [github.com/Erastova-group/DynDen](https://github.com/Erastova-group/DynDen)

## Assumed Background

Before starting, you should be comfortable with:

- Running GROMACS workflows from the command line
- Preparing and editing topology files
- Running energy minimisation, equilibration, and production simulations
- Submitting jobs to an HPC system
- Inspecting structures and trajectories in VMD
- Plotting and interpreting standard MD outputs

These tutorials build on that foundation to focus on the additional challenges specific to biochar: chemically heterogeneous surfaces, porous carbon-rich structures, experimental validation targets, and interfacial adsorption analysis.

## Citation

If you use these tutorials, models, or workflows in published work, please cite the relevant papers listed in the tutorial reference sections.
