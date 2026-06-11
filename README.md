# Biochar Molecular Modelling Tutorials

By Valentina Erastova, University of Edinburgh, June 2026.

Contact: valentina.erastova@ed.ac.uk 

-----

These tutorials support the biochar modelling methodology introduced in the 2024 papers by Wood, Masek & Erastova ([10.1016/j.xcrp.2024.102036](https://doi.org/10.1016/j.xcrp.2024.102036), [10.1016/j.xcrp.2024.102037](https://doi.org/10.1016/j.xcrp.2024.102037)) and Ngambia, Masek & Erastova ([10.1016/j.biombioe.2024.107199](https://doi.org/10.1016/j.biombioe.2024.107199)). They show how to build biochar molecular models from experimental descriptors, and how to use those models in GROMACS adsorption studies, using a recent example from Wood, Masek & Erastova ([10.1016/j.biombioe.2026.109080](https://doi.org/10.1016/j.biombioe.2026.109080)).


## Assumed Background

These tutorial assume that you are comfortable with:

- Running GROMACS workflows from the command line
- Preparing and editing topology files
- Running energy minimisation, equilibration, and production simulations
- Submitting jobs to an HPC system
- Inspecting structures and trajectories in VMD
- Plotting and interpreting standard MD outputs

These tutorials build on that foundation to focus on the additional challenges specific to biochar: chemically heterogeneous surfaces, porous carbon-rich structures, experimental validation targets, and interfacial adsorption analysis. 

If you want to revise the above basics, MD Research Techniques workshop ([github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques](https://github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques)) is a good starting point. That workshop provides the GROMACS, HPC, topology, `.mdp`, visualisation, and analysis background assumed here. These tutorials do not repeat the full introductory setup.

------

## Tutorials

1. [`biochar_model_building_tutorial.md`](biochar_model_building_tutorial.md)

   Build a biochar molecular model from experimental targets using the iterative approach developed by Wood, Masek & Erastova and extended by Ngambia, Masek & Erastova. The tutorial covers target selection, building-block design, model assembly, condensation, validation against experimental properties, and preparation of a surface-exposed model.

2. [`biochar_solvated_gromacs_tutorial.md`](biochar_solvated_gromacs_tutorial.md)

   Set up and analyse adsorption simulations using published biochar surface models in GROMACS. The tutorial covers selecting a biochar model, solvating the surface, adding solutes and ions, running the simulation, checking convergence, and analysing adsorption at the biochar-water interface.


------

## Supporting Repositories

The tutorials refer to files, examples, and workflows from the following repositories:

- MD Research Techniques workshop: [github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques](https://github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques)
- Published biochar molecular models: [github.com/Erastova-group/Biochar_MolecularModels](https://github.com/Erastova-group/Biochar_MolecularModels)
- Example 2,4-D adsorption simulations: [github.com/Erastova-group/24D_biochar](https://github.com/Erastova-group/24D_biochar)
- Mn(II)-biochar simulation materials: [github.com/Erastova-group/Mn_Biochar](https://github.com/Erastova-group/Mn_Biochar)
- DynDen convergence analysis: [github.com/Erastova-group/DynDen](https://github.com/Erastova-group/DynDen)


---

## Citation

We hope these tutorials were useful for you! If you use our models or approach in your published work, **please cite our works** that enabled these tutorials.

Adsorption studies:

1. Wood R, Mašek O & Erastova V. *Unravelling 2,4-D – biochar interactions by molecular dynamics: adsorption modes and surface functionalities.* Biomas and Bioenergy 210, 109080 (2025). DOI: [10.1016/j.biombioe.2026.109080](https://doi.org/10.1016/j.biombioe.2026.109080)

2. Ngambia A, Gavrilova A, Huang H, Lyu Z, Mašek O, Graham M & Erastova V. *Decoupling Precipitation and Surface Complexation during Mn (II) Removal by Biochar via Experiments and Atomistic Simulations.* arXiv preprint (2026). DOI: [10.48550/arXiv.2603.22144](https://doi.org/10.48550/arXiv.2603.22144)

Biochar model development:

3. Wood R, Mašek O & Erastova V. *Developing a molecular-level understanding of biochar materials using public characterization data.* Cell Reports Physical Science 5(7), 102036 (2024). DOI: [10.1016/j.xcrp.2024.102036](https://doi.org/10.1016/j.xcrp.2024.102036)

4. Wood R, Mašek O & Erastova V. *Developing realistic molecular models of biochars.* Cell Reports Physical Science 5(7), 102037 (2024). DOI: [10.1016/j.xcrp.2024.102037](https://doi.org/10.1016/j.xcrp.2024.102037)

5. Ngambia A, Mašek O & Erastova V. *Development of biochar molecular models with controlled porosity.* Biomass and Bioenergy 184, 107199 (2024). DOI: [10.1016/j.biombioe.2024.107199](https://doi.org/10.1016/j.biombioe.2024.107199)

Convergence analysis:

6. Degiacomi MT, Tian S, Greenwell HC & Erastova V. *DynDen: Assessing convergence of molecular dynamics simulations of interfaces.* Computer Physics Communications 269 (2021). DOI: [10.1016/j.cpc.2021.108126](https://doi.org/10.1016/j.cpc.2021.108126)
