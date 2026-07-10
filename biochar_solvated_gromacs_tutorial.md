# Setting Up Biochar Adsorption Simulations in GROMACS

*Tutorial written by Valentina Erastova (valentina.erastova@ed.ac.uk), University of Edinburgh, in June 2026*

*Based on Wood, Mašek & Erastova (BiomBioe, 2026) and Ngambia, et al (ArXiv, 2026)*

---

## Who is this tutorial for?

You have completed the MD Research Techniques course ([github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques](https://github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques)) — you know how to set up a protein and clay simulations in GROMACS, run energy minimisation and an equilibration on HPC, and make sense of what comes out. Biochar is the next step. The material is more complex, the models are larger, and the preparation is less automated than what you are used to. This tutorial walks you through the entire process of taking the Erastova group's published biochar molecular models and getting them running in GROMACS, from downloading the files to checking that your system has equilibrated. If you want to build your own biochar model, we recommend you follow our "Building Biochar Molecular Models from Scratch" tutorial instead.

By the end, you will have a working adsorption simulation: a biochar slab submerged in a solution of water, a pollutant molecule, and charge-balancing counterions — the same setup used in the group's published studies of 2,4-D herbicide adsorption on biochar (Wood, Mašek & Erastova, 2025) and Mn(II) removal by biochar (ArXiv, in review).

---

## Background: What is biochar and why does it need a special model?

Biochar is a porous, carbon-rich solid produced by heating biomass (wood, straw, agricultural waste) in the absence of oxygen — a process called *pyrolysis*. The temperature at which this is done, referred to as the *highest treatment temperature* (HTT), has an enormous effect on the resulting material. Low-temperature biochars (around 400 °C) are chemically rich, with many oxygen-containing surface groups such as carboxyls, phenols, carbonyls, and ethers. High-temperature biochars (800 °C and above) are predominantly aromatic — large fused graphitic sheets — with few surface functional groups. Medium-temperature materials sit between these extremes.

This HTT-dependence matters deeply for applications. Biochar is actively studied as an inexpensive, environmentally friendly sorbent for removing pollutants from water and soil. The functional groups at its surface govern what it adsorbs and how strongly, so understanding those interactions at the molecular level is directly useful for designing better materials.

The problem, until recently, was the absence of suitable molecular models. Biochars were frequently approximated as pure graphene sheets in simulations — a severe oversimplification that misses the chemistry that makes real biochars interesting. We (Wood, Mašek & Erastova) addressed this by developing the first realistic molecular models of woody biochars, validated against experimental H/C and O/C ratios, aromaticity indices, true densities, and electron microscopy morphologies. These models, representing biochars produced at 400 °C (BC400), 600 °C (BC600), and 800 °C (BC800), are freely available and ready for immediate use in GROMACS. A further set of microporous models representative of 600–650 °C biochars was developed in Erastova group (Ngambia, Mašek & Erastova), extending coverage to materials with controlled pore size distributions.

---

## What you will need

These simulations are heavy, so access to a High-Performance computing resource is a must - you will not be able to do these simulations of your laptop.

**Software:**
- GROMACS 2022 or newer (the published simulations used GROMACS 2022 onwards)
- Python 3 with MDAnalysis, [pytim](https://marcello-sega.github.io/pytim/), and Matplotlib (for analysis, later)
- VMD (for visualisation, optional but strongly recommended)
- A text editor

**Files to download** (all freely available on our GitHub):

- Biochar models: [github.com/Erastova-group/Biochar_MolecularModels](https://github.com/Erastova-group/Biochar_MolecularModels)
- Example adsorption simulation files (2,4-D on biochar): [github.com/Erastova-group/24D_biochar](https://github.com/Erastova-group/24D_biochar)
- For Mn(II) simulations: [github.com/Erastova-group/Mn_Biochar](https://github.com/Erastova-group/Mn_Biochar)

---

## Step 1: Understand what is in the repository

Clone or download the `Biochar_MolecularModels` repository. The structure is:

```
Biochar_MolecularModels/
├── building_blocks/      # Individual molecular structures (.mol, .pdb) + chemical descriptors
├── exper_data/           # Experimental datasets used to validate the models
└── models/               # Ready-to-run GROMACS files
    ├── BC400/
    │   ├── BC400_bulk.gro         # Bulk (periodic) model
    │   ├── BC400_surf.gro         # Surface-exposed model
    │   ├── BC400.top              # Topology file
    │   └── oplsaa.ff/             # OPLS-AA force field directory
    ├── BC600/
    └── BC800/
```

The key distinction is between the *bulk* model and the *surface* model. The bulk model is periodic in all three directions (x, y, and z) and represents the interior of the biochar material. The surface model has been prepared by expanding the simulation box in the z-direction, leaving empty space above and below the biochar slab — this is the form you need for adsorption studies, because it creates an accessible surface where molecules can bind.

A word on the force field: all models use OPLS-AA (Optimised Potentials for Liquid Simulations, All-Atom). The `oplsaa.ff` directory must remain in the same folder as your topology file when you run any GROMACS command — GROMACS will look for it there. Do not move the `.top` file without bringing `oplsaa.ff` with it.

---

## Step 2: Inspect the surface model before you start

Before adding any solvent or solute, always look at what you have. Open the surface `.gro` file in VMD:

```bash
vmd BC400_surf.gro

```

You should see a slab of biochar material occupying the middle of the simulation box, with empty space above and below. The simulation box dimensions are written in the last line of the `.gro` file (in nm). For BC400, this is approximately 6.5 × 6.5 × 15.2 nm. The z-dimension is significantly larger than x and y precisely to accommodate the solution layer on each side.

If the structure looks fragmented or does not display properly, check that you have the complete `.gro` file — sometimes partially downloaded files cause problems.

---

## Step 3: Add solvent and solute — the system you want to simulate

For an adsorption study, you want the biochar surface submerged in a solution. The general approach is:

1. Solvate the empty space in the simulation box with water
2. Add solute molecules (the pollutant) at your desired starting concentration
3. Add counterions to charge-balance the system

The published simulations (Wood, Mašek & Erastova, 2025) used SPC water and OPLS-AA parameters throughout — this is important for consistency with the force field already assigned to the biochar.

### 3a. Solvating with GROMACS

```bash
gmx solvate -cp BC400_surf.gro -cs spc216.gro -o BC400_solvated.gro -p BC400.top

```

Here `-cp` specifies the configuration (your biochar surface), `-cs` specifies the solvent configuration to use (`spc216.gro` is the standard pre-equilibrated SPC water box that ships with GROMACS), `-o` is the output, and `-p` updates the topology file with the number of water molecules added.

Check how many water molecules were added — GROMACS will print this to the terminal. Open `BC400_solvated.gro` in VMD to confirm the box looks sensible: the biochar slab should be surrounded by water on both faces.

### 3b. Preparing your solute

If you are working with the 2,4-D example from the `24D_biochar` repository, the solute structures and topologies are already provided. If you are introducing your own molecule, you need to:

1. Draw (in MarvinSketch as an example) or download somewhere the molecule structure in a supported format (`.mol`, `.pdb`, `.sdf`)
2. Generate OPLS-AA parameters using [LigParGen](http://zarbi.chem.yale.edu/ligpargen) or, for polymeric structures, [PolyParGen](http://polypargen.com)
3. Download the resulting `.gro`/`.itp` files

LigParGen accepts a SMILES string or uploaded structure file and returns a GROMACS-compatible topology (`.itp`) and coordinate file (`.gro`). For example, 2,4-D in the SMILES is `OC(=O)COc1ccc(Cl)cc1Cl` — in its deprotonated anionic form at environmental pH, you need the carboxylate: `[O-]C(=O)COc1ccc(Cl)cc1Cl`.

Once you have the solute `.gro` and `.itp` files, include the `.itp` in your topology by adding this line near the top of `BC400.top`, after the force field include but before the `[ system ]` directive:

```
#include "24D.itp"

```

### 3c. Inserting solute molecules

Use `gmx insert-molecules` to place your solute at the desired number of copies:

```bash
gmx insert-molecules -f BC400_solvated.gro -ci 24D.gro -nmol 100 -o BC400_with_24D.gro

```

GROMACS will attempt to place 100 copies of 2,4-D randomly in the available solvent space, avoiding clashes with the biochar and existing water. It will report how many it managed to insert — if it could not place all of them due to space constraints, try a smaller number or a lower initial concentration.

Update the `[ molecules ]` section at the bottom of your topology file to reflect the number of solute molecules inserted.

### 3d. Adding counterions

Since 2,4-D is anionic (charge −1 per molecule), 100 molecules require 100 Na⁺ to make the system neutral. GROMACS can add these automatically, replacing water molecules:

```bash

# First, do a compilation step to prepare, use em.mdp file
gmx grompp -f em.mdp -c BC400_with_24D.gro -p BC400.top -o ions.tpr

# Neutralise the system
gmx genion -s ions.tpr -o BC400_neutralised.gro -p BC400.top -pname NA -nname CL -neutral
```

When prompted, choose group `SOL` (water) to replace with ions. The system is now charge-neutral.

---

## Step 4: Energy minimisation

Before any dynamics, you must resolve any clashes introduced by the packing of molecules. The energy minimisation settings used in the published simulations are:

```ini
; em.mdp — Energy minimisation parameters
integrator      = steep         ; Steepest descent algorithm
emtol           = 500.0         ; Stop when max force < 500 kJ/mol/nm
emstep          = 0.01          ; Initial step size (nm)
nsteps          = 50000         ; Maximum number of steps

; Non-bonded settings
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.4           ; Coulomb cut-off (nm)
rvdw            = 1.4           ; vdW cut-off (nm)
coulombtype     = PME           ; Particle Mesh Ewald electrostatics
pbc             = xyz           ; Periodic boundary conditions in all directions
```

Run it:

```bash
# compile
gmx grompp -f em.mdp -c BC400_neutralised.gro -p BC400.top -o em.tpr

# MD simulation
gmx mdrun -v -deffnm em

```

The `-v` flag prints progress to the terminal. Energy minimisation is successful if the potential energy decreases monotonically and the maximum force drops below the `emtol` threshold (500 kJ mol⁻¹ nm⁻¹). If it does not converge, the most common cause is a severe clash — check the structure in VMD and look for overlapping molecules.

---

## Step 5: NPT equilibration and production run

Simulations should be performed in the NPT (constant number of particles, pressure, temperature) ensemble, representing atmospheric conditions at room temperature.

### The .mdp file

The key parameters to use, following the published protocol:

```ini
; npt.mdp — NPT equilibration and production
integrator      = md
nsteps          = 25000000      ; 50 ns at 2 fs timestep
dt              = 0.002         ; 2 fs timestep

; Output
nstxout         = 50000         ; Save coordinates every 100 ps
nstvout         = 50000
nstenergy       = 5000
nstlog          = 5000

; Bonds
constraints     = h-bonds       ; Constrain H-bond lengths (LINCS algorithm)
constraint-algorithm = lincs

; Non-bonded
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.4
rvdw            = 1.4
coulombtype     = PME
pbc             = xyz

; Temperature coupling — Velocity-rescale thermostat
tcoupl          = V-rescale
tc-grps         = System
tau-t           = 0.1           ; Coupling constant (ps)
ref-t           = 300           ; Target temperature (K)

; Pressure coupling — C-rescale barostat, semi-isotropic
pcoupl          = c-rescale     ; previously Berendsen
pcoupltype      = semiisotropic ; Independent coupling in xy and z
tau-p           = 1.0           ; Coupling constant (ps)
ref-p           = 1.0  1.0      ; Reference pressure (bar), xy then z
compressibility = 4.5e-5  4.5e-5
```

The semi-isotropic pressure coupling is important here. Because the biochar slab occupies the centre of the box, you want the box to be able to change its size independently in the xy-plane (parallel to the surface) and in z (perpendicular to it). Using isotropic coupling would distort the slab geometry.

Run the simulation:

```bash
gmx grompp -f npt.mdp -c em.gro -p BC400.top -o npt.tpr

gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8   # adjust thread numbers for your system as needed

```

50 ns will not be possible on the laptop, you should move the simulation onto HPC cluster. The 24D_biochar repository provides example submission scripts for SLURM-based HPC systems.

---

## Step 6: Checking convergence — did the simulation equilibrate?

This is the most critical step, and the one students most often skip. Running time and equilibration time are not the same thing. A simulation that ran for 50 ns may have equilibrated in the first 15 ns or may not have equilibrated at all — you must check.

### RMSD

The most common first check is the Root Mean Square Deviation (RMSD) of the most mobile component in your system — in an adsorption study, that is the solute molecules. Calculate it with:

```bash
gmx rms -s npt.tpr -f npt.xtc -o rmsd_24D.xvg
```

Select the molecule group when prompted. Plot the output with `xmgrace rmsd_24D.xvg` or in Python. A plateau in the RMSD curve indicates that the system has reached a stable state. In our 2,4-D study, most systems plateaued within 30–35 ns.

### DynDen

RMSD alone is not sufficient for interfacial systems, where different components (the slab, the adsorbed layer, the bulk solution) equilibrate on very different timescales. To this end, we developed *DynDen* ([github.com/Erastova-group/DynDen](https://github.com/Erastova-group/DynDen)). *DynDen* computes pairwise cross-correlation between frames in a trajectory, so that you can see when the structural distribution of each component has converged, independently of the others.

Install and run *DynDen* following its own documentation. It produces a cross-correlation matrix: when the matrix shows a stable, consistent colour from a certain frame onwards (typically visible as a block of uniform colour), the system is equilibrated from that point. Only the frames after the equilibration point should be used for analysis.

---

## Step 7: Basic analysis — where are the molecules?

Once you have confirmed equilibration, a useful first analysis is the linear density profile along the z-axis. This tells you the distribution of each species perpendicular to the biochar surface:

```bash
gmx density -f npt.xtc -s npt.tpr -o density.xvg -sl 100 -b 30000  # start from 30 ns, use 100 slices along z
```

Select each group in turn when prompted (biochar, solute, water, ions). Plot the output. You should see:

- A high-density region in the centre of the box — the biochar slab
- A depletion of water near the surface
- An enrichment of solute molecules near the surface — this is the adsorbed layer
- Bulk solution composition far from the surface

If you see a roughly uniform distribution of solute throughout the box without enrichment near the surface, either adsorption is not occurring (which may be physically reasonable for some systems) or the simulation has not equilibrated.


### Counting adsorbed molecules

To quantify adsorption, count how many solute molecules lie within a defined distance of the biochar surface. The surface itself must be identified carefully: for a slab model, measuring distance from the centre of mass of the whole biochar points to the middle of the slab, not to either interface, and approximating the surface as a single z-plane ignores surface roughness and functional groups.

We use **[pytim](https://marcello-sega.github.io/pytim/) ITIM** (Identification of Truly Interfacial Molecules) to locate biochar atoms at the solvent interface, then measure the minimum distance from each solute molecule's centre of mass (COM) to the nearest ITIM surface atom. In our 2,4-D work, a molecule is counted as adsorbed if this distance is below **0.6 nm** — corresponding to the first peak in the radial distribution function between the solute's aromatic ring and surface carbon atoms. The same workflow is used in our more recent atrazine adsorption analysis.

The analysis has two steps: (1) identify surface atoms with ITIM, (2) count adsorbed solute molecules frame by frame.

#### Step 1 — Identify biochar surface atoms with pytim ITIM

ITIM rolls a probe sphere (radius `alpha`) over the biochar atoms and labels those in the outermost solvent-facing layer. Hydrogen needs an explicit van der Waals radius — pytim's default gives hydrogen zero radius and would exclude most surface H atoms.

```python
import MDAnalysis as mda
import numpy as np
import pytim
from MDAnalysis.lib import distances

STRUCTURE = "npt.tpr"          # or a solvated .gro from the equilibrated run
TRAJECTORY = "npt.xtc"
BIOCHAR_SELECTION = "resname 4G12"   # BC400; use "resname ?G*" for any biochar model

# ITIM probe radius (Angstrom) and layer to keep as "surface"
PYTIM_ALPHA = 1.3
SURFACE_LAYER = 1
BIOCHAR_RADII = {"C": 1.70, "O": 1.52, "H": 1.20}  # Angstrom; required for surface H

u = mda.Universe(STRUCTURE, TRAJECTORY)
biochar = u.select_atoms(BIOCHAR_SELECTION)

def run_itim(universe, biochar_ag):
    pytim.ITIM(
        universe, biochar_ag,
        alpha=PYTIM_ALPHA, normal="z", molecular=False,
        max_layers=3, radii_dict=BIOCHAR_RADII, include_zero_radius=True,
    )

def add_bonded_hydrogens(universe, biochar_ag, surface_ag, cutoff_nm=0.12):
    """Add H atoms bonded to surface C/O that ITIM may still miss."""
    heavy = surface_ag.select_atoms("name C* O*")
    hydrogens = biochar_ag.select_atoms("name H*")
    if heavy.n_atoms == 0 or hydrogens.n_atoms == 0:
        return surface_ag
    d = distances.distance_array(heavy.positions, hydrogens.positions, box=universe.dimensions)
    bonded_h = hydrogens[np.any(d < cutoff_nm, axis=0)]
    idx = np.unique(np.concatenate([surface_ag.indices, bonded_h.indices]))
    return biochar_ag[np.isin(biochar_ag.indices, idx)]

# Union of layer-1 atoms over equilibrated frames (atoms on either face at any point)
surface_indices = set()
for ts in u.trajectory[1500:]:   # from 30 ns onward; adjust to your equilibration point
    run_itim(u, biochar)
    surface = biochar[biochar.layers == SURFACE_LAYER]
    surface = add_bonded_hydrogens(u, biochar, surface)
    surface_indices.update(surface.indices)

surface_atoms = biochar[np.isin(biochar.indices, list(surface_indices))]
print(f"ITIM surface atoms: {surface_atoms.n_atoms}")

# Write a GROMACS index group for reuse (e.g. gmx rdf -n surface_ITIM.ndx)
with open("surface_ITIM.ndx", "w") as f:
    f.write("[ surface ]\n")
    gmx_ids = [a.index + 1 for a in surface_atoms]
    for i in range(0, len(gmx_ids), 15):
        f.write(" ".join(str(x) for x in gmx_ids[i:i + 15]) + "\n")
```

#### Step 2 — Count adsorbed solute molecules

For each frame, compute the minimum distance from each solute COM to any ITIM surface atom. Because the biochar slab is periodic in xy but not through the vacuum gap in z, use minimum-image convention in the xy-plane only (set `box[2] = 0` for distance routines).

```python
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib import distances

u = mda.Universe("npt.tpr", "npt.xtc")
solute = u.select_atoms("resname DCPAA")   # 2,4-D residue name

# Load pre-computed ITIM surface indices (1-based GROMACS atom numbers)
surface_ids = set()
with open("surface_ITIM.ndx") as f:
    in_group = False
    for line in f:
        line = line.strip()
        if line.startswith("[") and line.endswith("]"):
            in_group = line == "[ surface ]"
        elif in_group and line:
            surface_ids.update(int(x) for x in line.split())

surface_atoms = u.atoms[[a.index + 1 in surface_ids for a in u.atoms]]

ADSORPTION_CUTOFF_NM = 0.6   # from solute–surface C RDF first peak

def analysis_box(box):
    """Periodic in xy; direct z (no wrap through vacuum)."""
    b = np.asarray(box, dtype=float).copy()
    b[2] = 0.0
    return b

def com_dist_to_surface(com, surface_pos, box):
    com = np.asarray(com, dtype=float).reshape(1, 3)
    d = distances.distance_array(com, surface_pos, box=analysis_box(box))
    return float(d.min())

adsorbed_counts = []
for ts in u.trajectory[1500:]:
    box = u.dimensions
    surf_pos = surface_atoms.positions
    n_ads = 0
    for res in solute.residues:
        com = res.atoms.center_of_mass()
        if com_dist_to_surface(com, surf_pos, box) < ADSORPTION_CUTOFF_NM:
            n_ads += 1
    adsorbed_counts.append(n_ads)

print(f"Mean adsorbed molecules: {np.mean(adsorbed_counts):.1f} ± {np.std(adsorbed_counts):.1f}")
print(f"Mean adsorbed fraction: {100 * np.mean(adsorbed_counts) / solute.n_residues:.1f}%")
```

This gives the mean number (and fraction) of adsorbed molecules over the equilibrated window. Divide by the solvent-accessible surface area of the biochar (from `gmx sasa`; see [Building Biochar Molecular Models from Scratch](biochar_model_building_tutorial.md), Section 7.3) to obtain adsorption per unit area. Choose the adsorption cutoff from the relevant solute–surface RDF — for 2,4-D aromatic ring to surface carbon, 0.6 nm is appropriate; other solutes may need a different value (our atrazine work used 0.85 nm).

---

## Step 8: Visualisation in VMD

VMD is invaluable for building intuition about what is happening. To view the trajectory:

```bash
vmd BC400_neutralised.gro npt.xtc -stride 100

```
Loading of the whole trajectory can take some time, and you are unlikely to need to visualise every single step. Instead, choose `-stride 100` option to skip every 100 steps.

Useful representations for biochar systems:
- Biochar atoms: **Licorice** representation, coloured by element (C = cyan, O = red, H = white)
- Solute molecules: **CPK** (van der Waals spheres), coloured distinctly
- Ions: **VDW** spheres
- Water: **QuickSurf** or **Transparent Surface** — making it partially transparent lets you see through the water to the interface

To colour different adsorption modes, you can use the Tcl console to select atoms by distance from the surface and colour them accordingly.

---

## Common problems and how to fix them

**"Fatal error: number of coordinates in coordinate file does not match topology"**
The number of atoms in your `.gro` file does not match what is listed in the `[ molecules ]` section of your `.top` file. This usually happens after `gmx insert-molecules` or `gmx genion`. Check the terminal output of those commands and update the `[ molecules ]` section manually.

**Energy minimisation does not converge**
There is likely a severe steric clash in the starting configuration. Do not ignore and proceed to the next step! Open the output `.gro` of the last minimisation step in VMD and look for atoms that are sitting on top of each other. The most common cause is solute molecules placed inside the biochar material by `gmx insert-molecules`. Re-run `gmx insert-molecules` with a smaller number of molecules, or verify the vdW exclusions are working correctly.

**Simulation crashes immediately with a "LINCS warning"**
Hydrogen bonds are flying apart. This usually means the timestep is too large for the current geometry, or that there is still a clash. Try running a brief NVT step (constant volume, no pressure coupling) at a very short timest
ep (0.5 fs) before the full NPT run to let the system relax.

**No adsorption visible after 5 ns**
This is usually just insufficient simulation time. Equilibration in these systems takes 30+ ns. Run longer. 

---

## What comes next

Once you have a working simulation and confirmed equilibration, we suggest the following analysis framework:

- **Adsorption isotherms** by counting adsorbed molecules across a series of concentrations and plotting adsorbed amount vs equilibrium concentration (Wood, Mašek & Erastova 2025, Figure 3)
- **Adsorption geometry** from distances between solute functional groups and surface atoms, and orientation angles relative to the surface normal
- **Cation bridging analysis** via radial distribution functions between counterions and surface/solute oxygen atoms
- **Convergence validation** with *DynDen*
- **Diffusion in the pore space** informing on the impact of porosity on solute retention


---

## Referencing

We hope this tutorial was useful for you! If you use our models or approach in your published work, **please cite our works** that enabled this tutorial.

Adsorption studies:

1. Wood R, Mašek O & Erastova V. *Unravelling 2,4-D – biochar interactions by molecular dynamics: adsorption modes and surface functionalities.* Biomas and Bioenergy 210, 109080 (2025). DOI: [10.1016/j.biombioe.2026.109080](https://doi.org/10.1016/j.biombioe.2026.109080)

2. Ngambia A, Gavrilova A, Huang H, Lyu Z, Mašek O, Graham M & Erastova V. *Decoupling Precipitation and Surface Complexation during Mn (II) Removal by Biochar via Experiments and Atomistic Simulations.* arXiv preprint (2026). DOI: [10.48550/arXiv.2603.22144](https://doi.org/10.48550/arXiv.2603.22144)

Biochar model development:

3. Wood R, Mašek O & Erastova V. *Developing a molecular-level understanding of biochar materials using public characterization data.* Cell Reports Physical Science 5(7), 102036 (2024). DOI: [10.1016/j.xcrp.2024.102036](https://doi.org/10.1016/j.xcrp.2024.102036)

4. Wood R, Mašek O & Erastova V. *Developing realistic molecular models of biochars.* Cell Reports Physical Science 5(7), 102037 (2024). DOI: [10.1016/j.xcrp.2024.102037](https://doi.org/10.1016/j.xcrp.2024.102037)

5. Ngambia A, Mašek O & Erastova V. *Development of biochar molecular models with controlled porosity.* Biomass and Bioenergy 184, 107199 (2024). DOI: [10.1016/j.biombioe.2024.107199](https://doi.org/10.1016/j.biombioe.2024.107199)

Convergence analysis:

6. Degiacomi MT, Tian S, Greenwell HC & Erastova V. *DynDen: Assessing convergence of molecular dynamics simulations of interfaces.* Computer Physics Communications 269 (2021). DOI: [10.1016/j.cpc.2021.108126](https://doi.org/10.1016/j.cpc.2021.108126)

PyTIM:

7. Sega M, Hantal G, Fábián B & Jedlovszky P. *Pytim: A python package for the interfacial analysis of molecular simulations.* Journal of Computational Chemistry 39(25), 2118–2125 (2018). DOI: [10.1002/jcc.25384](https://doi.org/10.1002/jcc.25384). Documentation: [marcello-sega.github.io/pytim](https://marcello-sega.github.io/pytim/)

**GitHub repositories (Erastova group):**

- [github.com/Erastova-group/Biochar_MolecularModels](https://github.com/Erastova-group/Biochar_MolecularModels) — BC400, BC600, BC800 models + experimental database
- [github.com/Erastova-group/Porous_Biochars_Models](https://github.com/Erastova-group/Porous_Biochars_Models) — microporous biochar models
- [github.com/Erastova-group/24D_biochar](https://github.com/Erastova-group/24D_biochar) — complete simulation files for 2,4-D adsorption study
- [github.com/Erastova-group/Mn_Biochar](https://github.com/Erastova-group/Mn_Biochar) — Mn(II) interactions study
- [github.com/Erastova-group/DynDen](https://github.com/Erastova-group/DynDen) — convergence analysis tool
- [github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques](https://github.com/Edinburgh-Chemistry-Teaching/MD_ResearchTechniques) - workshop materials for setting up MD simulations on HPC

---
## Happy simulating!
