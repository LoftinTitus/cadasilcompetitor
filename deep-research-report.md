# ML‑Driven Multiscale Pipeline for HS‑Selective Peptide Design in a Simulated Neurovascular Environment

## Executive summary

A computational pipeline to design a peptide with **selective binding to heparan sulfate (HS)** in a **neurovascular (BBB/NVU) extracellular matrix (ECM) context** must be **multiscale**: molecular binding (affinity, selectivity, kinetics) is necessary but not sufficient, because the functional question is whether the peptide achieves **adequate HS site occupancy over time** under **diffusion, washout/clearance, and barrier constraints**. 

Two realities drive the architecture:

- **HS is not a single structure**. HS chains vary in **sulfation pattern, uronic acid epimerization, and domain organization**, and these features alter protein/ligand interactions; rare yet biologically important motifs such as **3‑O‑sulfation** can dominate specific interactions (including anticoagulant biology). 
- In neurovascular and neurodegenerative contexts, **HSPGs/HS can mediate binding and cellular uptake of proteopathic aggregates** (strong evidence exists for tau aggregate binding/uptake via HSPGs), so HS‑binding peptides require **explicit negative design** against misfolded aggregate panels. 

A rigorous, scalable approach is a **Python‑first orchestration layer** that automates: (a) a curated HS/CS/DS oligo panel, (b) GAG‑aware docking as triage with benchmark calibration (protein–GAG docking performance varies widely across tools), (c) explicit‑solvent/ion MD refinement (salt sensitivity is central), (d) uncertainty‑calibrated ML surrogate models trained on simulation + sparse experimental kinetics labels, (e) multiobjective Bayesian optimization to maintain a Pareto frontier, and (f) a compartment or reaction–diffusion neurovascular transport model that maps kinetics into time‑dependent occupancy. 

Because you **do not currently have a paper**, this report treats all “paper‑provided” fields as **TBD**, supplies a **plausible default v0 input set** (explicitly labeled as assumptions), and frames implementation as a **config‑driven GitHub‑ready project skeleton** that can be updated when a real paper/spec becomes available. 

## Objectives and deliverables

### Core objectives

**Design objective:** identify peptide candidates that bind HS variants strongly enough to achieve target engagement but selectively enough to avoid broad polyanion binding (CS/DS, nucleic acids, membranes) and without driving adverse biological effects. 

**System objective:** maximize **time‑dependent HS site occupancy** in a specified neurovascular compartment (glycocalyx vs vascular basement membrane vs perivascular space), accounting for diffusion and clearance. BBB transport modeling literature emphasizes that delivery/uptake is governed by transport mechanisms and can be represented with mathematical (compartment/vessel) models that couple kinetics and transport. 

**Negative design objective:** minimize binding/adsorption to misfolded aggregates (e.g., tau ± Aβ ± α‑syn), motivated by evidence that HSPGs can be critical mediators of aggregate binding/uptake (tau aggregates in particular). 

### Hard constraint gates

- **BBB/NVU safety gate:** no barrier integrity compromise (TEER/permeability non‑worsening), with TEER measurement known to be highly sensitive to technical parameters and requiring standardized reporting. 
- **Anticoagulant‑risk gate:** explicitly avoid unintended interactions with anticoagulant HS/heparin motifs; 3‑O‑sulfated glucosamine is key to high‑affinity antithrombin binding and anticoagulant biology, so motifs/variants containing 3‑O‑sulfation should be included in screening as “do‑not‑engage” negatives (unless the therapeutic mechanism intends otherwise). 
- **Developability gate:** reduce peptide self‑aggregation and insolubility using established predictors (TANGO, AGGRESCAN, CamSol) as first‑pass filters before expensive modeling; these tools are widely cited for aggregation/solubility screening.  
- **Toxicity gate (computational filter):** use an explicit toxicity/hemolysis risk predictor as a screening penalty (not a substitute for wet lab), supported by literature demonstrating modern sequence+structure ML toxicity models. 

### Deliverables appropriate for a GitHub README

- A **config‑driven pipeline** that generates candidate peptides, computes sequence descriptors, runs multi‑fidelity scoring (cheap → expensive), trains a surrogate with uncertainty, and proposes next candidates via multiobjective Bayesian optimization. 
- Reproducible **data manifests** and versioned panels for HS/CS/DS oligos, off‑target proteins, and aggregate structures using programmatic access to primary/official sources. 

## Inputs and recommended default starting set

### Inputs that must be extracted from the user’s paper

The following are **TBD until a paper/spec exists** and should be represented as required fields in configuration files:

- Peptide sequence(s), modifications, length bounds, prohibited residues/chemistries  
- Target neurovascular compartment and access route (luminal vs abluminal)  
- Target HS structural hypothesis (which sulfation motifs and chain lengths matter)  
- Required concentration and exposure time scales for modeling and assays  
- Intended mechanism (competition/displacement vs blockade of aggregate docking vs occupancy thresholding)

These are not “nice to have”; they determine HS panel composition, modeling geometry, and evaluation metrics. 

### Recommended default v0 inputs for initial computational work

Everything below is an **assumption** meant to enable early pipeline development and benchmarking; replace once a real paper/spec exists.

#### Candidate peptide sequences and variants

Net charge is shown as an **approximate integer** assuming Lys/Arg = +1, Asp/Glu = −1, and His treated as ~0 at pH 7.4 (a simplification suitable for early screening but not for final electrostatics). *(Assumption.)* 

| ID | Sequence | Length | Approx. net charge | Modifications | Rationale (assumption) |
|---|---:|---:|---:|---|---|
| PaperLead‑1 | **TBD (from paper)** | TBD | TBD | TBD | Replace with paper-defined lead(s). |
| P1 | AKRKRQGK | 8 | +5 | none | Short HS‑binding motif seed with polar spacer (Q) to reduce pure charge clustering. |
| P2 | GRRGRKQK | 8 | +5 | none | Similar charge with altered spacing; probes sensitivity to motif ordering. |
| P3 | KRGKRRQA | 8 | +5 | none | Adds alanine spacer; tests compact basic patch behavior. |
| P4 | RRKQGRKR | 8 | +6 | none | Slightly higher charge density; use as “too sticky?” control. |
| P5 | AKKQKAKKQKAKKQKA | 16 | +9 | none | Longer repeat with spaced basics; stresses selectivity and uptake risk penalties. |
| P6 | GKKRGRRKRRK | 11 | +9 | none | Dense basic cluster positive control for HS affinity but high off‑target risk. |
| P7 | RKRGRQKRRKA | 11 | +8 | none | Similar charge with polar insertions for tunable binding. |
| P8 | GKRRKAKRGRR | 11 | +8 | none | Alternating K/R blocks; tests electrostatic patterning effects. |
| P9 | KRGRRKQKRGR | 11 | +8 | none | Basic‑rich motif with Q spacer to reduce nonspecific adsorption. |
| P10 | GKKRGKKEGKKRGK | 14 | +8 | none | Introduces single acidic residue (E) to probe specificity vs “global polyanion binding.” |
| P11 | GKKQGKQEGKKQGKQ | 15 | +5 | none | Moderate charge; tests whether lower charge can preserve HS selectivity. |
| P12 | GHHKHKQKHKQKHKQK | 16 | +7 | none | Histidine‑rich pH‑tunable design (screen across pH as needed). |
| P13 | GSSGSSGSSGSSGSS | 15 | 0 | none | Neutral negative control for nonspecific adsorption and assay baselines. |
| P14 | AKRGRKRRKQGA | 12 | +7 | none | Mid‑length candidate; balanced charge and spacing. |
| P15 | AKRKRGRRKQGAKR | 14 | +9 | none | Higher charge mid‑length candidate for Pareto tradeoff tests. |
| P16 | AKRGRKQKRRQGAKRGRK | 18 | +11 | none | Long, high-charge “upper envelope” stress test (should be penalized by safety/selectivity). |
| P17 | AKRGRQKQERKQGRQKQERK | 20 | +8 | none | Longer candidate with acidic residues controlling nonspecific binding. |

**Why motif lengths can be short but peptides must be longer:** empirical motif descriptions (e.g., Cardin–Weintraub consensus patterns) identify short basic patterns, but binding specificity depends on 3D presentation, spacing, and context; reviews reiterate that motifs are useful heuristics, not guarantees of specificity. 

#### HS/CS/DS oligo panel

Because IDs depend on exact structures, the **GlyTouCan accession IDs are marked TBD and should be resolved programmatically** during data ingestion by registering or querying precise WURCS strings. 

| Internal ID | GlyTouCan accession | GAG class | Length (dp) | Sulfation pattern | Rationale |
|---|---|---|---:|---|---|
| HS‑dp4‑NAc | TBD | HS | 4 | no sulfation (heparosan‑like) | Negative control for “charge‑driven” binding. |
| HS‑dp4‑NS | TBD | HS | 4 | N‑sulfation only | Minimal sulfation motif; tests N‑S dependence. |
| HS‑dp4‑NS‑2S | TBD | HS | 4 | N‑S + 2‑O‑S | Probes IdoA2S‑linked binding preference. |
| HS‑dp4‑NS‑6S | TBD | HS | 4 | N‑S + 6‑O‑S | Probes 6‑O‑S–driven binding preferences. |
| HS‑dp6‑mixed | TBD | HS | 6 | mixed GlcA/IdoA and mixed sulfation | Representative HS heterogeneity (not heparin‑like). |
| Heparin‑dp6‑highS | TBD | heparin (Hp) | 6 | highly sulfated (Hp‑like) | High-affinity positive control; known to differ structurally from HS and may overestimate binding. |
| HS‑dp5‑AT‑motif | TBD | HS | 5 | includes 3‑O‑sulfation | Anticoagulant‑risk sentinel; 3‑O‑S is key for antithrombin binding. |
| CS‑A‑dp4 | TBD | CS | 4 | GalNAc4S dominant (CS‑A) | Negative panel for selectivity.  |
| CS‑C‑dp4 | TBD | CS | 4 | GalNAc6S dominant (CS‑C) | Tests sulfation‑position sensitivity.  |
| CS‑E‑dp4 | TBD | CS | 4 | GalNAc4S6S (CS‑E) | Highly sulfated CS negative; detects “high sulfate = bind” failure modes. |
| DS‑dp4 | TBD | DS | 4 | IdoA + GalNAc4S | DS negative; tests iduronate‑driven conformational effects.  |

#### Neurovascular compartments and ECM components

The vascular basement membrane is described as a network primarily composed of **laminin, collagen IV, nidogen, and heparan sulfate proteoglycans**, and it supports interactions between brain endothelial cells, pericytes, and astrocyte endfeet—making it a biologically grounded target compartment for HS‑mediated interventions. 

| Compartment | Cell types present | Core ECM components (assumption) | HS accessibility | Modeling notes |
|---|---|---|---|---|
| Blood lumen | none | none | indirect | Include convection/washout boundary if flow is modeled. |
| Endothelial luminal glycocalyx | brain endothelial cells | HSPGs + glycocalyx matrix | high | Screen uptake/internalization risk because HS/HSPGs can mediate binding and entry for cationic cargos.  |
| Endothelial barrier layer | brain endothelial cells | tight junction system + basal lamina interface | moderate | Parameterize permeability using TEER/permeability constraints. |
| Vascular basement membrane | endothelium + pericytes + astrocyte endfeet adjacency | laminin, collagen IV, nidogen, HSPGs (perlecan/agrin) | high | Prime site for diffusion‑limited binding; model fixed binding sites + diffusion in ECM. |
| Perivascular space | astrocytes + pericytes nearby | interstitial ECM | variable | Apply safety constraints for unintended parenchymal exposure. |

### Default concentration/time scale assumptions

Because you requested “no invented values,” the repo should store these as **configurable priors** rather than hard-coded constants. For a v0 computational project, a reasonable practice is to define **log‑spaced concentration grids** and simulate **minutes→hours** occupancy under a small set of clearance regimes; BBB transport modeling literature supports using mechanistic models with parameter estimation rather than fixed universal constants. 

## Modeling scope and computational methods

### Molecular modeling stack

**Docking (triage):** protein–GAG docking performance is highly tool‑dependent; a benchmarking study evaluated eight docking programs on 28 protein–GAG complexes and motivates using docking primarily to generate poses and hypotheses, not as final rankers. 

**GAG‑aware docking:** specialized methods like **GAG‑Dock** were developed/validated for predicting poses of protein‑bound GAGs, providing a better starting point than generic docking for sulfated polysaccharides. 

**Atomistic MD refinement:** explicit salt/water MD is central because electrostatics and ion pairing dominate GAG interfaces; MD setup for glycans is supported by **CHARMM‑GUI Glycan Modeler**, and engine choice can be Python‑controlled with high performance via OpenMM. 

**Force field and sulfate/ion pairing sensitivity:** sulfated groups and counter‑ion interactions can be force‑field sensitive; published work extends GLYCAM parameters for heparin‑like GAGs and benchmarks sulfate/sulfamate parameterization and ion pairing behavior. 

**Heparin ≠ HS:** solution scattering studies conclude that HS adopts conformations significantly distinct from heparin (HS can be longer/more bent and shows different flexibility), so relying on heparin alone can bias binding expectations.

**Enhanced sampling / kinetics (finalists):** add enhanced sampling when unbinding/binding events are rare on standard MD timescales; reviews cover umbrella sampling/metadynamics families, and weighted ensemble methods are reviewed for rare-event kinetics estimation. 

image_group{"layout":"carousel","aspect_ratio":"16:9","query":["heparan sulfate sulfation pattern diagram","vascular basement membrane laminin collagen IV nidogen heparan sulfate proteoglycans schematic","blood brain barrier neurovascular unit schematic"],"num_per_query":1}

### Neurovascular transport and occupancy modeling

A recommended progression is:

- **Tier A (compartment ODE):** fast sweeps for identifiability and “which parameter matters?” reasoning.  
- **Tier B (1D reaction–diffusion PDE):** vessel wall thickness direction with fixed binding sites, diffusion, and clearance; enables prediction of gradients and diffusion-limited binding regimes. 
- **Tier C (optional ABM/3D):** if heterogeneous ECM patches or cell remodeling are central, use agent-based modeling; PhysiCell provides an open framework for multicellular simulation with diffusion fields.

### Multiscale coupling principle

The coupling mechanics should be explicit and testable:

1. Molecular binding estimates (KD/kon/koff or surrogate proxies) → effective binding parameters per HS variant. 
2. Transport model uses these parameters to compute **HS occupancy vs time** in the target compartment under diffusion/clearance. 
3. Optimization loop targets occupancy and selectivity simultaneously (multiobjective scoring).

### Candidate method comparison

| Method family | Primary question answered | Strengths | Key failure modes | When to use |
|---|---|---|---|---|
| Docking (benchmark‑calibrated) | pose hypotheses | high throughput | unreliable rank ordering for flexible, highly charged GAGs | early triage |
| GAG‑Dock (specialized) | GAG pose prediction | validated on GAG complexes | still needs calibration to your chemistry | early–mid |
| Atomistic MD | stability + ion/water mediation | salt/solvent realism; mechanistic features  | sampling cost; force-field sensitivity | refinement |
| Enhanced sampling / WE kinetics | ΔΔG and/or rates | resolves rare events | method complexity | finalists |
| Reaction–diffusion PDE | occupancy over time | interpretable, scalable; sensitivity analysis  | parameter identifiability | always |
| Agent‑based (PhysiCell) | heterogeneity, cell behaviors | emergent spatial effects | parameter explosion | only if justified |

## ML optimization loop and uncertainty management

### ML tasks and model classes

**Sequence → property predictors:** predict HS affinity proxy, HS vs CS/DS selectivity, aggregate-binding penalty, and developability/toxicity filters. Aggregation predictors (TANGO, AGGRESCAN) and solubility predictors (CamSol) are established choices for large-scale filtering. 

**Transfer learning:** use protein language model representations; ESM‑2 demonstrates strong sequence-to-structure signal capture at scale, and PLMFit benchmarks transfer learning strategies for protein engineering. 

**Generative models:** use constrained generation rather than unconstrained sequence sampling; if you later incorporate structure context, ligand/context-conditioned design methods (e.g., LigandMPNN) support conditioning on non-protein components—conceptually relevant to peptide–HS contexts (though HS is a polymeric ligand). 

### Active learning and Bayesian optimization

Bayesian optimization is appropriate because experimental labels (SPR kinetics, selectivity panels, TEER/permeability) are expensive; BoTorch provides an efficient MC-based BO framework, and qEHVI provides a practical acquisition function for parallel, constrained multiobjective optimization.

### Uncertainty quantification

Uncertainty must be evaluated, not assumed. A benchmark study on uncertainty quantification methods for protein engineering evaluates UQ methods on FLIP tasks and assesses utility under distribution shift and for active learning/optimization. 

### Suggested Pareto trade-off plots (no data required)

Design Pareto charts around decisions you actually need to make:

- HS affinity proxy vs HS/CS selectivity proxy  
- HS affinity proxy vs aggregate-binding penalty (tau panel)  
- Predicted HS occupancy time-above-threshold vs TEER/permeability risk penalty  
- Residence time proxy (1/koff) vs washout sensitivity (from transport model)

These make “best peptide” a transparent multiobjective argument rather than a single-score gamble.

## Data sources, evaluation metrics, and validation strategy

### Priority data sources

Use primary/official sources as the backbone:

- entity["organization","RCSB Protein Data Bank","structure database"] APIs (Search/Data) for structures, ligands, and metadata.  
- entity["organization","UniProt Knowledgebase (UniProtKB)","protein sequence database"] API endpoints for protein sequences and annotations. 
- entity["organization","GlyTouCan","international glycan repository"] for stable glycan accessions and WURCS-encoded structure registration/query. 
- entity["organization","GlyGen","glycoinformatics knowledgebase"] for integrated glyco/protein relationships and harmonized metadata. 

#### Data source comparison table

| Source | What to pull | Why it matters | Practical note |
|---|---|---|---|
| PDB APIs | protein–GAG complexes, fibrils/aggregates | benchmarking + negative design structures  | prefer API-based manifests for reproducibility  |
| UniProtKB APIs | sequences/annotations of HS-binding proteins | build off-target protein panels | store accessions + versioned query results |
| GlyTouCan | accession IDs + WURCS | standardize HS/CS/DS library  | resolve IDs programmatically from exact structures |
| GlyGen | integrated glyco/protein metadata | faster joins and curated context | use as enrichment layer, not sole truth |

### Evaluation metrics

Minimum viable metric suite:

- **Affinity/selectivity:** KD (or proxies) across HS panel and CS/DS panel (ratio/ΔΔG). 
- **Kinetics:** kon/koff and derived residence time; critical for occupancy under washout.
- **Aggregate binding penalty:** adsorption/binding against tau aggregate panels (and others if specified).  
- **Developability:** aggregation propensity and solubility metrics (TANGO, AGGRESCAN, CamSol).  
- **Safety proxies:** TEER/permeability penalty; aPTT penalty; MMP interference penalty. 
- **System-level:** simulated HS site occupancy vs time (time above threshold). 

### Validation assays and what they calibrate

| Assay | Primary output | Pipeline component calibrated | Notes |
|---|---|---|---|
| SPR | kon/koff/KD | molecular kinetics + surrogate ground truth | supports transport model parameterization |
| Displacement assay | competition curves | mechanism validation | confirms “competitive occupancy” hypothesis |
| GAG selectivity panel | HS vs CS vs DS ratio | selectivity metric | prevents “polyanion binder” winners  |
| TEER/permeability | barrier integrity | BBB safety and permeability parameters | TEER is sensitive to setup; follow critical review guidance  |
| Cytokine panels | inflammatory activation | immunotoxicity gate | interpret with NVU model choice |
| aPTT | coagulation impact | anticoagulant-risk gate | PTT reflects intrinsic/common pathways; used clinically to assess coagulation |
| MMP activity assay | protease activity change | protease interference gate | fluorogenic peptide substrate methods are standard  |
| Aggregate-binding assay | adsorption/colocalization | negative design objective | motivated by HSPG-mediated tau uptake evidence  |

## Implementation plan for GitHub, compute, and milestones

### Python-first implementation guidance

Python-first is the correct default because: (a) primary data ingestion is API-driven, (b) BO/UQ stacks are Python-native, and (c) major simulation/PDE stacks provide Python entry points with optimized backends. 

Suggested repo backbone (README-friendly):

- `configs/` (all assumptions and paper-derived values live here)
- `data/manifests/` (versioned panels: HS/CS/DS, off-target proteins, aggregate structures)
- `peptides/` (candidate generation + filtering + descriptors)
- `scoring/` (multi-fidelity scoring wrappers; objective function)
- `models/` (ML surrogates, UQ, BO)
- `transport/` (Tier A ODE, Tier B PDE specs)
- `reports/` (figures, Pareto plots, run summaries)
- `tests/` (unit tests for configs, descriptors, filtering)
- `docs/` (assumptions and decision logs)

### Compute requirements (expressed as scalable constraints)

- Docking: CPU-parallel; scale with candidates × oligos; calibrate using protein–GAG docking benchmarks before trusting scores. 
- MD: GPU recommended; explicit salt/water; allocate additional budget for sulfate/ion pairing sensitivity tests. 
- PDE/ODE (Tier A/B): lightweight compared to MD; suitable for large sweeps and sensitivity analyses. 

