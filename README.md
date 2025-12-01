![](bannerv3.png)
# PRECISE

Contemporary drug-screening pipelines struggle to perform two operations simultaneously: (a) leveraging extremely large databases of synthesizable ligands such as ZINC and Enamine, and (b) conducting targeted screening against a specific protein receptor site. Here, we introduce PRECISE, a machine-learning framework that achieves both by reimagining protein–ligand binding as an association between quantized ligand features and a protein-surface representation. PRECISE accepts a protein structure (via PDB) and a ligand SMILES string, internally converts the structure into a protein surface representation enriched with physicochemical and geometric features, and transforms the ligand into a discretized representation, ultimately integrating them to produce a per-vertex binding likelihood. For each surface-mesh vertex, PRECISE computes four key properties: (a) hydropathy, (b) hydrogen-bond donor–acceptor characteristics, (c) surface electrostatics, and (d) a shape index capturing local curvature. Ligand discretization is performed using CoNCISE, an ML model developed by our lab in 2024 (Erden et al., 2024), which enables PRECISE to efficiently screen entire ligand libraries at scale. This fundamental operation supports higher-level operations such as (a) ligand-specific binding-site prediction on proteins and (b) site-specific identification of potential ligand binders for any target, and it further exhibits robust zero-shot generalization to challenging systems including metalloproteins and multi-chain complexes.

In the PRECISE paper, we demonstrate that PRECISE performs on par with state-of-the-art co-folding and ML-based docking approaches in identifying ligand-binding sites. We also show that our screening pipeline, which employs an efficient Markov Chain Tree Search (MCTS) strategy, can capitalize on the scale and chemical diversity of large ligand databases to identify Vina-inferred strong binders while requiring only minimal docking.

### Getting Started

To install the environment for inference, you need to install via the following codes:

```
micromamba create -f environment.yml
micromamba activate precise

pip install torch-geometric \
  torch-scatter \
  torch-sparse \
  torch-cluster \
  torch-spline-conv \
  -f https://data.pyg.org/whl/torch-2.5.1+cu121.html

pip install -e .
```

Also, for the surface calculation, you need to install the tools via the following codes:

```
precise install-tools
```

### Examples

Here are some examples of how to use our tools:
```
# If you want to calculate the surface of a single chain protein
precise create-surface --chain-id A --pdb-path ./example/2QCS.pdb

# If you want to calculate the surface of a multi-chain protein
precise create-surface --chain-id A,B --pdb-path ./example/2QCS.pdb

# For the pocket detect function
precise detect-pockets --pdb-path ./example/3WBB.pdb --smiles 'CCO'

# For the virtual screening function
# Example with SDF:
precise virtual-screen \
    --pdb-path ./example/3WBB.pdb \
    --output-dir results/ \
    --center-sdf ./example/3WBB_NAP.sdf \
    --search-mode wide \
    --no-codes-to-consider 500 \
    --no-smiles-per-codes 100
# --db-path zinc.db 

# Example with coordinates:
precise virtual-screen \
    --pdb-path ./example/3WBB.pdb \
    --output-dir results/ \
    --center-x 10.5 \
    --center-y 20.3 \
    --center-z 15.7 \
# --db-path zinc.db
```
