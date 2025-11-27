# PRECISE

Official implementation of 'Learning a PRECISE model for small-molecule binding'

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
