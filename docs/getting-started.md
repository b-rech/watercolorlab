# Getting started

## Installation with Conda

Clone the repository and create the project environment:

```bash
git clone https://github.com/b-rech/watercolorlab.git
cd watercolorlab
conda env create -f environment.yml
conda activate watercolorlab
pip install .
```

Alternatively, install the package from GitHub after installing its
dependencies:

```bash
pip install git+https://github.com/b-rech/watercolorlab.git
```

## Examples

The repository includes Jupyter notebooks demonstrating:

- [CDOM processing](https://github.com/b-rech/watercolorlab/blob/main/examples/cdom_processing_wpi.ipynb)
- [Particulate processing](https://github.com/b-rech/watercolorlab/blob/main/examples/tr_processing_perkinelmer.ipynb)

Install JupyterLab in the project environment if needed:

```bash
conda install -c conda-forge jupyterlab
jupyter lab
```
