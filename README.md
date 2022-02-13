# Install
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal
conda install pip
pip install Jinja2
pip install -U pytest

# Run the demo data
snakemake -j 4 --use-conda

# Unit tests
snakemake -j 4 --use-conda --generate-unit-tests 
pytest .tests/
