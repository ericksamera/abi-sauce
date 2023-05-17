# abi-sauce
This repo is intended to hold the scripts for my "introductory bioinformatics" toolset.

## Installation
### Requirements
* `biopython` >= `3.11`
* `streamlit` >= `1.22.0`
* `streamlit-ext` >= `0.1.7`
* `plotly` >= `5.14.1`

Installation using only `conda` has some errors. I'd recommend this order of installation:
```
conda create -n streamlit -c bioconda biopython plotly pip
conda activate streamlit
pip install streamlit streamlit-ext 
```

## Usage
Like any other Streamlit app, running the apps in this repo use the following command:
```
streamlit run main.py
```
