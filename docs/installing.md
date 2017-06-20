# How to Use this Project

## Contents 

- [Dependancies](#Dependancies)
- [Setting Up](#setting-up)

<hr>

## Dependancies

In order to use this project a list of dependancies will first need to be installed or setup. This also depends on what you plan to use the project for. Procomp can be broken down into the following categories with their associated dependancies:

- Wrangler ( querying data from ensembl )
    - virtualenv or conda
    - pycogent
    - sqlalchemy
    - mysql-connector-python
- Process ( cleaning datasets and preparing them for analysis )
    - pandas
- Procomp ( analysis on the cleaned datasets )
    - biopython

all packages used are available to install using pip (a common python package installer)

<hr>

## Setting Up

To use procomp at the base level of setup involves: 

<h4>1. Downloading the project from github to a desired folder</h4>
    
```bash
$ cd "path-to-folder"
$ git clone "link-to-repo"
```
<h4>2. Adding the project directory to PATH.</h4>
    
To do this locate what shell you are using
```bash
$ echo $0
```
Then add the following lines to the rc file for that shell. (i.e. the .bashrc file is located in the home directory ~/)
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/procomp
```
    
<h4>3. Using Procomp</h4>
To Use procomp with the previous steps completed simply write scripts in a similar way to those provided in the examples folder. Our method for importing procomp follows the convention:

```python
import procomp as pc
import process as ps
import wrangler as wg
```

<hr>
