# How to Use this Project

## Contents 

- [Dependancies](#Dependancies)
- [Setting Up](#setting-up)
    1. [Downloading the Project](#1-downloading-the-project-from-github-to-a-desired-folder)
    2. [Adding the project directory to PATH](#2-adding-the-project-directory-to-path)
    3. [Basics of Using Procomp](#3-basics-of-using-procomp)

<hr>

## Dependancies

In order to use this project a list of dependancies will first need to be installed or setup. This primarily depends on what you plan to use the workflow for. Procomp can be broken down into the following categories with their associated dependancies:

- Example 1-Wrangling ( querying data from ensembl )
    - virtualenv or conda
    - pycogent
    - sqlalchemy
    - mysql-connector-python
- Example 2-Process ( cleaning datasets and preparing them for analysis )
    - pandas
- Example 3-Procomp ( analysis on the cleaned datasets )
    - biopython
- Example 4-Presentaion ( Visualization of data)
    - matplotlib
    - pandas

all packages used are available to install using pip (a common python package installer)

<hr>

## Setting Up

To use procomp at the base level of setup involves: 

#### 1. Downloading the project from github to a desired folder
    
```bash
$ cd "path-to-folder"
$ git clone "link-to-repo"
```
#### 2. Adding the project directory to PATH
    
To do this locate what shell you are using
```bash
$ echo $0
```
Then add the following lines to the rc file for that shell. (i.e. the .bashrc file is located in the home directory ~/)
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/procomp
```
    
#### 3. Basics of Using Procomp
To Use procomp with the previous steps completed simply write scripts in a similar way to those provided in the examples folder. Our method for importing procomp follows the convention:

```python
import procomp as pc
import process as ps
import wrangler as wg
```

<hr>
