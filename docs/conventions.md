# Coding Conventions

## Contents 

- [Introduction](#introduction)
- [Conventions](#conventions)
    - [Naming](#naming)
    - [Variables](#variables)
    - [Docstrings](#docstrings)

<hr>

## Introduction

In this document we highlight some of theconventions that are seen in the code which make the project more coherent as well easier to understand when looking at the code.

## Conventions
This project uses a few core principles to organize and maintain consistency within the project.
- various uses of the project are shown in the examples folder where each experiment contains the high level use of our code as well as the results that experiment returns.
- Folders and Results are capitalized and use camel case.
- variables, function names, resource files and scripts are lowercase and underscore delimited

#### Naming
|Value    |Meaning |
|---------|--------|
|L        |any variable with "L" in the name represents a python list |
|spe      |any variable with "spe" in the name represents a species |
|_x       |underscore in front of a variable or function indicates a private object which should not be called by the user |
|gn |any variable with "gn" in the name represents a gene |
|tr |any variable with "gn" in the name represents a transcript |
|pr |any variable with "gn" in the name represents a protein |
|R |represents a "regenerating species" |
|NR |represents a "non-regenerating species" |
|fl |represents a file |
|dir |represents a directory |
|cond |represents a condition for a particular outcome | 


#### Variables
|Value    |Meaning |Type |
|---------|--------|-----|
|ret_val  | any function that returns a value, will return this variable |Any |
|fl_out |variable designated for the output file a function writes to if applicable |File / String |
|fl_in |variable designated for the input file a function reads from |File / String |


#### Docstrings
In our project each method's docsring contains an:
- Overview
- Use
- Notes (optional)

Each example should also contain a description of the experiment in a seperate markdown file breifly outlining:
- Goals
- Methods
- Results
<hr>