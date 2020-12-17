#!/bin/bash

# Complete List of Pathways
wget https://reactome.org/download/current/ReactomePathways.txt

# Pathways hierarchy relationship
# The "Pathway hierarchy relationship" file consists of two columns of Reactome Stable identifiers (ST_ID), defining the relationship between pathways within the pathway hierarchy. The first column provides the parent pathway stable identifier, whereas the second column provides the child pathway stable identifier.
wget https://reactome.org/download/current/ReactomePathwaysRelation.txt

# Reactome Pathways Gene Set
wget https://reactome.org/download/current/ReactomePathways.gmt.zip
unzip ReactomePathways.gmt.zip
rm -f ReactomePathways.gmt.zip

# Human reactions in SBML  (level 3, version 1) format
wget https://reactome.org/download/current/homo_sapiens.3.1.sbml.tgz

# Protein-Protein Interactions derived from Reactome pathways, human, PSI-MITAB format
# Interactions are inferred between all protein components of a complex. To limit the creation of false-positive inferences, only complexes with four or fewer protein components are considered.
# More information: https://github.com/reactome/interaction-exporter/wiki/Methods
wget https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.psi-mitab.txt

# date mark
echo "\n# Data last downloaded on:" $(date +%F) >> dnload.sh

# Data last downloaded on: 2020-12-17