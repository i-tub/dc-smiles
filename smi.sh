#!/bin/bash

echo "$1" | rev | od -t d1 -An | dc -f - minified_smiles.dc
