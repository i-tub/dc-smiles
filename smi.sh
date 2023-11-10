#!/bin/bash

echo "$1" | rev | od -t d1 -An | dc -f - smi.dc
