dc-based SMILES parser
======================

I have come to realize that using RDKit to parse SMILES has been a
distraction, because as it turns out, UNIX has shipped with a built-in SMILES
parser since the 1970s! It is a special mode of the `dc` program; one only
needs to invoke dc with the right one-word command to enable SMILES-parsing
mode.

Given the legendary stability of UNIX, it only makes sense to switch back to
using the 50-year old SMILES parser that comes with dc instead of relative
newcomers such as RDKit.

Fun historical tidbit: according to Doug McIlroy, dc is the "senior language
on UNIX systems":

> The PDP-11 assembler, a desk calculator dc, and B itself were written in B
> to bootstrap the system to the PDP-11. Because it could run before the disk
> had arrived, dc-—not the assembler—-became the first language to run on our
> PDP-11. Soon revised to handle arbitrary precision numbers (v1), dc was
> taken over by Bob Morris and Lorinda Cherry. It now ranks as the senior
> language on UNIX systems.

https://www.cs.dartmouth.edu/~doug/reader.pdf

Other that stripping of whitespace, the SMILES-parsing dc command is not
deliberately obfuscated (choice of language notwithstanding), and in fact it
is as structured and readable as I could make it, even sacrificing brevity.
For the interested reader, I have provided an annotated version, which
includes whitespace and extensive comments.

Requirements
------------

Tested with:

- dc (GNU bc 1.06.95) 1.3.95.
- Python 3.8.10
- RDKit 2023.09.1
- od (GNU coreutils) 8.22
- rev from util-linux 2.23.2.

Contents
--------

- Makefile: for generating the minified version and running tests
- annotated_smiles.dc: the human "readable" version
- minified_smiles.dc: same, after stripping out comments and whitespace
- make_test_inline.py: generate a version of test.py with the dc code inline
- test.py: run tests using minified_smi.dc
- test_inline.py: run tests with the dc code inlined
