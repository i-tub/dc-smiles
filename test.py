"""
Test the dc-based SMILES parser using RDKit.
"""

import json
import unittest
import subprocess

from rdkit import Chem


class TestSmilesParser(unittest.TestCase):

    def check_smiles(self, smiles, useHCounts=False):
        """
        Parse a SMILES using dc, producing a JSON representation. Load the
        latter into RDKit, sanitize it, and check that the SMILES generated from
        it is the same as the one we get from a simple round-trip using RDKit.

        :param useHCounts: whether to use the hydrogen counts produced by the
                           parser or not. The parser can't guess the number of
                           hydrogens for "organic" (non-bracket) atoms, so this
                           flag should only be enabled when all the atoms in the
                           SMILES are bracketed.
        """
        cmd = f'echo "{smiles}" | rev | od -An -t d1 | ' \
            """dc -f - minified_smiles.dc"""
        output = subprocess.check_output(cmd, shell=True, encoding='utf8')

        p = Chem.JSONParseParameters()
        p.useHCounts = useHCounts

        mol = Chem.JSONToMols(output, p)[0]
        Chem.SanitizeMol(mol)
        smiles2 = Chem.MolToSmiles(mol)
        smiles1 = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        self.assertEqual(smiles1, smiles2, f'{smiles} -> {output}')


    def test_empty(self):
        self.check_smiles('')

    def test_monatomic(self):
        self.check_smiles('C')

    def test_two_letter_symbol(self):
        self.check_smiles('Cl')

    def test_linear(self):
        self.check_smiles('CCO')

    def test_bond_orders(self):
        self.check_smiles('C#CC=CC-C')

    def test_dot(self):
        self.check_smiles('CC.CC')

    def test_single_branch(self):
        self.check_smiles('CC(O)C')

    def test_double_branch(self):
        self.check_smiles('CC(F)(O)C')

    def test_nested_branch(self):
        self.check_smiles('CC(C(F)I)O')

    def test_ring(self):
        self.check_smiles('C1CC1')

    def test_kekule_benzene(self):
        self.check_smiles('C1C=CC=CC=1')

    def test_two_rings(self):
        self.check_smiles('C1CC1C1CCC1')

    def test_cubane(self):
        self.check_smiles('C12C3C4C1C5C2C3C45')

    def test_positive_charge_single(self):
        self.check_smiles('[CH3+]')

    def test_positive_charge_double(self):
        self.check_smiles('[CH2++]')

    def test_positive_charge_numeric(self):
        self.check_smiles('[CH2+2]')

    def test_negative_charge_single(self):
        self.check_smiles('[CH3-]')

    def test_negative_charge_double(self):
        self.check_smiles('[CH2--]')

    def test_negative_charge_numeric(self):
        self.check_smiles('[CH2-2]')

    def test_isotope(self):
        self.check_smiles('[13CH4]')

    def test_bracket_with_the_works(self):
        self.check_smiles('[81BrH2+]')

    def test_hcount_single(self):
        self.check_smiles('[OH]', useHCounts=True)

    def test_hcount_numeric(self):
        self.check_smiles('[NH2]', useHCounts=True)

if __name__ == '__main__':
    unittest.main()
