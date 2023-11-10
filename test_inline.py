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
            """dc -f - -e '[0Sbls0:bln1:blo2:blm1+sm]sB[0sylYxlyln:a_1ln:hln1+snls_1!=Blnss1so]sO[lOxq]sA[szlsSsq]s([szLsszq]s)[sz1soq]s-[sz2soq]s=[sz3soq]s%[sz_1ssq]s.[0Sb0:bln1:blo2:blm1+sm0r:r3Q]sR[d;rd0!=Rszlnr:rq]sD[48-lx10*+sxd58r<I]sI[q]sV[ly256*+syz0=Vd96r>Y]sY[szd64=@]s@[sz1sh[[48-sh]SHd47r>HLHsz]SHd58r<HLHsz]sH[szlq1-sqd45=W[48-_1*sq]SWd58r<WLWsz]sW[szlq1+sqd43=X[48-sq]SXd58r<XLXsz]sX[sz0Sxd58r<I0sylYxd64=@0Shd72=H0Sqd45=Wd43=XszlylOxln1-ddLxr:xLqr:qLhr:hq]sK[d35=%d40=(d41=)d45=-d46=.d47=-d92=-d61==d91=Kd64r>Ad58r<Dse]sC[lCxz0!=L]sL[44P32P]sQ[44P10P]sN[32P32P["error":]P32P[true,]P10P]sE[rszq]s>[d1r72=>szd5r66=>szd6r67=>szd7r78=>szd8r79=>szd9r70=>szd15r80=>szd16r83=>szd17r17260=>szd35r17010=>szd53r73=>szsz0]sF[lQx["impHs":]P32Pli;hn]sG[32P32P32P123P["z":]P32Pli;alFxnlQx["isotope":]P32Pli;xnlQx["chg":]P32Pli;qnli;h_1!=Gli1+ddsi125Pln!=Nln!=T]sT[32P32P32P123P["atoms":]P32P91P0;b1-nlQx1;b1-n93PlQx["bo":]P32P2;bn125PLbszli1+ddsilm!=Nlm!=U]sU[123P10Ple_1!=E32P["commonchem":]P32P123P["version":]P32P10n125P44P10P32P["defaults":]P32P123P10P32P32P["atom":{"stereo":"unspecified"},]P10P32P32P["bond":{"stereo":"unspecified"}]P10P32P125P44P10P32P["molecules":]P32P91P123P10P32P32P["atoms":]P32P91P10P0siln0!=T10P32P32P93PlNx32P32P["bonds":]P32P91P10P0silm0!=U10P32P32P93P10P32P125P93P10P125P10P]sJ[sz]sZ[0sn0sm1so_1ss_1sez0!=Zz0!=LlJx]sMlMx'"""
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
