import unittest
from expand_dataset import usable, structure

class ExpansionTests(unittest.TestCase):
    def test_no_censored_invalid_or_other_endpoint(self):
        a=dict(standard_relation='=',standard_type='IC50',standard_units='nM',data_validity_comment=None,potential_duplicate=False)
        self.assertTrue(usable(a))
        for field,value in [('standard_relation','>'),('standard_type','ED50'),('standard_units','%'),('potential_duplicate',1),('data_validity_comment','Outside typical range')]:
            b=dict(a);b[field]=value;self.assertFalse(usable(b))
    def test_salts_and_stereochemistry_excluded(self):
        self.assertEqual(structure('CC(=O)[O-].[Na+]')[1],structure('CC(=O)O')[1])
        self.assertEqual(structure('C[C@H](O)C(=O)O')[1],structure('C[C@@H](O)C(=O)O')[1])
if __name__=='__main__':unittest.main()
