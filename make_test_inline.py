
with open('minified_smiles.dc') as fh:
    dc = fh.read()
dc = dc.strip()

with open('test.py') as fh:
    test = fh.read()

test = test.replace('minified_smiles.dc', f"-e '{dc}'")
with open('test_inline.py', 'w') as fh:
    fh.write(test) 
