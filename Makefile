
minified_smiles.dc: annotated_smiles.dc
	sed 's/#.*//; s/  *//g' annotated_smiles.dc | awk '{printf $$0}' > minified_smiles.dc

test: minified_smiles.dc
	python3 test.py

test_inline.py: test.py minified_smiles.dc make_test_inline.py
	python3 make_test_inline.py

test_inline: test_inline.py
	python3 test_inline.py
	
