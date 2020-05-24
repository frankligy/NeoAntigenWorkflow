# NeoAntigenWorkflow

### Dependency
* python3 (Author developed on python3.7.0, but python3+ should work)
* Biopython
* requests
* xmltodict
* pathos
* numpy
* pandas
* matplotlib
```
pip install -r prerequsite.txt
```

### Module1: MHC-binding protein prediction
Normal usage:
```
python3 mhcPresent.py --help
```
If you are running on cchmc cluster:
```
bsub < sub_py.bat
```
### Module2: GTEx viewer
```
python3 queryGTEx.py --help
```

### Module3: Receptor protein prediction
```
python3 NeoEpitopePredictor.py --help
```
