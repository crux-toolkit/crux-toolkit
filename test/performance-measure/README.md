# Compare the performance of two versions of Crux.
This script will systematically compare the performance of two versions of crux to each other, ideally corresponding to a stable version and an updated version of the program. The script runs a series of tests specified in a given json file.
```bash
python3 run-tests.py [TEST_JSON] [PATH_TO_CRUX_STABLE] [PATH_TO_CRUX_NEW]
```
`run-tests.py` would run all the tests in json file. The parameters for each test you could also modify in json file.
## JSON schema
Each `.json` file represents test cases of tide-search on the same tide index.
### Index config
The index config should be written to test config as a dict(params - value). Also, there should be fasta file parameter as `"fasta"`

```json
{
    ...
    "index": {
        "fasta": "/path/to/fasta", // required
        ... // parameters
    }
    ...
}
```

## Test cases
All the test cases should be placed into `"tests"` parameter. Each test case represents two tide-search runs(stable and new). For each test case there should be written tide-search parameters in `"parameters"`. Also, for each test there should be `"test_name"` with the name of test case.
```json
{
    ...
    "tests": [
        { // one test case
            "test_name": "some name",
            "parameters": {
                ... // tide-search parameters
            }
        }
    ]
}
```


