import json, os, argparse
class TestCaseParam:
    def __init__(self, stable, new, threads, data_file, index_dir):
        self.stable = stable
        self.new = new
        self.threads = threads
        self.data_file = data_file
        self.index_dir = index_dir

def prepare_index(case: TestCaseParam, index):
    query = f"{case.stable} tide-index "
    for key, value in index.items():
        if key == "fasta":
            continue
        query += f"--{key} {value} "
    fasta_file = index["fasta"]
    query += f" {fasta_file} {case.index_dir}"
    print(query)
    os.system(query)

parser = argparse.ArgumentParser(description="run all tests")

parser.add_argument(
    "test",
    help="test file(.json)"
)
parser.add_argument(
    "stable",
    help="stable version"
)
parser.add_argument(
    "new",
    help="new version"
)
args, _ = parser.parse_known_args()
with open(args.test, 'r') as file:
    data = json.load(file)

case = TestCaseParam(
    args.stable,
    args.new,
    data["threads"],
    data["data_file"],
    data["index_dir"]
)

if "index" in data:
        prepare_index(case, data["index"])

for test in data["tests"]:
    query = f"python3 compare_perf.py {case.stable} {case.new}"
    threads = " ".join(str(v) for v in case.threads)
    if "output_dir" in test:
        output_dir = test["output_dir"]
    else:
        output_dir = test["test_name"] + "_output"

    if "crux_logs" in test:
        crux_logs = test["crux_logs"]
    else:
        crux_logs = test["test_name"] + "_logs"
    query += f" --threads {threads} --output_dir {output_dir} --data_file {case.data_file} --index_dir {case.index_dir} --logs_dir {crux_logs}"
    query += " --"
    for key, value in test["params"].items():
        query += f" --{key} {value}"

    os.system(query)
