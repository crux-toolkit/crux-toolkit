import json, os, argparse
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
for test in data["tests"]:
    query = f"python3 compare_perf.py {args.stable} {args.new}"
    t = test["threads"]
    threads = " ".join(str(v) for v in t)
    query += f" --threads {threads} --output_dir {test["output_dir"]} --data_file {test["data_file"]} --index_dir {test["index_dir"]}"
    query += " --"
    for key, value in test["params"].items():
        query += f" --{key} {value}"

    query += f" 2> {test["crux-logs"]}"
    os.system(query)