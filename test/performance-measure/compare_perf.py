import argparse
import os
import time

def run_comparison(args, extra):
    print(args.logs_dir)
    if not os.path.isdir(args.logs_dir):
        os.makedirs(args.logs_dir)
    for t in args.threads:
        command = f"{args.stable} tide-search {" ".join(extra[1:])} --num-threads {t} --output-dir {args.output_dir} {args.data_file} {args.index_dir}"
        start_time = time.time()
        command += f" 2> {args.logs_dir}/stable_{t} 1> /dev/null"
        os.system(command)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Stable elapsed with {t} threads: {elapsed_time}")

        command = f"{args.new} tide-search {" ".join(extra[1:])} --num-threads {t} --output-dir {args.output_dir} {args.data_file} {args.index_dir}"
        start_time = time.time()
        command += f" 2> {args.logs_dir}/new_{t} 1> /dev/null"
        os.system(command)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"New version elapsed with {t} threads: {elapsed_time}")

parser = argparse.ArgumentParser(description="Run performance comparison test.")

parser.add_argument(
    "stable",
    help="stable crux executable",
    type=str,
)

parser.add_argument(
    "new",
    help="new version crux executable",
    type=str
)

parser.add_argument(
    "--threads",
    nargs='+',
    help="Number of threads to run comparison"
)

parser.add_argument(
    "--output_dir",
    type=str,
    default="./crux-output-dir"
)

parser.add_argument(
    "--data_file",
    type=str,
)

parser.add_argument(
    "--index_dir",
    type=str,
    default="./crux-index-dir"
)

parser.add_argument(
    "--fileroot",
    type=str,
    default="./crux-fileroot"
)

parser.add_argument(
    "--logs_dir",
    type=str
)

args, extra = parser.parse_known_args()
run_comparison(args, extra)
