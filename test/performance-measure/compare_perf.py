import argparse
import os
import time
import csv

def run_comparison(args, extra):
    print(args.logs_dir)
    if not os.path.isdir(args.logs_dir):
        os.makedirs(args.logs_dir)
    
    stats = []
    for t in args.threads:
        stat = [t]
        command = f"{args.stable} tide-search {" ".join(extra[1:])} --num-threads {t} --output-dir {args.output_dir} {args.data_file} {args.index_dir}"
        start_time = time.time()
        command += f" 2> {args.logs_dir}/stable_{t} 1> /dev/null"
        os.system(command)
        end_time = time.time()
        elapsed_time = end_time - start_time
        stat.append(elapsed_time)
        print(f"Stable elapsed with {t} threads: {elapsed_time}")

        command = f"{args.new} tide-search {" ".join(extra[1:])} --num-threads {t} --output-dir {args.output_dir} {args.data_file} {args.index_dir}"
        start_time = time.time()
        command += f" 2> {args.logs_dir}/new_{t} 1> /dev/null"
        os.system(command)
        end_time = time.time()
        elapsed_time = end_time - start_time
        stat.append(elapsed_time)
        print(f"New version elapsed with {t} threads: {elapsed_time}")
        stats.append(stat)

    with open('performance.csv', 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar="|", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['threads', 'Stable', 'New'])
        for stat in stats:
            writer.writerow(stat)
        

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
