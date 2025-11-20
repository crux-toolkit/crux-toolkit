import argparse
import os
import time

def perf_record(args, extra):
    command = f"perf record -F {args.hertz} -g -a {args.file} {" ".join(extra[1:])} > {args.perf_output}"
    print(f"Run command: {command}")
    os.system(command)

    draw_flamegraph = f""
    

parser = argparse.ArgumentParser(description="Run perf test.")

parser.add_argument(
    "file",
    help="crux executable",
    type=str,
)

parser.add_argument(
    "--output",
    help="output flamegraph file",
    type=str,
    default="output.svg",
)

parser.add_argument(
    "--hertz",
    default=99,
    type=int,
    help="hetrz number in perf record command"
)

parser.add_argument(
    "--perf_output",
    default="out.perf_stack",
    nargs=1,
    type=str,
    help="perf record output stack file"
)

args, extra = parser.parse_known_args()
print(args.file)
print(extra)
perf_record(args, extra)
