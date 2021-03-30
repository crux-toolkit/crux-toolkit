#!/usr/bin/env python3
"""Lint the C++ files under the src/app directory

Note that to run this script, you'll need the following dependencies:
  - Python 3.6+
  - git
  - cpplint (Install with 'pip install cpplint')
"""
import os
import sys
import logging
import subprocess
from pathlib import Path
from argparse import ArgumentParser

LOGGER = logging.getLogger(__name__)


def get_args():
    """Create the parser and command line arguments.

    This is where you should add additionaly options, should they be needed.

    Returns
    -------
    ArguementParser
        The parsed command line arguments.
    """
    desc = "Lint the C++ source code of the crux mass spectrometry toolkit."
    parser = ArgumentParser(description=desc)
    parser.add_argument(
        "-d",
        "--diff-only",
        default=False,
        action="store_true",
        help=(
            "Lint only modified C++ files as found by performing a "
            "'git diff' against 'origin/master'."
        ),
    )

    parser.add_argument(
        "-l",
        "--lenient",
        default=False,
        action="store_true",
        help=(
            "Enable lenient linting. This excludes rules from"
            "'--lenient-include'"
        ),
    )
    parser.add_argument(
        "-e",
        "--exclude",
        default="exclude_rules.txt",
        type=str,
        help="A file indicating cpplint rules to exclude, one per line.",
    )

    parser.add_argument(
        "-i",
        "--include",
        default="include_lax_rules.txt",
        type=str,
        help="A file indicating cpplint rules to always include, one per line.",
    )

    parser.add_argument(
        "-n",
        "--lenient-include",
        default="include_rules.txt",
        type=str,
        help=(
            "A file indicating cpplint rules to include, one per line. "
            "These rules are excluded when '--lenient' is used."
        ),
    )

    return parser.parse_args()


def parse_rules(rule_file, include):
    """Parse the cpplint rules from a rule file.

    Parameters
    ----------
    rule_file : Path
        The file containing cpplint rules, one per line.
    include : bool
        Should the rules be included during linting?

    Returns
    -------
    str
        The parsed rules as a single string in a list.
    """
    prefix = {True: "+", False: "-"}

    with rule_file.open() as rule_ref:
        rules = rule_ref.read().splitlines()

    return ",".join([prefix[include] + r for r in rules])


def setwd_to_git_root():
    """Set the working directory to the root of the git project

    Using this makes this script modular---it should work from anywhere
    within the project.

    Returns
    -------
    Path
        The path of the project root.
    """
    cmd = ["git", "rev-parse", "--show-toplevel"]
    root = subprocess.run(cmd, capture_output=True, check=True).stdout.decode()
    root = Path(root.rstrip()).resolve()
    os.chdir(str(root))
    return root


def git_diff():
    """Perform a git diff against origin/master

    Returns
    -------
    List of Paths
        A list of the changed C++ files as Path objects.
    """
    cmd = ["git", "diff", "origin/master", "--name-only"]
    res = subprocess.run(cmd, capture_output=True, check=True).stdout.decode()
    diff_files = [Path(f) for f in sorted(res.splitlines())]
    return diff_files


def find_cpp_files(diff_files=None):
    """Find C++ files that are in the src directory

    Parameters
    ----------
    diff_files : List of Paths
        List of potential C++ files to filter. If None, this function will
        merely return a single path to the src directory.

    Returns
    -------
    List of Paths
        C++ files in the src directory to lint.
    """
    src = Path("src/app").resolve()
    if diff_files is None:
        return [src]

    # Tuple of extensions associated with C++ code.
    # Change this if you don't want to consider certain file types.
    cpp_exts = (
        ".c",
        ".c++",
        ".cc",
        ".cpp",
        ".cu",
        ".cuh",
        ".cxx",
        ".h",
        ".h++",
        ".hh",
        ".hpp",
        ".hxx",
    )

    cpp_files = []
    for diff_file in diff_files:
        is_cpp = diff_file.name.lower().endswith(cpp_exts)
        in_src = src in diff_file.resolve().parents
        if is_cpp and in_src:
            cpp_files.append(diff_file)

    return cpp_files


def cpplint(cpp_files, filters=None):
    """Run cpplint on the specified files

    Parameters
    ----------
    cpp_files: list of Paths
        The C++ files to lint.
    filters: str
        The rules to include and exclude for linting.
    """
    cmd = ["cpplint", "--recursive"]
    if filters is not None:
        cmd.append(f"--filter={filters}")

    subprocess.run(cmd + list(cpp_files), check=True, text=True)


def main():
    """The main function"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    # Parse command line arguments:
    args = get_args()

    # Parse the linting rules:
    rules = [
        parse_rules(Path(args.include), True),
        parse_rules(Path(args.lenient_include), (not args.lenient)),
        parse_rules(Path(args.exclude), False),
    ]
    rules = ",".join(rules)

    # Move to the project root:
    setwd_to_git_root()

    # Get the modified files:
    if args.diff_only:
        cpp_files = find_cpp_files(git_diff())
        LOGGER.info("Linting %i modified C++ files...", len(cpp_files))
    else:
        cpp_files = find_cpp_files()
        LOGGER.info("Linting all C++ files in src tree...")

    # If no files are modified and modified_files is not None:
    if not cpp_files:
        LOGGER.info("There was nothing to lint with this change.")
        sys.exit(0)

    # Run cpplint
    cpplint(cpp_files, filters=rules)
    LOGGER.info("Great job! Linting was completed successfully :D")


if __name__ == "__main__":
    main()
