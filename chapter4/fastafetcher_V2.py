#!/usr/bin/env python

# Extract fasta files by their descriptors stored in a separate file.
# Requires biopython
# use -h flag for use.

# TODO:
# - Create more sophisticated logic for matching IDs/Descriptions/Partial matches etc.
#    - Create a mode variable to encapsulate invert/partial/description/id etc?
from Bio import SeqIO
import sys
import argparse

def get_keys(args):
    """Turns the input key file into a list. May be memory intensive."""
    with open(args.keyfile, "r") as kfh:
        keys = [line.rstrip("\n").lstrip(">") for line in kfh]
    return keys

def get_args():
    try:
        parser = argparse.ArgumentParser(description="Retrieve one or more fastas from a given multifasta.",)
        parser.add_argument("-f", "--fasta", action="store", required=True, help="The multifasta to search.",)
        parser.add_argument("-k", "--keyfile", action="store", help="A file of header strings to search the multifasta for. Must be one per line.",)
        parser.add_argument("-s", "--string", action="store", help="Provide a string to look for directly, instead of a file (can accept a comma separated list of strings).",)
        parser.add_argument("-o","--outfile", action="store", required=True, help="Output file to store the new fasta sequences in.",)
        parser.add_argument("-knf", "--keysnotfound", action="store", default="keys_not_found.txt", help="Output file to store the unfound header strings in. Useful for debugging.",)
        parser.add_argument("-v", "--verbose", action="store_true", help="Set whether to print the key list out before the fasta sequences. Useful for debugging.",)
        parser.add_argument("-i", "--invert", action="store_true", help="Invert the search, and retrieve all sequences NOT specified in the keyfile.",)
        parser.add_argument("-m", "--method", action="store", choices=["exact", "partial"], default="exact", help="Search the headers as exact matches, or as partial substring matches.",)
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occured with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return parser.parse_args()

def main():
    """Takes a string or list of strings in a text file (one per line) and retreives them and their sequences from a provided multifasta."""
    args = get_args()
    # Call getKeys() to create the list of keys from the provided file:
    if not (args.keyfile or args.string):
        sys.stderr.write("No key source provided. Exiting.")
        sys.exit(1)
    if args.keyfile:
        keys = get_keys(args)
    else:
        keys = args.string.split(",")

    if args.verbose:
        if args.invert is False:
            sys.stderr.write("Fetching the following keys:\n")
            for key in keys:
                sys.stderr.write(key + "\n")
        elif args.invert is True:
            sys.stderr.write(
                "Ignoring the following keys, and retreiving everything else from: {}\n".format(
                    args.fasta
                )
            )
            for key in keys:
                sys.stderr.write(key + "\n")
        sys.stderr.write(
            "-" * 80 + "\n"
        )

    # found_keys is a list of keys that have matches when parsing the multifasta
    found_keys = []
    # Parse in the multifasta and assign an iterable variable:
    to_write = []
    for rec in SeqIO.parse(args.fasta, "fasta"):
        match_found = False
        if args.invert is False:
            if args.method == "exact" and rec.id in keys:
                match_found = True
            elif args.method == "partial" and any(key in rec.description for key in keys):
                match_found = True
        else:
            if args.method == "exact" and rec.id not in keys:
                match_found = True
            elif args.method == "partial" and all(key not in rec.description for key in keys):
                match_found = True
        if match_found:
            to_write.append(rec)
            if args.verbose:
                print(rec.format("fasta"))
            if args.method == "exact":
                found_keys.append(rec.id)
            else:
                found_keys.extend([key for key in keys if key in rec.description])

    # number of found keys
    found_count = len(set(found_keys))
    # unfound_keys is a set of all not found keys when parsing the multifasta
    unfound_keys = set(keys) - set(found_keys)
    # number of unfound keys
    unfound_count = len(unfound_keys)
    
    # Ensure summary is always printed
    print(f"Number of matches found: {found_count}")
    print(f"Number of keys not found: {unfound_count}")

    if args.outfile:
        SeqIO.write(to_write, args.outfile, "fasta")

    # Write unfound keys to unfound_keys.txt if there are more than one
    if len(unfound_keys) > 1:
        print(f'Found multiple keys not found in search, writing findings to {args.keysnotfound}...')
        with open(args.keysnotfound, "w") as unfound_fh:
            for key in unfound_keys:
                unfound_fh.write(f"{key}\n")

if __name__ == "__main__":
    main()
