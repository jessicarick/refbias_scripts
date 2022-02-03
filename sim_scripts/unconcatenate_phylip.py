#!/usr/bin/env python3

# By Jonathan Chang - September 2017

import sys, argparse, logging, re, os


def get_partitions(path):
    """Given a path, returns a dictionary with partition names as keys and a list of
    length 2, representing a range of characters."""
    with open(path) as handle:
        res = {}
        for line in handle:
            line = line.strip()
            if line:
                parts = re.split("\s*[,=\-]\s*", line)
                assert len(parts) == 4, "part needs [type, name, start, end], got %r" % parts
                res["_".join(parts[0:2])] = [int(x) for x in parts[2:]]
        return res


def parse_phylip(path):
    """Given a path, tries to parse that file as a phylip format character matrix. Does some
    basic error checking."""
    with open(path) as handle:
        for no, line in enumerate(handle, start=1):
            line = line.strip()
            if no == 1:
                ntax, nchar = [int(x) for x in line.split()]
                continue
            name, seq = line.split(None, 1)
            if len(seq) != nchar:
                raise Exception("{} wants {} characters, got {}".format(name, nchar, len(seq)))
            yield name, seq
        if (no - 1 != ntax):
            raise Exception("{} wants {} taxa, got {}".format(path, no-1, ntax))

def get_phylip_ntax(path):
    with open(path) as rfile:
        line = rfile.readline()
        ntax, _ = line.split()
    return int(ntax)


def append_sequence(path, name, sequence, kind="phylip"):
    """Writes a given sequence to a sequence file specified by `path`."""
    if kind == "phylip":
        fmt = "{} {}\n"
    else:
        fmt = ">{}\n{}\n"
    with open(path, "a") as wfile:
        wfile.write(fmt.format(name, sequence))

def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    partitions = get_partitions(args.partition)
    if args.prefix:
        base = args.prefix
        fn_format = base + "{}.{}"
    else:
        base = os.path.splitext(args.matrix)[0]
        fn_format = base + "_{}.{}"

    ntax = get_phylip_ntax(args.matrix)
    logging.info("Working on {} taxa".format(ntax))

    for no, (identifier, sequence) in enumerate(parse_phylip(args.matrix), start=1):
        for part_name, [start, end] in partitions.items():
            filename = fn_format.format(part_name, args.type)
            chars = sequence[slice(start-1, end)]
            if args.trim:
                chars = chars.strip("-")
            if args.remove_missing:
                tmp = chars.strip("-").strip("?")
                if not tmp:
                    chars = tmp
            if chars:
                if no == 1 and args.type == "phylip":
                    with open(filename, "w") as wfile:
                        wfile.write("{} {}\n".format(ntax, len(chars)))

                append_sequence(filename, identifier, chars, kind=args.type)

        logging.debug("{} {}".format(no, identifier))

    for part_name in partitions.keys():
        logging.info("Wrote to " + fn_format.format(part_name, args.type))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description = "Unconcatenates a raxml-style phylip/partitions file",
            epilog = "Arguments can also be placed in a file, one per line, and specified like '%(prog)s @params.conf'.",
            fromfile_prefix_chars = '@' )

    parser.add_argument("matrix", help="matrix to unconcatenate")
    parser.add_argument("partition", help="raxml-formatted partitions file")
    parser.add_argument("-t", "--type", help="output type", choices=["fasta", "phylip"], default="phylip")
    parser.add_argument("-P", "--prefix", help="change default output prefix")
    parser.add_argument("--verbose", "-v", action="store_true", help="make lots of noise while working")
    parser.add_argument("--trim", action="store_true", help="trim gaps and remove sequences with only gaps")
    parser.add_argument("--remove-missing", action="store_true", help="remove taxa that are only missing")

    args = parser.parse_args()

    if args.remove_missing and args.type == "phylip":
        logging.warn("Setting --remove-missing and --type=phylip will produce incorrect species counts.")

    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    main(args, loglevel)
