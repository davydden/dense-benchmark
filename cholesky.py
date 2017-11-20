# a short script to parse Cholesky benchmark
import argparse
import re


# define command line arguments
parser = argparse.ArgumentParser(
    description='A short script to parse '
                'terminal output of Cholesky benchmark.')
parser.add_argument('files', metavar='files', nargs='+',
                    help='Files with terminal output.')
args = parser.parse_args()

# ready to start...
input_files = args.files

# we want to process the output into tables
# - for a pair of
sizes = ['5000', '20000']
#   and block
blocks = ['16', '32', '64']
# - get data in Column for solvers
#   where rows is the number of cores
#
# cores p q lapack scalapack elemental

# list of 2-tuples size-block used as key in dictionary
size_blocks = [(x, y) for x in sizes for y in blocks]

# so far empty dictionary for output
tables = {key: [] for key in size_blocks}

pattern = r'[+\-]?(?:[0-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'

for f in input_files:
    fin = open(f, 'r')
    for line in fin:
        if 'MPI' in line:
            # get number of MPI cores
            n = re.findall(pattern, line)[0]
        elif '--' not in line:
            # get actual data
            numbers = re.findall(pattern, line)
            s = numbers[0]
            b = numbers[1]
            p = numbers[2]
            q = numbers[3]
            lapack = numbers[4]
            scalapack = numbers[5]
            elemental = numbers[6]
            key = (s, b)
            if key in size_blocks:
                row = [n, p, q, lapack, scalapack, elemental]
                tables[key].append(row)

# write the data
for key in size_blocks:
    output_file = 'processed_{}_{}.txt'.format(key[0], key[1])
    print 'Writing output to {}'.format(output_file)
    fout = open(output_file, 'w')
    for row in tables[key]:
        s = ' '.join(row)
        fout.write('{}\n'.format(s))
    fout.close()

print "done!"
