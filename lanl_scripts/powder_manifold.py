import matplotlib as mpl
mpl.use('Qt5Agg')
import argparse
import powder_tSNE

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('csv',nargs = '*',)

    args = parser.parse_args()
    for f in args.csv:
        print(f)
        system = powder_tSNE.visualize()
