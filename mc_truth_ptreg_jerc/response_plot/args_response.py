import argparse

parser = argparse.ArgumentParser(description="Run the jme analysis")
parser.add_argument(
    "-u",
    "--unbinned",
    help="Binned or unbinned",
    action="store_true",
    default=False,
)
parser.add_argument(
    "-l",
    "--load",
    action="store_true",
    help="Load medians from file",
    default=False,
)
parser.add_argument(
    "-d",
    "--dir",
    type=str,
    help="Input dir",
)
parser.add_argument(
    "--cartesian",
    action="store_true",
    help="Run cartesian multicuts",
    default=False,
)
parser.add_argument(
    "--histo",
    action="store_true",
    help="Plot the histograms",
    default=False,
)
parser.add_argument(
    "--full",
    action="store_true",
    help="Run full cartesian analysis in all eta bins and all flavours sequentially",
    default=False,
)
parser.add_argument(
    "--test",
    action="store_true",
    help="Run test",
    default=False,
)
parser.add_argument(
    "--central",
    action="store_true",
    help="Run central eta bin",
    default=False,
)
parser.add_argument(
    "-a",
    "--abs-eta-inclusive",
    action="store_true",
    help="Run over inclusive abs eta bins",
    default=False,
)
parser.add_argument(
    "-n",
    "--num-processes",
    type=int,
    help="Number of processes",
    default=16,
)
parser.add_argument(
    "-p",
    "--num-params",
    type=int,
    help="Num param fit polynomial + 2 for the jet pt range. This is forced to 13 for extendedPT for uparT and PNET",
    default=9,
)
parser.add_argument(
    "--no-plot",
    action="store_true",
    help="Do not plot",
    default=False,
)
parser.add_argument(
    "-c",
    "--choose-plot",
    action="store_true",
    help="Choose the plot",
    default=False,
)
parser.add_argument(
    "--flav",
    help="Flavour",
    type=str,
    default="inclusive",
)
parser.add_argument(
    "-y",
    "--year",
    help="Year",
    type=str,
    default="",
)
parser.add_argument(
    "--all-flavs",
    help="Do all flavours",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--upart",
    help="Use UparT Regression instead of PNet one",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--ptmin-fit-inv-median",
    help="Minimum pT for fit in inv median",
    type=float,
    default=8.0,       #17.0 for PNET, 15.0 for UparT ?
)

args = parser.parse_args()
