import argparse

parser = argparse.ArgumentParser(description="Run the jme analysis")
parser.add_argument(
    "--inclusive-eta",
    "-i",
    action="store_true",
    help="Run over eta bins",
    default=False,
)
parser.add_argument(
    "--kill", "-k", action="store_true", help="Kill all tmux sessions", default=False
)
parser.add_argument(
    "--cartesian",
    action="store_true",
    help="Run cartesian multicuts",
    default=False,
)
parser.add_argument(
    "-p", "--parallel", action="store_true", help="Run parallel eta bins", default=False
)
parser.add_argument(
    "-s",
    "--sign",
    help="Sign of eta bins",
    type=str,
    default="neg3",
)
parser.add_argument(
    "-fs",
    "--flavsplit",
    action="store_true",
    help="Flavour split",
    default=False,
)
parser.add_argument(
    "-pnet",
    "--pnet",
    action="store_true",
    help="Use ParticleNet regression",
    default=False,
)
parser.add_argument(
    "-upart",
    "--upart",
    action="store_true",
    help="Use UparT regression",
    default=False,
)
parser.add_argument(
    "-t",
    "--test",
    action="store_true",
    help="Test run",
    default=False,
)
parser.add_argument(
    "-d",
    "--dir",
    help="Output directory",
    type=str,
    default="",
)
parser.add_argument(
    "--suffix",
    help="Suffix",
    type=str,
    default="",
)
parser.add_argument(
    "-y",
    "--year",
    help="Year",
    type=str,
    default="2023_preBPix",
)
parser.add_argument(
    "-f",
    "--flav",
    help="Flavour",
    type=str,
    default="inclusive",
)
parser.add_argument(
    "--full",
    action="store_true",
    help="Run full cartesian analysis in all eta bins and all flavours sequentially",
    default=False,
)
parser.add_argument(
    "--plot",
    action="store_true",
    help="Make plots",
    default=False,
)
parser.add_argument(
    "--central",
    action="store_true",
    help="Central eta bin (-1.3, 1.3)",
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
    "-c",
    "--closure",
    action="store_true",
    help="Produce closure test",
    default=False,
)
parser.add_argument(
    "--pnet-reg-15",
    action="store_true",
    help="Evaluate ParticleNet regression also for jet with pT < 15 GeV",
    default=False,
)
parser.add_argument(
    "--split-pnet-reg-15",
    action="store_true",
    help="Evaluate ParticleNet regression also for jet with pT < 15 GeV and slit between < and > 15 GeV",
    default=False,
)
parser.add_argument(
    "--extended-pt-bins",
    action="store_true",
    help="Evaluate ParticleNet or UparT regression for jetPT down to 8 GeV",
    default=False,
)
parser.add_argument(
    "--neutrino",
    help="Sum neutrino pT to GenJet pT",
    default=-1,
    type=int,
)
parser.add_argument(
    "--lxplus",
    action="store_true",
    help="Run on lxplus",
    default=False,
)
parser.add_argument(
    "--overwrite",
    action="store_true",
    help="Overwrite existing files",
    default=False,
)
args = parser.parse_args()


# cast the bool to int to be used as flags
args.flavsplit = int(args.flavsplit)
args.pnet = int(args.pnet)
args.upart = int(args.upart)
args.central = int(args.central)
args.closure = int(args.closure)
args.pnet_reg_15 = int(args.pnet_reg_15)
args.split_pnet_reg_15 = int(args.split_pnet_reg_15)
args.extended_pt_bins = int(args.extended_pt_bins)
args.neutrino = int(args.neutrino)
args.abs_eta_inclusive = int(args.abs_eta_inclusive)
