import uproot
import pandas as pd

file_path = "/project/ctb-stelzer/hamza95/photo_gen_files_PbWO4/photons_1x1x5cm_1to6GeV_PbWO4.root:photon_sim"

tree = uproot.open(file_path)

# Read selected branches into pandas
df = tree.arrays(
    ["EventID", "primaryE", "x", "y", "z", "dE", "t_ns", "material"],
    library="pd"
)

print(df.head(10))