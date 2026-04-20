import uproot

tree = uproot.open("/project/ctb-stelzer/hamza95/photo_gen_files_PbF2/photons_1x1x5cm_1to6GeV_PbF2.root:photon_sim")

count = 0

for batch in tree.iterate(["EventID"], library="np", step_size=100000):
    count += (batch["EventID"] == 2).sum()

print(count)