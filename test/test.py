import uproot

with uproot.open("pAu.root") as f:
    tree = f["nt_p_Au"]
    nu_array = tree["Ncoll"].array(library="np")

print(nu_array.mean())  # 应该 ~3.6 左右


#runAndSaveNucleons(1000, "Au", "Au", 42.0, 0.0, 0.4, 0, 0.0, 20.0, "AuAu_nucleons.root")