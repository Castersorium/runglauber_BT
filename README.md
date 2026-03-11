# Workflow Guidance

## Generate Data

```bash
# inside bash
root runglauber_v3.3.C

# inside root

## p+p
runAndSaveNucleons(1000, "p", "p", 42.0, 0.0, 0.4, 0, 0.0, 5.0, "nucleonGeneratedData/pp200_nucleons_1k.root")

## A+A
runAndSaveNucleons(100000, "Pbpnrw", "Pbpnrw", 31.70, 0.0, 0.4, 0, 0.0, 20.0, "PbpnrwPbpnrw_17p3_nucleons_10k.root")
runAndSaveNucleons(100000, "Au197pnHFB14", "Au197pnHFB14", 35.75, 0.0, 0.4, 0, 0.0, 20.0, "AuAu197pnHFB14_62p4_nucleons_10k.root")
runAndSaveNucleons(100000, "Au197pnHFB14", "Au197pnHFB14", 41.29, 0.0, 0.4, 0, 0.0, 20.0, "AuAu197pnHFB14_200_nucleons_10k.root")

```

### Parameters
1. Event Number: Adjust if needed.
2. Target: "p", "Au197pnHFB14", "Pbpnrw". [3]
3. Projectile: "p", "Au197pnHFB14", "Pbpnrw". [3]
4. Inelastic cross section, usually determined by energy: 41.29 for 200 GeV; 35.75 for 62.4 GeV; 31.70 for 17.3 GeV.
    - Note: Could use `getSigmaNN()`; for below 10 GeV, use  `_Bystricky()` [4]
5. cross section width
6. Minimum distance between nucleons (fm)
7. verbose
8. b min
9. b max
10. Output file name

## Read Data & Assign Rapidity Loss

- Manually change the input & output file name
- Manually change the centrality
- Select system

```bash
./saveMultipleTree.sh
# parameters can also be tuned inside bash script
```
### saveTreeAA.C 
[My suggestion is to change it]
Apply the rapidity loss and save in a ROOT Tree.
The data is now stored in `rapidityTree/`

### readTree.C
QA plot generation. Stored in `rapidity_dNdy/`

> [!TIP]
> **Naming Scheme**
> - Data generated: `<System>_<Energy>_nucleons_<events>.root`
> - Rapidity assigned:`<System>_<Energy>_rapidityloss_<centrality>.root`

## Analyis

- Experiment data in `include/Exp_dNdy.h`

- Goto `cumulant/`

## Reference
1. [My Own Obsidian][obsidian://open?vault=Cerebrum%20extra&file=4%20-%20Workflows%2FTech%20Notes%2FTGlauberMC%20Model%20Usage]
2. https://tglaubermc.web.cern.ch/html/runglauber__v3_83_8C_source.html
3. https://tglaubermc.web.cern.ch/html/classTGlauNucleus.html#a7f8ad270b7cdc1b0f9fe811ad465e5b1
4. https://tglaubermc.web.cern.ch/html/runglauber__v3_83_8C.html#a5ef700860747f3a3dedc50765519c677

<!-- # Folder
```
├── BT
│   ├── 3Dplot/
│   ├── JupyNotes/
│   ├── dat/
│   ├── old_test/
└── README.md
```
- `3Dplot` just draws some 3D picture of glauber event
- `JupyNotes` is for data analysis
- `dat` is for TGlauberMC
- ignore the `old_test` -->