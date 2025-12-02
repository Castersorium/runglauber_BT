# Workflow Guidance

## Generate Data

```bash
# inside bash
root runglauber_v3.3.C

# inside root
runAndSaveNucleons(1000, "p", "p", 42.0, 0.0, 0.4, 0, 0.0, 5.0, "pp200_nucleons_1k.root")
```


## Read Data & Assign Rapidity Loss

- Manually change the input & output file name
- Manually change the centrality
- Select system

```bash
root readGlauberpA.C
```

> [!TIP]
> **Naming Scheme**
> - Data generated: `<System><Energy>_nucleons_<events>.root`
> - Rapidity assigned:`<System><Energy>_rapidityloss_<centrality>.root`

## Analyis

- Goto `JupyNotes/`

## Reference
1. [My Own Obsidian][obsidian://open?vault=Cerebrum%20extra&file=4%20-%20Workflows%2FTech%20Notes%2FTGlauberMC%20Model%20Usage]
2. https://tglaubermc.web.cern.ch/html/runglauber__v3_83_8C_source.html



# Folder
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
- ignore the `old_test`