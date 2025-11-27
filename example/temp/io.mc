##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue Nov 25 23:34:25 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path 2QCS_A
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (161, 161, 161)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 5527 atoms
Valist_getStatistics:  Max atom coordinate:  (102.672, 2.14, 13.113)
Valist_getStatistics:  Min atom coordinate:  (46.759, -54.531, -44.622)
Valist_getStatistics:  Molecule center:  (74.7155, -26.1955, -15.7545)
NOsh_setupCalcMGAUTO(/home/runner/work/apbs/apbs/src/generic/nosh.c, 1868):  coarse grid center = 74.7155 -26.1955 -15.7545
NOsh_setupCalcMGAUTO(/home/runner/work/apbs/apbs/src/generic/nosh.c, 1873):  fine grid center = 74.7155 -26.1955 -15.7545
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1885):  Coarse grid spacing = 0.619257, 0.632358, 0.638711
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1887):  Fine grid spacing = 0.489269, 0.496975, 0.500713
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1889):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.79009, 0.785908, 0.783942 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1983):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1985):  coarse mesh center = 74.7155 -26.1955 -15.7545
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1990):  coarse mesh upper corner = 124.256 24.3931 35.3424
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1995):  coarse mesh lower corner = 25.1749 -76.7841 -66.8514
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 2000):  initial fine mesh upper corner = 113.857 13.5625 24.3025
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 2005):  initial fine mesh lower corner = 35.574 -65.9535 -55.8115
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 2066):  final fine mesh upper corner = 113.857 13.5625 24.3025
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 2071):  final fine mesh lower corner = 35.574 -65.9535 -55.8115
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 37.6995
Vpbe_ctor2:  solute dimensions = 58.283 x 59.516 x 60.114
Vpbe_ctor2:  solute charge = 4
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (66.989, 67.747, 68.811)
Vclist_setupGrid:  Grid lower corner = (41.221, -60.069, -50.16)
Vclist_assignAtoms:  Have 3631970 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.304374
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 8.832490e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (161, 161, 161)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.076160e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (081, 081, 081)
Vbuildops: Galer: (041, 041, 041)
Vbuildops: Galer: (021, 021, 021)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 4.164190e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.491554e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.119031e-01
Vprtstp: contraction number = 1.119031e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.433248e-02
Vprtstp: contraction number = 1.280794e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.019890e-03
Vprtstp: contraction number = 1.409309e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 3.040137e-04
Vprtstp: contraction number = 1.505101e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 4.842455e-05
Vprtstp: contraction number = 1.592841e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 8.256344e-06
Vprtstp: contraction number = 1.704991e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.648695e-06
Vprtstp: contraction number = 1.996883e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 4.255403e-07
Vprtstp: contraction number = 2.581073e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.568549e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 2.154527e+00
Vpmg_setPart:  lower corner = (25.1749, -76.7841, -66.8514)
Vpmg_setPart:  upper corner = (124.256, 24.3931, 35.3424)
Vpmg_setPart:  actual minima = (25.1749, -76.7841, -66.8514)
Vpmg_setPart:  actual maxima = (124.256, 24.3931, 35.3424)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 1.212320139207E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 5.696000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 37.6995
Vpbe_ctor2:  solute dimensions = 58.283 x 59.516 x 60.114
Vpbe_ctor2:  solute charge = 4
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (66.989, 67.747, 68.811)
Vclist_setupGrid:  Grid lower corner = (41.221, -60.069, -50.16)
Vclist_assignAtoms:  Have 3631970 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = 35.574, -65.9535, -55.8115
VPMG::focusFillBound -- New mesh maxs = 113.857, 13.5625, 24.3025
VPMG::focusFillBound -- Old mesh mins = 25.1749, -76.7841, -66.8514
VPMG::focusFillBound -- Old mesh maxs = 124.256, 24.3931, 35.3424
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (35.574, -65.9535, -55.8115)
Vpmg_setPart:  upper corner = (113.857, 13.5625, 24.3025)
Vpmg_setPart:  actual minima = (25.1749, -76.7841, -66.8514)
Vpmg_setPart:  actual maxima = (124.256, 24.3931, 35.3424)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (35.574, -65.9535, -55.8115)
VPMG::extEnergy    Disj part upper corner = (113.857, 13.5625, 24.3025)
VPMG::extEnergy    Old lower corner = (25.1749, -76.7841, -66.8514)
VPMG::extEnergy    Old upper corner = (124.256, 24.3931, 35.3424)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 1.2274 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.302746
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.154494e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (161, 161, 161)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.108270e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (081, 081, 081)
Vbuildops: Galer: (041, 041, 041)
Vbuildops: Galer: (021, 021, 021)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 4.154800e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 4.824117e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.350556e-01
Vprtstp: contraction number = 1.350556e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.788746e-02
Vprtstp: contraction number = 1.324451e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.592651e-03
Vprtstp: contraction number = 1.449424e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 4.013806e-04
Vprtstp: contraction number = 1.548148e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 6.528868e-05
Vprtstp: contraction number = 1.626603e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 1.146561e-05
Vprtstp: contraction number = 1.756141e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 2.364482e-06
Vprtstp: contraction number = 2.062238e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 6.034136e-07
Vprtstp: contraction number = 2.551990e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.561826e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 2.148707e+00
Vpmg_setPart:  lower corner = (35.574, -65.9535, -55.8115)
Vpmg_setPart:  upper corner = (113.857, 13.5625, 24.3025)
Vpmg_setPart:  actual minima = (35.574, -65.9535, -55.8115)
Vpmg_setPart:  actual maxima = (113.857, 13.5625, 24.3025)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 1.896056854674E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 5.678000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
printEnergy:  Performing global reduction (sum)
Vcom_reduce:  Not compiled with MPI, doing simple copy.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 7.383600e+00
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue Nov 25 23:34:33 2025

##############################################################################
Vgrid_readDX:  Grid dimensions 161 x 161 x 161 grid
Vgrid_readDX:  Grid origin = (35.574, -65.9535, -55.8115)
Vgrid_readDX:  Grid spacings = (0.489269, 0.496975, 0.500713)
Vgrid_readDX:  allocating 161 x 161 x 161 doubles for storage
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue Nov 25 23:38:52 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path RECP_A
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (161, 161, 161)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 4560 atoms
Valist_getStatistics:  Max atom coordinate:  (72.297, 54.113, 62.889)
Valist_getStatistics:  Min atom coordinate:  (6.935, -1.346, 2.25)
Valist_getStatistics:  Molecule center:  (39.616, 26.3835, 32.5695)
NOsh_setupCalcMGAUTO(/home/runner/work/apbs/apbs/src/generic/nosh.c, 1868):  coarse grid center = 39.616 26.3835 32.5695
NOsh_setupCalcMGAUTO(/home/runner/work/apbs/apbs/src/generic/nosh.c, 1873):  fine grid center = 39.616 26.3835 32.5695
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1885):  Coarse grid spacing = 0.719971, 0.619002, 0.674581
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1887):  Fine grid spacing = 0.548512, 0.489119, 0.521813
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1889):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.761853, 0.790173, 0.773535 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1983):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1985):  coarse mesh center = 39.616 26.3835 32.5695
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1990):  coarse mesh upper corner = 97.2137 75.9036 86.536
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 1995):  coarse mesh lower corner = -17.9817 -23.1367 -21.397
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 2000):  initial fine mesh upper corner = 83.497 65.513 74.3145
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 2005):  initial fine mesh lower corner = -4.265 -12.746 -9.1755
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 2066):  final fine mesh upper corner = 83.497 65.513 74.3145
NOsh_setupCalcMGAUTO (/home/runner/work/apbs/apbs/src/generic/nosh.c, 2071):  final fine mesh lower corner = -4.265 -12.746 -9.1755
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 41.3535
Vpbe_ctor2:  solute dimensions = 67.762 x 58.259 x 63.49
Vpbe_ctor2:  solute charge = 1.63203e-14
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (76.438, 66.535, 71.715)
Vclist_setupGrid:  Grid lower corner = (1.397, -6.884, -3.288)
Vclist_assignAtoms:  Have 2626625 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.253046
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 7.866680e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (161, 161, 161)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.074590e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (081, 081, 081)
Vbuildops: Galer: (041, 041, 041)
Vbuildops: Galer: (021, 021, 021)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 4.136580e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.386465e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.017335e-01
Vprtstp: contraction number = 1.017335e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.379920e-02
Vprtstp: contraction number = 1.356406e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.123136e-03
Vprtstp: contraction number = 1.538594e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 3.558071e-04
Vprtstp: contraction number = 1.675857e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 6.365553e-05
Vprtstp: contraction number = 1.789046e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 1.194353e-05
Vprtstp: contraction number = 1.876275e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 2.322532e-06
Vprtstp: contraction number = 1.944594e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 4.659896e-07
Vprtstp: contraction number = 2.006387e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.584769e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 2.167219e+00
Vpmg_setPart:  lower corner = (-17.9817, -23.1367, -21.397)
Vpmg_setPart:  upper corner = (97.2137, 75.9036, 86.536)
Vpmg_setPart:  actual minima = (-17.9817, -23.1367, -21.397)
Vpmg_setPart:  actual maxima = (97.2137, 75.9036, 86.536)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 8.954660786116E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 5.734000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 41.3535
Vpbe_ctor2:  solute dimensions = 67.762 x 58.259 x 63.49
Vpbe_ctor2:  solute charge = 1.63203e-14
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (76.438, 66.535, 71.715)
Vclist_setupGrid:  Grid lower corner = (1.397, -6.884, -3.288)
Vclist_assignAtoms:  Have 2626625 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -4.265, -12.746, -9.1755
VPMG::focusFillBound -- New mesh maxs = 83.497, 65.513, 74.3145
VPMG::focusFillBound -- Old mesh mins = -17.9817, -23.1367, -21.397
VPMG::focusFillBound -- Old mesh maxs = 97.2137, 75.9036, 86.536
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-4.265, -12.746, -9.1755)
Vpmg_setPart:  upper corner = (83.497, 65.513, 74.3145)
Vpmg_setPart:  actual minima = (-17.9817, -23.1367, -21.397)
Vpmg_setPart:  actual maxima = (97.2137, 75.9036, 86.536)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-4.265, -12.746, -9.1755)
VPMG::extEnergy    Disj part upper corner = (83.497, 65.513, 74.3145)
VPMG::extEnergy    Old lower corner = (-17.9817, -23.1367, -21.397)
VPMG::extEnergy    Old upper corner = (97.2137, 75.9036, 86.536)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 0.325524 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.250849
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.024008e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (161, 161, 161)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.067770e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (081, 081, 081)
Vbuildops: Galer: (041, 041, 041)
Vbuildops: Galer: (021, 021, 021)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 4.198790e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 4.605141e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.334860e-01
Vprtstp: contraction number = 1.334860e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.839602e-02
Vprtstp: contraction number = 1.378124e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.836942e-03
Vprtstp: contraction number = 1.542150e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 4.750981e-04
Vprtstp: contraction number = 1.674684e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 8.516721e-05
Vprtstp: contraction number = 1.792624e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 1.614303e-05
Vprtstp: contraction number = 1.895451e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 3.204703e-06
Vprtstp: contraction number = 1.985193e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 6.605342e-07
Vprtstp: contraction number = 2.061140e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.578441e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 2.165697e+00
Vpmg_setPart:  lower corner = (-4.265, -12.746, -9.1755)
Vpmg_setPart:  upper corner = (83.497, 65.513, 74.3145)
Vpmg_setPart:  actual minima = (-4.265, -12.746, -9.1755)
Vpmg_setPart:  actual maxima = (83.497, 65.513, 74.3145)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 1.486605963327E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 5.776000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
printEnergy:  Performing global reduction (sum)
Vcom_reduce:  Not compiled with MPI, doing simple copy.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 7.237532e+00
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue Nov 25 23:39:00 2025

##############################################################################
Vgrid_readDX:  Grid dimensions 161 x 161 x 161 grid
Vgrid_readDX:  Grid origin = (-4.265, -12.746, -9.1755)
Vgrid_readDX:  Grid spacings = (0.548512, 0.489119, 0.521813)
Vgrid_readDX:  allocating 161 x 161 x 161 doubles for storage
