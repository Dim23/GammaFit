# Quick information with step-by-step instruction.


To start the fitting process, dowland script from `https://github.com/Dim23/GammaFit.git` and run `GammaFit.C` :

    root -b GammaFit.C(fileadres, current_mult, outadres, minNch)

Where the arguments are:

`fileadres` - input root file.

`current_mult` - TH1D Histogram with multiplisity.

`outadres` - output file from fitting process.

`minNch` - lower value of the fitting area ( `20` by default).

### OUTPUT

Resulting file contains TCanvas showing fit results and data-to-fit ratio - `Canvas_0f_fit_result` .

The resulting fit function of the multiplicity distribution - `fit_func` .

TGraphErrors of impact parametr as a function of centrality - `fit_B_Mean` .

Histograms with the distribution of the impact parameter in centrality class - `ImpactParametDist_CENT*_*`.
