# Quick information with step-by-step instruction.

Below is a detailed description of the program used to reconstruct the distribution of the impact parameter.

## 1. Fitting model data

In the case of fitting model data, it is sufficient to have only a histogram with a multiplicity distribution.

To start the fitting process, dowland script from `https://github.com/Dim23/GammaFit.git` and run `GammaFit.C` :

    root -b GammaFit.C(fileadres, current_mult, outadres, minNch)

Where the arguments are:

`fileadres` - input root file.

`current_mult` - Histogram of multiplisity.

`outadres` - output file from fitting process.

`minNch` - lower value of the fitting area ( `20` by default).

## 1.1 OUTPUT

Resulting file `outadres` will contain TCanvas with fit results and data-to-fit ratio - `Canvas_0f_fit_result` .

The resulting fit function of the multiplicity distribution - `fit_func` .

TGraphErrors of impact parametr as a function of centrality - `fit_B_Mean` .

Histograms with the distribution of the impact parameter in centrality class - `ImpactParametDist_CENT*_*` .

## 2. Fitting reconstracted data

When fitting the reconstructed data, it is necessary to take into account the efficiency of the detectors. To take these features into account, normalization is made to non-reconstructed model data.

To start the fitting process, run `GammaFit.C` with next option:

    root -b GammaFit.C(fileadres, current_mult, outadres, minNch, efficiencyFit,fileadres2, current_mult2)

Where the arguments are:

`fileadres` - input root file with multiplisity of the reconstructed data.

`current_mult` - Histogram with multiplisity.

`outadres` - output file from fitting process.

`minNch` - lower value of the fitting area ( `20` by default).

`efficiencyFit` - should be set `true` ( `false` by default).

`fileadres2` - input root file with multiplisity of non-reconstructed model data.

`current_mult2` - Histogram with multiplisity of non-reconstructed model data.


## 2.1 OUTPUT

Resulting file `outadres` will contain TCanvas with fit results and data-to-fit ratio - `Canvas_0f_fit_result` .

The resulting fit function of the multiplicity distribution - `fit_func` .

TGraphErrors of impact parametr as a function of centrality - `fit_B_Mean` .

Histograms with the distribution of the impact parameter in centrality class - `ImpactParametDist_CENT*_*` .
