
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

void printFinal(TString inFileName = "/home/dim/FIT/FIToutGamma/urqmd_7_fitGamma_reco16.root", TString outFileName = "/home/dim/FIT/FIToutGamma/urqmd_7reco16_fitGammaOut.tex")
{
    if (inFileName == "")
        return;

    Bool_t isTypeSimple = false; // true;
    Bool_t isTypeLatex = true;
    Bool_t isTypeScv = false;
    Bool_t isTypeCpp = false;

    // Define output file format
    if (outFileName == "")
        isTypeSimple = true;
    else if (outFileName.Contains(".tex"))
        isTypeLatex = true;
    else if (outFileName.Contains(".csv"))
        isTypeScv = true;
    else if (outFileName.Contains(".C"))
        isTypeCpp = true;
    else if (outFileName.Contains(".cpp"))
        isTypeCpp = true;
    else if (outFileName.Contains(".h"))
        isTypeCpp = true;
    else if (outFileName.Contains(".hpp"))
        isTypeCpp = true;

    TFile *fi = new TFile(inFileName.Data(), "read");
    if (!fi)
        return;

    TGraphErrors *hBavg = (TGraphErrors *)fi->Get("Fit_B_Mean");

    TTree *Result = (TTree *)fi->Get("Result");
    if (!Result)
        return;

    //Int_t Ncc;
    Float_t MinPercent;
    Float_t MaxPercent;
    Int_t MinBorder;
    Int_t MaxBorder;

    //Result->SetBranchAddress("Ncc", &Ncc);
    Result->SetBranchAddress("MinPercent", &MinPercent);
    Result->SetBranchAddress("MaxPercent", &MaxPercent);
    Result->SetBranchAddress("MinBorder", &MinBorder);
    Result->SetBranchAddress("MaxBorder", &MaxBorder);

    std::vector<std::pair<Float_t, Float_t>> vCent;
    std::vector<std::pair<Int_t, Int_t>> vBorders;
    std::vector<std::pair<Float_t, Float_t>> vBimp;

    std::vector<Float_t> vBavg;
    std::vector<Float_t> vBavgRMS;

    Int_t Nclasses = Result->GetEntries();

    // Read TTree
    for (int i = 0; i < Nclasses; i++)
    {
        Result->GetEntry(i);
        vCent.push_back({MinPercent, MaxPercent});
        vBorders.push_back({MinBorder, MaxBorder});
    }
    Double_t *vy_gr1 = hBavg->GetY();
    Double_t *ey_gr1 = hBavg->GetEY();
    // Read averaged values
    for (int i = 0; i < hBavg->GetN(); i++)
    {
        if (vy_gr1[i] != 0)
            vBavg.push_back(vy_gr1[i]);
        vBavgRMS.push_back(ey_gr1[i]);
    }

    // Fitting mean values with pol5 function
    TF1 *fitB = new TF1("fitB", "pol5", 0., 100.);

    hBavg->Fit(fitB, "R0Q");

    // Extract min/max values for B
    Float_t fmin, fmax;
    for (int i = 0; i < vBavg.size(); i++)
    {
        fmin = fitB->Eval(vCent.at(i).first);
        fmax = fitB->Eval(vCent.at(i).second);
        vBimp.push_back({fmin, fmax});
    }

    Int_t NreasonableClasses = vBavg.size();

    TTree *FitResult = (TTree *)fi->Get("FitResult");
    if (!FitResult)
        return;

    Float_t a1, a2, a3, teta, n_knee, chi2, NDF;
    Int_t minNch;

    FitResult->SetBranchAddress("teta", &teta);
    FitResult->SetBranchAddress("n_knee", &n_knee);
    FitResult->SetBranchAddress("a1", &a1);
    FitResult->SetBranchAddress("a2", &a2);
    FitResult->SetBranchAddress("a3", &a3);
    FitResult->SetBranchAddress("chi2", &chi2);
    FitResult->SetBranchAddress("NDF", &NDF);
    FitResult->SetBranchAddress("minNch", &minNch);

    FitResult->GetEntry(0);

    if (!isTypeSimple && !isTypeLatex && !isTypeScv && !isTypeCpp)
    {
        std::cerr << "Output format is not known!" << std::endl;
        return;
    }

    // Output for isTypeSimple case
    if (isTypeSimple)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        std::cout << std::endl;
        std::cout << "Cent, %   | Mult_min | Mult_max | <b>, fm |   RMS   | bmin, fm | bmax, fm |" << std::endl;
        std::cout << "----------|----------|----------|---------|---------|----------|----------|" << std::endl;
        for (int i = 0; i < NreasonableClasses; i++)
        {
            std::cout << Form("%3.0f - %3.0f | %8i | %8i | %7.2f | %7.2f | %8.2f | %8.2f |",
                              vCent.at(i).first, vCent.at(i).second,
                              vBorders.at(i).first, vBorders.at(i).second,
                              vBavg.at(i), vBavgRMS.at(i), vBimp.at(i).first, vBimp.at(i).second)
                      << std::endl;
        }
        std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
        std::cout << "  teta  | n_knee |   a1   |   a2   |   a3   |  chi2  |   NDF   | minNch |" << std::endl;
        std::cout << "--------|--------|--------|--------|--------|--------|---------|--------|" << std::endl;

        std::cout << Form(" %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %7.2f | %6i |", teta, n_knee, a1, a2, a3, chi2, NDF, minNch) << std::endl
                  << std::endl;
    }

    // Output for isTypeLatex case
    ofstream myfile;
    if (isTypeLatex)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        myfile.open(outFileName.Data());
        myfile << "\\documentclass[11pt]{article}\n";
        myfile << "\\usepackage[utf8]{inputenc}\n";
        myfile << "\\usepackage{geometry}";
        // myfile << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm}\n" << std::endl;
        myfile << "\\geometry{legalpaper, landscape, margin=2in}\n"
               << std::endl;
        myfile << "\\begin{document}\n"
               << std::endl;
        myfile << "Generated from a file:\n";
        myfile << "\\begin{verbatim*}\n";
        myfile << inFileName.Data() << std::endl;
        myfile << "\\end{verbatim*}\n";
        myfile << "\\begin{center}\n";
        myfile << "\\begin{tabular}{ |c|c|c|c|c|c|c|c| }\n";
        myfile << "\t\\hline\n";
        myfile << "\t Centrality, \\% & $N_{ch}^{min}$ & $N_{ch}^{max}$ & $\\langle b \\rangle$, fm & RMS & $b_{min}$, fm & $b_{max}$, fm  \\\\\n";
        for (int i = 0; i < NreasonableClasses; i++)
        {
            myfile << "\t\\hline\n";
            myfile << Form("\t%.0f - %.0f & %i & %i & %.2f & %.2f & %.2f & %.2f  \\\\\n",
                           vCent.at(i).first, vCent.at(i).second,
                           vBorders.at(i).first, vBorders.at(i).second,
                           vBavg.at(i), vBavgRMS.at(i), vBimp.at(i).first, vBimp.at(i).second);
        }
        myfile << "\t\\hline\n";
        myfile << "\\end{tabular}\n";
        myfile << "\\end{center}\n";

        //print fit result ------------------------

        myfile << "\\begin{center}\n";
        myfile << "\\begin{tabular}{ |c|c|c|c|c|c|c|c| }\n";
        myfile << "\t\\hline\n";

        myfile << "\t $\\theta$ & $N_{knee}$ & $a_{1}$ & $a_{2}$ & $a_{3}$ & $\\chi^{2}$ & $NDF$ & $N_{ch}^{fit}(min)$  \\\\\n";

        myfile << "\t\\hline\n";
        myfile << Form("\t%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %i \\\\\n", teta, n_knee, a1, a2, a3, chi2, NDF, minNch);
        myfile << "\t\\hline\n";
        myfile << "\\end{tabular}\n";
        myfile << "\\end{center}\n";

        //end printing fit result -----------------

        myfile << "\\end{document}";
        std::cout << "Output file " << outFileName.Data() << " is created." << std::endl;
        myfile.close();
    }
    /*
    // Outout for isTypeScv case
    if (isTypeScv)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        myfile.open(outFileName.Data());

        myfile << "Generated from a file: " << inFileName.Data() << "\n";
        myfile << "Centrality class,Mult_min,Mult_max,Mean b,b_min,b_max,Mean Npart,Npart_min,Npart_max,Mean Ncoll,Ncoll_min,Ncoll_max,\n";
        for (int i=0; i<NreasonableClasses; i++)
        {
            myfile << Form("%.0f - %.0f,%i,%i,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,\n",
                vCent.at(i).first, vCent.at(i).second,
                vBorders.at(i).first, vBorders.at(i).second,
                vBavg.at(i), vBimp.at(i).first, vBimp.at(i).second,
                vNpartavg.at(i), vNpart.at(i).second, vNpart.at(i).first,
                vNcollavg.at(i), vNcoll.at(i).second, vNcoll.at(i).first);
        }
        std::cout << "Output file " << outFileName.Data() << " is created." << std::endl;
        myfile.close();
    }

    // Output for isTypeCpp case
    if (isTypeCpp)
    {
        std::cout << "File: " << inFileName.Data() << "." << std::endl;
        myfile.open(outFileName.Data());

        myfile << "// Generated from a file: " << inFileName.Data() << "\n";
        myfile << "Float_t minCentPercent [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vCent.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vCent.at(i).first << "};\n";
        }
        myfile << "Float_t maxCentPercent [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vCent.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vCent.at(i).second << "};\n";
        }
        myfile << "Int_t minMult [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBorders.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vBorders.at(i).first << "};\n";
        }
        myfile << "Int_t maxMult [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBorders.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vBorders.at(i).second << "};\n";
        }
        myfile << "Float_t meanB [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBavg.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vBavg.at(i) << "};\n";
        }
        myfile << "Float_t minB [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBimp.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vBimp.at(i).first << "};\n";
        }
        myfile << "Float_t maxB [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vBimp.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vBimp.at(i).second << "};\n";
        }
        myfile << "Float_t meanNpart [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNpartavg.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vNpartavg.at(i) << "};\n";
        }
        myfile << "Float_t minNpart [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNpart.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vNpart.at(i).first << "};\n";
        }
        myfile << "Float_t maxNpart [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNpart.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vNpart.at(i).second << "};\n";
        }
        myfile << "Float_t meanNcoll [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNcollavg.at(i) << ", ";
            if (i==NreasonableClasses-1) myfile << vNcollavg.at(i) << "};\n";
        }
        myfile << "Float_t minNcoll [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNcoll.at(i).first << ", ";
            if (i==NreasonableClasses-1) myfile << vNcoll.at(i).first << "};\n";
        }
        myfile << "Float_t maxNcoll [" << NreasonableClasses << "] = { ";
        for (int i=0; i<NreasonableClasses; i++)
        {
            if (i!=NreasonableClasses-1) myfile << vNcoll.at(i).second << ", ";
            if (i==NreasonableClasses-1) myfile << vNcoll.at(i).second << "};\n";
        }
        

        std::cout << "Output file " << outFileName.Data() << " is created." << std::endl;
        myfile.close();
    }*/
}
