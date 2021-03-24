#include <TH2.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>
#include <TF1.h>
#include <TLine.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/IntegratorOptions.h>
#include <TLegend.h>
#include <TPaveText.h>

const double sigma = 685, pi = TMath::Pi(), bmax = 20;
Color_t color[10] = {kRed + 2, kBlue + 1, 14, kGreen + 3, kMagenta + 3, kGreen + 1, kYellow + 2, 46, kBlue - 9, kViolet + 8};
Float_t a1 = -4.125, a2 = 1.51, a3 = -3.29, teta = 1.41, n_knee = 205.5, chi2, NDF, chi2_NDF;
int bin_cent[11];
TFile *d_outfile;

double Scale(TH1D *HistTrue, TH1D *Hist)
{
    TF1 *f1 = new TF1("f1", "gaus", 0.75 * HistTrue->FindLastBinAbove(), 0.98 * HistTrue->FindLastBinAbove());
    HistTrue->Fit(f1, "RMN");
    double MaxX1 = f1->GetParameter(1) + 3 * f1->GetParameter(2);

    TF1 *f2 = new TF1("f2", "gaus", 0.75 * Hist->FindLastBinAbove(), 0.98 * Hist->FindLastBinAbove());
    Hist->Fit(f2, "RMN");
    double MaxX2 = f2->GetParameter(1) + 3 * f2->GetParameter(2);
    cout << "MaxX2 / MaxX1 " << MaxX2 / MaxX1 << endl;

    return MaxX2 / MaxX1;
}

TGraphErrors *RatioGr(TGraphErrors *const &gr1, TGraphErrors *const &gr2, double Xmin, double Ymin, double Xmax, double Ymax)
{
    // Read points
    Double_t *vx_gr1 = gr1->GetX();
    Double_t *vy_gr1 = gr1->GetY();
    Double_t *vx_gr2 = gr2->GetX();
    Double_t *vy_gr2 = gr2->GetY();
    // Read errors
    Double_t *ex_gr1 = gr1->GetEX();
    Double_t *ey_gr1 = gr1->GetEY();
    Double_t *ex_gr2 = gr2->GetEX();
    Double_t *ey_gr2 = gr2->GetEY();

    int n1bins = gr1->GetN();
    int n2bins = gr2->GetN();
    if (n2bins < n1bins)
    {
        n1bins = n2bins;
    }
    Double_t vx_gr3[n1bins], vy_gr3[n1bins], ex_gr3[n1bins], ey_gr3[n1bins];
    for (int i = 0; i < n1bins; i++)
    {
        vx_gr3[i] = vx_gr1[i];
        ex_gr3[i] = ex_gr1[i];
        vy_gr3[i] = vy_gr1[i] / vy_gr2[i];
        ey_gr3[i] = sqrt(pow(ey_gr1[i] / vy_gr2[i], 2) + pow(vy_gr1[i] * ey_gr2[i] / (vy_gr2[i] * vy_gr2[i]), 2));
    }

    TGraphErrors *grRatio = new TGraphErrors(n1bins, vx_gr3, vy_gr3, ex_gr3, ey_gr3);
    grRatio->GetXaxis()->SetLimits(Xmin, Xmax);
    grRatio->GetYaxis()->SetRangeUser(Ymin, Ymax);

    return grRatio;
}

//Уссловная вероятность;
double PnbGamma(double *x, double *par)
{
    double Cb = x[0], n = par[0], teta = par[1], n_knee = par[2], a1 = par[3], a2 = par[4], a3 = par[5];
    double fn = n_knee * exp(a1 * Cb + a2 * pow(Cb, 2) + a3 * pow(Cb, 3)) / teta;
    return ROOT::Math::gamma_pdf(n, fn, teta);
};

double PnbGamma2(double *x, double *par)
{
    double n = x[0], Cb = par[0], teta = par[1], n_knee = par[2], a1 = par[3], a2 = par[4], a3 = par[5];
    double fn = n_knee * exp(a1 * Cb + a2 * pow(Cb, 2) + a3 * pow(Cb, 3)) / teta;
    return ROOT::Math::gamma_pdf(n, fn, teta);
};

double ftPn(double *x, double *par)
{
    // Parameters
    double teta = par[0];
    double n_knee = par[1];
    double a1 = par[2];
    double a2 = par[3];
    double a3 = par[4];

    // Variables
    double n = x[0];
    // Function
    TF1 *f = new TF1("f", PnbGamma, 0, 10, 6);
    f->SetParameters(n, teta, n_knee, a1, a2, a3);
    double func = f->Integral(0, 1);
    return func;
}

void Start(const char *fileadres, const char *current_mult, const char *outadres, int minNch, bool efficiencyFit, const char *fileadres2, const char *current_mult2)
{

    int Nch0 = minNch;
    TFile *file = new TFile(fileadres);
    TH1D *Gev = (TH1D *)file->Get(current_mult);
    Gev->Scale(1 / Gev->Integral(1, Gev->GetNbinsX(), "width"));
    bin_cent[0] = 1.05 * Gev->FindLastBinAbove();

    if (efficiencyFit == true)
    {
        int mediumNn = 0.5 * Gev->FindLastBinAbove();
        TFile *fileGl = new TFile(fileadres2);
        TH1D *Ideal = (TH1D *)fileGl->Get(current_mult2);
        Ideal->Scale(1 / Ideal->Integral(1, Ideal->GetNbinsX(), "width"));
        float Integr = Ideal->Integral(mediumNn, Ideal->GetNbinsX(), "width");
        int EfmediumNn = mediumNn * Scale(Ideal, Gev);
        Gev->Scale(Integr / Gev->Integral(EfmediumNn, Gev->GetNbinsX(), "width"));
        cout << "start3" << endl;
    }

    Gev->SetTitle("");
    Gev->GetYaxis()->SetTitle("1/N dN_{ch}/dN");
    Gev->GetXaxis()->SetTitle("N_{ch}");
    Gev->GetXaxis()->SetRangeUser(0, bin_cent[0]);
    Gev->GetYaxis()->SetRangeUser(0.5 * 1e-6, 0.5);

    // Y axis plot settings
    Gev->GetYaxis()->SetTitleSize(14);
    Gev->GetYaxis()->SetLabelSize(14);
    Gev->GetYaxis()->SetTitleFont(43);
    Gev->GetYaxis()->SetLabelFont(43);
    Gev->GetYaxis()->SetTitleOffset(1.5);
    Gev->GetYaxis()->SetLabelOffset(0.006);
    // X axis plot settings
    Gev->GetXaxis()->SetTitleSize(14);
    Gev->GetXaxis()->SetLabelSize(14);
    Gev->GetXaxis()->SetTitleFont(43);
    Gev->GetXaxis()->SetLabelFont(43);
    Gev->GetXaxis()->SetTitleOffset(3.6);
    Gev->GetXaxis()->SetLabelOffset(0.025);

    TH1F *GevC = (TH1F *)Gev->Clone();
    GevC->GetYaxis()->SetTitle("Data/Fit");
    GevC->GetYaxis()->CenterTitle(true);
    GevC->GetXaxis()->SetTitle("N_{ch}");
    GevC->GetXaxis()->SetTickLength(0.08);
    GevC->GetYaxis()->SetLabelOffset(0.015);

    //Задаем функцию для фитирования Gaus и строим верхнию половину графика
    TF1 *fc22 = new TF1("fit_func_full", ftPn, 0.1, bin_cent[0], 5);
    TF1 *fc2 = new TF1("fit_func", ftPn, Nch0, bin_cent[0], 5);
    fc2->SetParameters(teta, bin_cent[0] * 0.6, a1, a2, a3);

    TCanvas *c2 = new TCanvas("Canvas_0f_fit_result", "ratio data/fit", 650, 500);
    TPad *pad12 = new TPad("pad12", "pad12", 0, 0.3, 1, 1.0);
    pad12->SetBottomMargin(0.0);
    pad12->Draw();
    pad12->cd()->SetLogy();
    Gev->Fit(fc2, "RM");
    Gev->Draw();

    teta = fc2->GetParameter(0);
    n_knee = fc2->GetParameter(1);
    a1 = fc2->GetParameter(2);
    a2 = fc2->GetParameter(3);
    a3 = fc2->GetParameter(4);
    chi2 = fc2->GetChisquare();
    NDF = fc2->GetNDF();
    chi2_NDF = chi2 / NDF;

    fc22->SetParameters(teta, n_knee, a1, a2, a3);
    fc22->SetNpx(2000);
    fc22->Draw("SAME");
    cout << endl
         << "Fit results" << endl;
    cout << "teta=" << teta << ", n_knee=" << n_knee << ", a1=" << a1 << " ,a2=" << a2 << ", a3=" << a3 << ", chi^2/NDF=" << chi2_NDF << endl
         << endl;
    TLine *line2 = new TLine(n_knee, 0.5 * 1e-7, n_knee, 0.5);
    line2->Draw("SAME");

    TLegend *legdif = new TLegend(0.18, 0.15, 0.33, 0.35);
    legdif->SetBorderSize(0);
    legdif->AddEntry(Gev, "Data", "pl");
    legdif->AddEntry(fc2, "Fit", "pl");
    legdif->Draw();
    c2->cd(); // Go back to the main canvas before defining pad2
    TPad *pad22 = new TPad("pad22", "pad22", 0, 0, 1, 0.3);
    pad22->SetBottomMargin(0.3);
    pad22->SetTopMargin(0.0);
    pad22->Draw();
    pad22->cd();

    //TF1 *NEWfc2 = new TF1;
    //NEWfc2 = Gev->GetFunction("fit_func");
    GevC->Divide(fc22);
    GevC->SetMinimum(0.45); // Define Y ..
    GevC->SetMaximum(1.55); // .. range
    GevC->SetStats(0);      // No statistics on lower plot
    GevC->Draw();
    TLine *line = new TLine(0, 1, bin_cent[0], 1);
    line->Draw("SAME");
    d_outfile = new TFile(outadres, "recreate");
    c2->Write();
    fc2->Write();
    fc22->Write();

    TTree *treeFit = new TTree("FitResult", "FitResult");
    treeFit->Branch("teta", &teta, "teta/F");
    treeFit->Branch("n_knee", &n_knee, "n_knee/F");
    treeFit->Branch("a1", &a1, "a1/F");
    treeFit->Branch("a2", &a2, "a2/F");
    treeFit->Branch("a3", &a3, "a3/F");
    treeFit->Branch("chi2", &chi2, "chi2/F");
    treeFit->Branch("NDF", &NDF, "NDF/F");
    treeFit->Fill();

    treeFit->Write();
}
double Pb(double *b, double *Nch)
{
    // Parameters
    double n0 = Nch[0];
    double nn = Nch[1];
    // Variables
    double cb = pi * b[0] * b[0] / sigma;
    // Function
    TF1 *Pnb = new TF1("Pnb", PnbGamma2, 0, bin_cent[0], 6);
    Pnb->SetParameters(cb, teta, n_knee, a1, a2, a3);
    double IPnb = Pnb->Integral(n0, nn);
    TF1 *Pn = new TF1("Pn", ftPn, 0, bin_cent[0], 5);
    Pn->SetParameters(teta, n_knee, a1, a2, a3);
    double IPn = Pn->Integral(n0, nn);

    return 2 * pi * b[0] * IPnb / (sigma * IPn);
}
double bmean(double n0, double nn)
{
    TF1 *fPb = new TF1("fPb", Pb, 0, bmax, 2);
    fPb->SetParameters(n0, nn);
    return fPb->Mean(0, bmax);
}

double bsigma(double n0, double nn)
{
    double mean_b = 0, mean_b2 = 0;
    TF1 *fPb = new TF1("fPb", Pb, 0, bmax, 2);
    fPb->SetParameters(n0, nn);
    mean_b = fPb->Mean(0, bmax);
    mean_b2 = fPb->Moment(2, 0, bmax);
    return sqrt(mean_b2 - (mean_b * mean_b));
}

double integ(double n0, double nn)
{
    TF1 *Pn = new TF1("Pn", ftPn, 0, bin_cent[0], 5);
    Pn->SetParameters(teta, n_knee, a1, a2, a3);
    double IPn = Pn->Integral(n0, nn);
    return IPn;
}

void Rebin(double N0)
{
    double cen[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    int n0 = bin_cent[0] - 1;
    double integr = 0, norm;
    norm = integ(N0, bin_cent[0]);
    bin_cent[10]=N0;
    for (int i = 0; i < 9; i++)
    {
        while (integr < cen[i])
        {
            n0 = n0 - 1;
            integr = integ(n0, bin_cent[0]) / norm;
            
        }
        bin_cent[i + 1] = n0+1;
    }
    cout << "Centrality classes in multiplicity" << endl;
    for (int i = 1; i < 11; i++)
    {
        cout << Form("%d%s", (10 - i) * 10, "% - ") << bin_cent[i] << " ,centrality " << integ(bin_cent[i], bin_cent[0]) / norm << endl;
    }
}

void PlotMeanb()
{
    d_outfile->cd();
    Rebin(1);
    TF1 *fPb[10];
    TH1F *Himpact[10];
    double Bmean[10], Bsigma[10], centbin[10], centEr[10];
    for (int cent = 0; cent < 10; cent++)
    {
        fPb[cent] = new TF1(Form("fPb%d", cent), Pb, 0, bmax, 2);
        fPb[cent]->SetParameters(bin_cent[cent + 1], bin_cent[cent]);
        Himpact[cent] = new TH1F(Form("%s%d_%d", "B_VS_CentralityClass_", cent * 10, (cent + 1) * 10), "ImpactParametDist", 200, 0, bmax);
        Himpact[cent]->FillRandom(Form("fPb%d", cent), 200000);
        Himpact[cent]->Scale(1 / Himpact[cent]->Integral(1, Himpact[cent]->GetNbinsX(), "width"));
        Himpact[cent]->Write();

        Bmean[cent] = Himpact[cent]->GetMean();
        Bsigma[cent] = Himpact[cent]->GetRMS();
        centbin[cent] = cent * 10 + 5;
        centEr[cent] = 0;
    }

    TGraphErrors *GrFit = new TGraphErrors(10, centbin, Bmean, centEr, Bsigma);
    GrFit->SetMarkerStyle(20);
    GrFit->SetMarkerSize(1.2);
    GrFit->SetMarkerColor(1);
    GrFit->SetLineColor(1);
    GrFit->SetLineWidth(2.);
    GrFit->SetName("Fit_B_Mean");
    GrFit->SetTitle("");
    GrFit->GetXaxis()->SetTitle("Centrality, %");
    GrFit->GetYaxis()->SetTitle("<b>, fm");
    GrFit->SetTitle("B vs Centraliry");
    GrFit->SetName("Fit_B_Mean");
    GrFit->Write();

    TTree *tree = new TTree("Result", "Result");
    Float_t MinPercent;
    Float_t MaxPercent;
    Int_t MinBorder;
    Int_t MaxBorder;
    tree->Branch("MinPercent", &MinPercent, "MinPercent/F");
    tree->Branch("MaxPercent", &MaxPercent, "MaxPercent/F");
    tree->Branch("MinBorder", &MinBorder, "MinBorder/I");
    tree->Branch("MaxBorder", &MaxBorder, "MaxBorder/I");

    for (int i = 0; i < 10; i++)
    {
        MinPercent = (i) * 10;
        MaxPercent = (i+1) * 10;
        MinBorder = bin_cent[i + 1];
        MaxBorder = bin_cent[i];
        tree->Fill();
    }

    tree->Write();
    d_outfile->Close();

    cout << endl
         << "Centrality classes in impact parameter" << endl;
    cout << "Centr. class, bmin-bmax , bmean" << endl;
    for (int i = 1; i < 10; i++)
    {
        cout << Form("%d%s%d", (10 - i) * 10, "% - ", (11 - i) * 10) << "%, " << GrFit->Eval((10 - i) * 10) << " - " << GrFit->Eval((11 - i) * 10) << " fm, " << GrFit->Eval((11 - i) * 10 - 5) << " fm" << endl;
    }
    cout << Form("%d%s%d", 0, "% - ", 10) << "%, " << 0. << " - " << GrFit->Eval(10) << " fm, " << GrFit->Eval(5) << " fm" << endl;
}

void GammaFit(const char *fileadres = "/home/dim/FIT/data/UrQMD/7.7Gev/refMult_UrQMD_7.7gev_500k.root", const char *current_mult = "hRefMultSTAR", const char *outadres = "/home/dim/FIT/FIToutGamma/urqmd_7_fitGamma.root", int minNch = 20, bool efficiencyFit = false, const char *fileadres2 = "/home/dim/FIT/data/UrQMD/7.7Gev/refMult_UrQMD_7.7gev_500k.root", const char *current_mult2 = "hRefMultSTAR")
{
    Start(fileadres, current_mult, outadres, minNch, efficiencyFit, fileadres2, current_mult2);
    PlotMeanb();
}
