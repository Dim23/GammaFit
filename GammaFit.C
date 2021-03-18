#include <TH2.h>
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

//Задаем постоянные;
const double sigma = 685, pi = TMath::Pi(), bmax = 20;
//int color[10] = {1, 2, 3, 4, 42, 6, 46, 8, 9, 12};
Color_t color[10] = {kRed + 2, kBlue + 1, 14, kGreen + 3, kMagenta + 3, kGreen + 1, kYellow + 2, 46, kBlue - 9, kViolet + 8};
double a1 = -4.125, a2 = 1.51, a3 = -3.29, teta = 1.41, n_knee = 205.5, chi2_NDF;
int bin_cent[11] = {400, 191, 134, 93, 63, 41, 25, 14, 7, 3, 2};
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

void Start(const char *fileadres , const char *current_mult , const char *outadres , int minNch , bool efficiencyFit , const char *fileadres2 , const char *current_mult2 )
{
    
    int Nch0 = minNch;
    TFile *file = new TFile(fileadres);
    TH1D *Gev = (TH1D *)file->Get(current_mult);
    Gev->Scale(1 / Gev->Integral(1, Gev->GetNbinsX(), "width"));
    bin_cent[0] = 1. * Gev->FindLastBinAbove();

    if (efficiencyFit == true)
    {
        int mediumNn = 0.5 * Gev->FindLastBinAbove();
        TFile *fileGl = new TFile(fileadres2);
        TH1D *Ideal = (TH1D *)fileGl->Get(current_mult2);
        Ideal->Scale(1 / Ideal->Integral(1, Ideal->GetNbinsX(), "width"));
        float Integr = Ideal->Integral(mediumNn, Ideal->GetNbinsX(), "width");
        int EfmediumNn = mediumNn * Scale(Ideal, Gev);
        Gev->Scale(Integr / Gev->Integral(EfmediumNn, Gev->GetNbinsX(), "width"));cout<<"start3"<<endl;
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
    TF1 *fc22 = new TF1("fit_func2", ftPn, 0.1, bin_cent[0], 5);
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
    chi2_NDF = fc2->GetChisquare() / fc2->GetNDF();

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
    for (int i = 0; i < 10; i++)
    {
        while (integr < cen[i])
        {
            integr = integ(n0, bin_cent[0]) / norm;
            n0 = n0 - 1;
        }
        bin_cent[i + 1] = n0 + 1;
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

    TF1 *fPb = new TF1("fPb", Pb, 0, bmax, 2);
    TH1F *Himpact[10];
    double Bmean[10], Bsigma[10], centbin[10], centEr[10];
    for (int cent = 0; cent < 10; cent++)
    {
        fPb->SetParameters(bin_cent[cent + 1], bin_cent[cent]);
        Himpact[cent] = new TH1F(Form("%sCENT%d_%d", "ImpactParametDist_", cent * 10, (cent + 1) * 10), "ImpactParametDist", 200, 0, bmax);
        Himpact[cent]->FillRandom("fPb", 200000);
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
    GrFit->SetMarkerColorAlpha(kRed + 1, 1);
    GrFit->SetLineColorAlpha(kRed + 1, 1);
    GrFit->SetLineWidth(2.);
    GrFit->SetName("Fit_B_Mean");
    GrFit->SetTitle("");
    GrFit->GetXaxis()->SetTitle("Centrality, %");
    GrFit->GetYaxis()->SetTitle("<b>, fm");
    GrFit->SetTitle("B vs Centraliry");
    GrFit->SetName("Fit_B_Mean");
    TCanvas *c = new TCanvas("canvas fit b", "B vs Centraliry", 650, 500);
    GrFit->Draw();
    GrFit->Write();
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

void PlotMeanbRatio(const char *fileadres = "/home/dim/FIT/data/UrQMD/7.7Gev/refMult_UrQMD_7.7gev_500k.root", const char *current_B_vs_mult = "hBvsRefMult")
{
    Rebin(1);
    double Bmean[10], Bsigma[10], cent[10], centEr[10];
    for (int i = 0; i < 10; i++)
    {
        Bmean[i] = bmean(bin_cent[i + 1], bin_cent[i]);
        Bsigma[i] = bsigma(bin_cent[i + 1], bin_cent[i]);
        cent[i] = i * 10 + 5;
        centEr[i] = 0;
    }
    TGraphErrors *GrFit = new TGraphErrors(10, cent, Bmean, centEr, Bsigma);
    GrFit->SetMarkerStyle(20);
    GrFit->SetMarkerSize(1.2);
    GrFit->SetMarkerColorAlpha(kRed + 1, 1);
    GrFit->SetLineColorAlpha(kRed + 1, 1);
    GrFit->SetLineWidth(2.);
    GrFit->SetName("Fit_B_Mean");

    cout << endl
         << "Centrality classes in impact parameter" << endl;
    cout << "Centr. class, bmin-bmax , bmean" << endl;
    for (int i = 1; i < 10; i++)
    {
        cout << Form("%d%s%d", (10 - i) * 10, "% - ", (11 - i) * 10) << "%, " << GrFit->Eval((10 - i) * 10) << " - " << GrFit->Eval((11 - i) * 10) << " fm, " << GrFit->Eval((11 - i) * 10 - 5) << " fm" << endl;
    }
    cout << Form("%d%s%d", 0, "% - ", 10) << "%, " << 0. << " - " << GrFit->Eval(10) << " fm, " << GrFit->Eval(5) << " fm" << endl;

    TFile *file_B_vs_mult = new TFile(fileadres);
    TH2F *MeanB;
    MeanB = (TH2F *)file_B_vs_mult->Get(current_B_vs_mult);
    TH1D *Hist[10];
    double y[10], ey[10];
    for (int i = 0; i < 10; i++)
    {
        Hist[i] = MeanB->ProjectionY(Form("%d", i), bin_cent[i + 1] + 1, bin_cent[i] + 1);
        y[i] = Hist[i]->GetMean();
        ey[i] = Hist[i]->GetStdDev();
    }
    TGraphErrors *grData;
    grData = new TGraphErrors(10, cent, y, centEr, ey);
    grData->GetYaxis()->SetRangeUser(-1, bmax);
    grData->GetXaxis()->SetLimits(0, 100);
    grData->GetXaxis()->SetNdivisions(505);
    grData->SetMarkerStyle(25);
    grData->SetMarkerSize(1.2);
    grData->SetMarkerColorAlpha(1, 1);
    grData->SetLineColorAlpha(1, 1);
    grData->GetYaxis()->SetTitle("<b>, fm");
    grData->SetTitle("");
    grData->SetName("model");

    TCanvas *RatCan = new TCanvas("RatCan", "RatCan", 600, 650);
    TPad *pad12 = new TPad("pad12", "pad12", 0, 0.3, 1, 1.0);
    pad12->SetBottomMargin(0.0);
    pad12->Draw();
    pad12->cd();
    grData->Draw("AP");
    GrFit->Draw("SAME P");
    TLegend *legdif = new TLegend(0.18, 0.7, 0.33, 0.85);
    legdif->SetBorderSize(0);
    legdif->AddEntry(grData, "Data", "p");
    legdif->AddEntry(GrFit, "Fit", "p");
    legdif->Draw();

    RatCan->cd(); // Go back to the main canvas before defining pad2
    TPad *pad22 = new TPad("pad22", "pad22", 0, 0, 1, 0.3);
    pad22->SetBottomMargin(0.3);
    pad22->SetTopMargin(0.0);
    pad22->Draw();
    pad22->cd();

    TGraphErrors *RatioMeanb = RatioGr(GrFit, grData, 0, 0.85, 100, 1.15);
    RatioMeanb->SetMarkerStyle(20);
    RatioMeanb->SetMarkerSize(1.2);
    RatioMeanb->SetMarkerColorAlpha(kRed + 1, 1);
    RatioMeanb->SetLineColorAlpha(kRed + 1, 1);
    RatioMeanb->SetTitle("");
    RatioMeanb->GetXaxis()->SetTitle("Centrality, %");
    RatioMeanb->GetXaxis()->SetLabelSize(0.08);
    RatioMeanb->GetXaxis()->SetTitleSize(0.09);
    RatioMeanb->GetYaxis()->SetTitle("Fit/Data");
    RatioMeanb->GetYaxis()->SetLabelSize(0.08);
    RatioMeanb->GetYaxis()->SetTitleSize(0.09);
    RatioMeanb->GetYaxis()->CenterTitle(true);
    RatioMeanb->GetYaxis()->SetNdivisions(505);
    RatioMeanb->GetXaxis()->SetNdivisions(505);
    RatioMeanb->GetYaxis()->SetTitleOffset(0.4);
    RatioMeanb->GetXaxis()->SetTickLength(0.08);
    RatioMeanb->SetLineWidth(2.);
    RatioMeanb->Draw("AP");
    TLine line;
    line.DrawLine(0, 1, 100, 1);
    line.SetLineStyle(2);
    line.DrawLine(0, 1.05, 100, 1.05);
    line.DrawLine(0, 0.95, 100, 0.95);
    d_outfile->cd();
    GrFit->Write();
    grData->Write();
    RatCan->Write();
    d_outfile->Close();
}

void Plot(int M, const char *current_fileadres, const char *current_B_vs_mult)
{
    TFile *file = new TFile(current_fileadres);
    TH2F *B_vs_mult;
    B_vs_mult = (TH2F *)file->Get(current_B_vs_mult);
    int NbinY = B_vs_mult->GetNbinsY();
    const int number = NbinY / M, numberY = bin_cent[0] / (1.1 * M);
    TH1D *Hist[M], *HistImpact[M], *histY;
    histY = B_vs_mult->ProjectionY("y", 10, 11);
    float width = histY->GetBinWidth(1), bi, bii;
    float meanN[M], sigmaN[M], ratio[M], b[M], cb[M], fn[M], sign[M], ration[M];
    float binNch[M], binEr[M], meanb[M], meanbFit[M];
    int binN0 = 1, binNn = binN0 + numberY;
    cout << "Nchmax=" << bin_cent[0] << endl;
    for (int i = 0; i < M; i++)
    {
        Hist[i] = B_vs_mult->ProjectionX(Form("%d_hist_", i), 1 + i * number, (i + 1) * number);
        meanN[i] = Hist[i]->GetMean();
        sigmaN[i] = Hist[i]->GetStdDev();
        ratio[i] = sigmaN[i] * sigmaN[i] / (meanN[i]);
        bi = (i)*number * width;

        bii = bi + number * width;
        b[i] = 0.5 * (bi + bii);
        cb[i] = (pi / sigma) * (bi * bi + bi * bii + bii * bii) / 3;

        fn[i] = n_knee * exp(a1 * cb[i] + a2 * pow(cb[i], 2) + a3 * pow(cb[i], 3));
        sign[i] = sqrt(fn[i] * teta);
        ration[i] = (sign[i] * sign[i]) / (fn[i]);

        HistImpact[i] = B_vs_mult->ProjectionY(Form("%d_hist_mult", i), 1 + binN0, 1 + binNn);
        meanb[i] = HistImpact[i]->GetMean();
        binNch[i] = 0.5 * (binN0 + binNn);
        meanbFit[i] = bmean(binNn, binN0);
        //cout<<"Nn="<<binNn<<endl;
        binEr[i] = 0;
        binNn += numberY;
        binN0 += numberY;
    }

    TCanvas *cdif = new TCanvas("cdif", "def Flow All", 900, 700);
    cdif->Divide(2, 2);
    cdif->cd(1);
    TGraphErrors *GrNchData = new TGraphErrors(M, b, meanN, binEr, binEr);
    GrNchData->SetTitle("");
    GrNchData->GetYaxis()->SetTitle("<N_{ch}>");
    GrNchData->GetXaxis()->SetTitle("b");
    GrNchData->GetXaxis()->SetLimits(0., bmax);
    GrNchData->SetMarkerStyle(25);
    GrNchData->Draw();
    TGraphErrors *GrNchFit = new TGraphErrors(M, b, fn, binEr, binEr);
    GrNchFit->SetMarkerStyle(20);
    GrNchFit->Draw("SAME P");
    TLegend *legdif2 = new TLegend(0.55, 0.65, 0.85, 0.85);
    legdif2->SetBorderSize(0);
    legdif2->AddEntry(GrNchData, "Data", "pe");
    legdif2->AddEntry(GrNchFit, "Fit", "pe");
    legdif2->Draw();
    cdif->cd(3);
    TGraphErrors *RatioNch = RatioGr(GrNchFit, GrNchData, 0, 0.75, bmax, 1.25);
    RatioNch->SetTitle("");
    RatioNch->GetYaxis()->SetTitle("Fit/Data");
    RatioNch->GetXaxis()->SetTitle("b");
    RatioNch->SetMarkerStyle(20);
    RatioNch->Draw("AP");
    TLine line;
    line.SetLineWidth(2);
    line.SetLineStyle(2);
    line.DrawLine(0, 1, bmax, 1);

    cdif->cd(2);
    TGraphErrors *GrNch_vs_b = new TGraphErrors(M, binNch, meanb, binEr, binEr);
    GrNch_vs_b->SetName("Nch_vs_b_data");
    GrNch_vs_b->GetYaxis()->SetTitle("<b>");
    GrNch_vs_b->GetXaxis()->SetTitle("N_{ch}");
    GrNch_vs_b->GetXaxis()->SetLimits(0., 1.1 * binNch[M - 1]);
    GrNch_vs_b->SetMarkerStyle(25);
    GrNch_vs_b->Draw();
    TGraphErrors *GrNch_vs_bFit = new TGraphErrors(M, binNch, meanbFit, binEr, binEr);
    GrNch_vs_bFit->SetMarkerStyle(20);
    GrNch_vs_bFit->SetName("Nch_vs_b_Fit");
    GrNch_vs_bFit->Draw("SAME P");
    TLegend *legdif4 = new TLegend(0.55, 0.65, 0.85, 0.85);
    legdif4->SetBorderSize(0);
    legdif4->AddEntry(GrNchData, "Data", "pe");
    legdif4->AddEntry(GrNch_vs_bFit, "Fit", "pe");
    legdif4->Draw();
    cdif->cd(4);
    TGraphErrors *RatioNch_vs_b = RatioGr(GrNch_vs_bFit, GrNch_vs_b, 0, 0.75, binNch[M - 1], 1.25);
    RatioNch_vs_b->SetTitle("");
    RatioNch_vs_b->GetYaxis()->SetTitle("Fit/Data");
    RatioNch_vs_b->GetXaxis()->SetTitle("N_{ch}");
    RatioNch_vs_b->SetMarkerStyle(20);
    RatioNch_vs_b->Draw("AP");
    line.DrawLine(0, 1, 1.1 * binNch[M - 1], 1);

    d_outfile->cd();
    GrNch_vs_b->Write();
    GrNch_vs_bFit->Write();
    d_outfile->Close();
    /*
    cdif->cd(3);
    TGraphErrors *GrSigmaData = new TGraphErrors(M, b, sigmaN,binEr,binEr);
    GrSigmaData->SetTitle("");
    GrSigmaData->GetYaxis()->SetTitle("#sigma");
    GrSigmaData->GetXaxis()->SetTitle("b");
    GrSigmaData->GetXaxis()->SetLimits(0., bmax);
    GrSigmaData->Draw();
    TGraph *GrSigmaFit = new TGraphErrors(M, b, sign,binEr,binEr);
    GrSigmaFit->SetMarkerStyle(20);
    GrSigmaFit->SetMarkerColor(2);
    GrSigmaFit->Draw("SAME P");
    TLegend *legdif1 = new TLegend(0.55, 0.65, 0.85, 0.85);
    legdif1->SetBorderSize(0);
    legdif1->AddEntry(GrSigmaData, "Data", "pe");
    legdif1->AddEntry(GrSigmaFit, "Fit", "pe");
    legdif1->Draw();

    cdif->cd(4);
    TGraphErrors *GrRatioData = new TGraphErrors(M, b, ratio,binEr,binEr);
    GrRatioData->SetTitle("");
    GrRatioData->GetYaxis()->SetTitle("#sigma^{2}/<N_{ch}>");
    GrRatioData->GetXaxis()->SetTitle("b");
    GrRatioData->GetXaxis()->SetLimits(0., bmax);
    GrRatioData->GetYaxis()->SetLimits(0., ratio[M - 1]);
    GrRatioData->Draw();
    TGraphErrors *GrRatioFit = new TGraphErrors(M, b, ration,binEr,binEr);
    GrRatioFit->SetMarkerStyle(20);
    GrRatioFit->SetMarkerColor(2);
    GrRatioFit->Draw("SAME P");
    TLegend *legdif3 = new TLegend(0.25, 0.65, 0.55, 0.85);
    legdif3->SetBorderSize(0);
    legdif3->AddEntry(GrRatioData, "Data", "pe");
    legdif3->AddEntry(GrRatioFit, "Fit", "pe");
    legdif3->Draw();*/
}

/*void GammaFit(const char *fileadres = "/home/dim/FIT/data/UrQMD/11.5Gev/refMult_UrQMD_11.5gev_500k.root", const char *current_mult = "hRefMultSTAR", const char *outadres = "/home/dim/FIT/FITout/urqmd_7_fitGamma.root", int minNch = 15)
{
    Start(fileadres, current_mult, outadres, minNch, false, fileadres, current_mult);
    //PlotMeanb();
    PlotMeanbRatio(fileadres, "hBvsRefMult");
}*/

void GammaFit(const char *fileadres = "/home/dim/FIT/data/UrQMD/11.5Gev/refMult_UrQMD_11.5gev_500k.root", const char *current_mult = "hRefMultSTAR", const char *outadres = "/home/dim/FIT/FITout/urqmd_7_fitGamma.root", int minNch = 20,bool efficiencyFit = false, const char *fileadres2 = "/home/dim/FIT/data/UrQMD/7.7Gev/refMult_UrQMD_7.7gev_500k.root", const char *current_mult2 = "hRefMultSTAR")
{
    Start(fileadres,current_mult,outadres,minNch,efficiencyFit,fileadres2,current_mult2);
    PlotMeanb();
}
