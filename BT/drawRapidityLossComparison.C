#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TGraphErrors.h>

#include <iostream>
#include <vector>

/*------------------------------------------
  统一设置 histogram 风格
------------------------------------------*/
void styleHist(TH1* h, int color, int marker) {
    //h->SetMarkerStyle(marker);
    h->SetMarkerSize(1.2);
    h->SetMarkerColor(color);
    h->SetFillColor(color);
    h->SetLineColor(color);
    h->SetLineWidth(2);
}

/*------------------------------------------
  主程序
------------------------------------------*/
void drawRapidityLossComparison() {

    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    /* ---------- 文件与参数 ---------- */
    const int N = 5;
    TString files[N] = {
        "AA62p4_rapidityloss_010pCen_a1.root",
        "AA62p4_rapidityloss_010pCen_a2.root",
        "AA62p4_rapidityloss_010pCen_a3.root",
        "AA62p4_rapidityloss_010pCen_a4.root",
        "AA62p4_rapidityloss_010pCen_a5.root"
    };

    TString labels[N] = {
        "MCM 62.4 GeV 0~10% , #alpha = 1.0",
        "MCM 62.4 GeV 0~10% , #alpha = 2.0",
        "MCM 62.4 GeV 0~10% , #alpha = 3.0",
        "MCM 62.4 GeV 0~10% , #alpha = 4.0",
        "MCM 62.4 GeV 0~10% , #alpha = 5.0"
    };

    int colors[N]  = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+2};
    int markers[N] = {20, 21, 22, 33, 29};

    /* ---------- Canvas ---------- */
    auto* c1 = new TCanvas("c1", "Rapidity loss comparison", 1200, 800);
    c1->SetTicks();
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.12);

    /* ---------- Frame ---------- */
    TH1F* frame = c1->DrawFrame(
        -7.5, 0.0,
         0.5, 0.7   // 根据你 h_dNdDyA 实际范围可再调
    );

    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("#Deltay=y-y_{beam}");
    frame->GetYaxis()->SetTitle("(N_{part}/2)dN^{#it{B-#bar{B}}}/d(y-y_{beam})");

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);

    /* ---------- Legend ---------- */
    auto* leg = new TLegend(0.15, 0.35, 0.48, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    /* ---------- 读文件并画图 ---------- */
    for (int i = 0; i < N; ++i) {

        TFile* f = TFile::Open(files[i], "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Cannot open file: " << files[i] << std::endl;
            continue;
        }

        TH1F* h = (TH1F*)f->Get("h_dNdDyA");
        if (!h) {
            std::cerr << "Histogram dNdDyA not found in "
                      << files[i] << std::endl;
            f->Close();
            continue;
        }

        // 防止文件关闭后 histogram 消失
        h->SetDirectory(0);

        styleHist(h, colors[i], markers[i]);
        h->Draw("E1 SAME");
 
        leg->AddEntry(h, labels[i], "l p");

        f->Close();
    }

    leg->Draw();

// TFile* fexp = TFile::Open("expData.root", "READ");
// if (!fexp || fexp->IsZombie()) {
//     std::cerr << "Cannot open expData.root" << std::endl;
// } else {

//     const int Nexp = 5;
//     TString expNames[Nexp] = {
//         "NA49_17p3GeV",
//         "BRAHMS_62p4GeV",
//         "BRAHMS_200GeV",
//         "STAR_62p4GeV",
//         "STAR_200GeV"
//     };

//     TString expLabels[Nexp] = {
//         "NA49 Pb+Pb 17.3 GeV",
//         "BRAHMS Au+Au 62.4 GeV",
//         "BRAHMS Au+Au 200 GeV",
//         "STAR Au+Au 62.4 GeV",
//         "STAR Au+Au 200 GeV"
//     };

//     for (int i = 0; i < Nexp; ++i) {

//         TGraphErrors* gr =
//             (TGraphErrors*)fexp->Get(expNames[i]);

//         if (!gr) {
//             std::cerr << "Cannot find " << expNames[i]
//                       << " in expData.root" << std::endl;
//             continue;
//         }


//         // 直接画，不改任何风格
//         gr->Draw("P SAME");

//         leg->AddEntry(gr, expLabels[i], "p");
//     }

//     fexp->Close();
// }

    


    c1->SaveAs("AA62p4_rapidityloss_010pCen_compare.pdf");
    c1->SaveAs("AA62p4_rapidityloss_010pCen_compare.root");
}
