#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>

#include <iostream>

/*------------------------------------------
  主函数
------------------------------------------*/
void draw_compare_BT() {

    gStyle->SetOptStat(0);

    /* ---------- 打开 NA49 实验数据 ---------- */
    TFile* fExp = TFile::Open(
        "./experimentData/NA49_2011_PbPb_17p3GeV_netP_dNdy.root",
        "READ"
    );
    if (!fExp || fExp->IsZombie()) {
        std::cerr << "Cannot open NA49 experiment file!" << std::endl;
        return;
    }

    /* ---------- centrality tags ---------- */
    const int Ncent = 5;
    const char* centTag[Ncent] = {"C0", "C1", "C2", "C3", "C4"};
    //const double scale[Ncent] =   {  0.5,    0.5,    0.5,   0.5,    0.5};// Half Isospin
    //const double scale[Ncent] = { 0.43, 0.4387,  0.454, 0.461,  0.596};// alpha = 4.0
    //const double scale[Ncent] = {0.510,  0.516, 0.5336, 0.536, 0.6989};// alpha = 5.0
    const char* centRange[Ncent] = {
        "0-5%", "5-12.5%", "12.5-23.5%", "23.5-33.5%", "33.5-43.5%"
    };

    /* ---------- alpha values ---------- */
    const int Nalpha = 3;
    const double alphaVal[Nalpha] = {3.0, 4.0, 5.0};
    const int alphaColor[Nalpha] = {kBlue+1, kGreen+2, kRed+1};

    /* ---------- Canvas ---------- */
    TCanvas* c = new TCanvas("c", "NA49 vs Model", 3000, 1000);
    c->Divide(5, 1, 0.001, 0.001);

    /* ---------- Loop over centrality ---------- */
    for (int ic = 0; ic < Ncent; ++ic) {

        c->cd(ic + 1);
        gPad->SetTicks();
        gPad->SetLeftMargin(0.25);
        gPad->SetRightMargin(0.005);
        gPad->SetBottomMargin(0.15);

        if (ic != 0) gPad->SetLeftMargin(0.001);

        /* ---------- 空 frame ---------- */
        TH1F* frame = gPad->DrawFrame(
            -4.0, 0.0,
             4.0, 0.2   // 你可以之后自己再调
        );
        frame->SetTitle(Form("Centrality %s", centRange[ic]));
        frame->GetXaxis()->SetTitle("y");
        frame->GetXaxis()->CenterTitle();
        frame->GetYaxis()->SetTitle("(1/<N_{w}>) dN/dy");

        frame->GetXaxis()->SetTitleSize(0.07);
        frame->GetYaxis()->SetTitleSize(0.07);
        frame->GetXaxis()->SetLabelSize(0.06);
        frame->GetYaxis()->SetLabelSize(0.06);

        /* ---------- 画 NA49 ---------- */
        TString grName = Form(
            "NA49_2011_17p3GeV_netP_dNdy_overNw_%s",
            centTag[ic]
        );

        auto* grExp =
            (TGraphErrors*) fExp->Get(grName);

        if (!grExp) {
            std::cerr << "Missing graph: " << grName << std::endl;
        } else {
            grExp->SetMarkerStyle(20);
            grExp->SetMarkerSize(1.4);
            grExp->SetMarkerColor(kRed);
            grExp->SetLineColor(kRed);
            grExp->Draw("P SAME");
        }

        /* ---------- 画模型 ---------- */
        for (int ia = 0; ia < Nalpha; ++ia) {

            TString fname = Form(
                "./rapidity_dNdy/"
                "dNdy_PbpnrwPbpnrw_17p3_2011_alpha%.2f_cent%s.root",
                alphaVal[ia],
                (ic==0) ? "0.0_5.0" :
                (ic==1) ? "5.0_12.5" :
                (ic==2) ? "12.5_23.5" :
                (ic==3) ? "23.5_33.5" :
                          "33.5_43.5"
            );

            TFile* fModel = TFile::Open(fname, "READ");
            if (!fModel || fModel->IsZombie()) {
                std::cerr << "Cannot open " << fname << std::endl;
                continue;
            }

            TH1F* h =
                (TH1F*) fModel->Get("h_dNdy_over_Nw");

            if (!h) {
                std::cerr << "Missing histogram in " << fname << std::endl;
                fModel->Close();
                continue;
            }
            h->SetDirectory(0);
            h->Scale(scale[ic]);
            h->SetLineColor(alphaColor[ia]);
            h->SetLineWidth(2);
            h->SetFillStyle(0);
            h->Draw("HIST C SAME");

            fModel->Close();
        }

        /* ---------- Legend（只在第一个 pad 画） ---------- */
        if (ic == 0) {
            TLegend* leg = new TLegend(0.30, 0.65, 0.90, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.055);

            leg->AddEntry(grExp, "NA49 (2011)", "p");
            TLegendEntry* e =
            leg->AddEntry((TObject*)0, "Scale: #alpha = 5 at y #approx 0 ", "");
            e->SetTextColor(kRed+1);
            leg->Draw();
        }
    }

    c->SaveAs("compare_NA49_vs_model_BT_ScaleA5.root");
}
