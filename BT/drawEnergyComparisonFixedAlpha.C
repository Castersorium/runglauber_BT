#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>

#include <iostream>

/*------------------------------------------
  统一设置 histogram 风格
------------------------------------------*/
void styleHist(TH1* h, int color, int marker) {
    //h->SetMarkerStyle(marker);
    h->SetMarkerSize(1.3);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetLineWidth(2);
}

/*------------------------------------------
  主程序
------------------------------------------*/
void drawEnergyComparisonFixedAlpha() {

    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    TFile* fout = new TFile(
        "dNdDyA_energyComparison_fixedAlpha.root",
        "RECREATE"
    );

    /* ---------- 能量信息 ---------- */
    const int Nenergy = 3;
    TString energyTag[Nenergy] = {
        "PbPb17p3",
        "AA62p4",
        "AA200"
    };

    TString energyLabel[Nenergy] = {
        "Pb+Pb 17.3 GeV",
        "Au+Au 62.4 GeV",
        "Au+Au 200 GeV"
    };

    int colors[Nenergy]  = {kBlack, kRed+1, kBlue+1};
    int markers[Nenergy] = {20, 21, 22};

    /* ---------- alpha loop ---------- */
    for (int a = 1; a <= 5; ++a) {

        TString cname;
        cname.Form("c_alpha%d", a);

        auto* c = new TCanvas(
            cname,
            Form("Fixed a = %d comparison", a),
            800, 600
        );
        c->SetTicks();
        c->SetLeftMargin(0.13);
        c->SetBottomMargin(0.12);

        /* ---------- Frame ---------- */
        TH1F* frame = c->DrawFrame(
            -7.5, 0.0,
             0.5, 0.7    // 按你模型最大值调整
        );

        frame->SetTitle("");
        frame->GetXaxis()->SetTitle("#Deltay=y-y_{beam}");
        frame->GetYaxis()->SetTitle("(N_{part}/2)dN^{#it{N-#bar{N}}}/d(y-y_{beam})");

        frame->GetXaxis()->SetTitleSize(0.045);
        frame->GetYaxis()->SetTitleSize(0.045);
        frame->GetXaxis()->SetLabelSize(0.040);
        frame->GetYaxis()->SetLabelSize(0.040);

        /* ---------- Legend ---------- */
        auto* leg = new TLegend(0.15, 0.65, 0.48, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);

        /* ---------- 能量 loop ---------- */
        for (int ie = 0; ie < Nenergy; ++ie) {

            TString fname;
            fname.Form(
                "%s_rapidityloss_05pCen_a%d.root",
                energyTag[ie].Data(),
                a
            );

            TFile* f = TFile::Open(fname, "READ");
            if (!f || f->IsZombie()) {
                std::cerr << "Cannot open " << fname << std::endl;
                continue;
            }

            TH1F* h = (TH1F*)f->Get("h_dNdDyA");
            if (!h) {
                std::cerr << "dNdDyA not found in "
                          << fname << std::endl;
                f->Close();
                continue;
            }

            // 防止 file 关闭后 histogram 消失
            h->SetDirectory(0);

            styleHist(h, colors[ie], markers[ie]);
            h->Draw("E1 SAME");

            leg->AddEntry(h, energyLabel[ie], "l p");


            


            f->Close();
        }

        leg->AddEntry((TObject*)0,
        Form("#alpha = %d.0", a),
        "");
        leg->Draw();
        fout->cd();
        c->Write();

        // /* ---------- 输出 ---------- */
        // c->SaveAs(Form("dNdDyA_fixedAlpha%d_energyCompare.pdf", a));
        // c->SaveAs(Form("dNdDyA_fixedAlpha%d_energyCompare.root", a));
    }
    fout->Write();
    fout->Close();
}
