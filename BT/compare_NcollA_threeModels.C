void compare_NcollA_threeModels()
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    // --- 文件名 ---
    const char* fnames[3] = {
                            "./rapidityLossData/AuAu_62p4_rapidityloss_7080pCen_a3.root",
                      "./rapidityLossData/Au2rwAu2rw_62p4_rapidityloss_7080pCen_a3.root",
        "./rapidityLossData/Au197pnHFB14Au197pnHFB14_62p4_rapidityloss_7080pCen_a3.root"
    };

    const char* titles[3] = {
        "Au + Au",
        "Au2rw + Au2rw",
        "Au197pn(HFB14) + Au197pn(HFB14)"
    };

    // --- 画布 ---
    TCanvas* c = new TCanvas("c_NcollA_compare",
                             "NcollA comparison",
                             2400, 800);

    TPad* pad1 = new TPad("pad1","", 0.00, 0.00, 0.33, 1.00);
    TPad* pad2 = new TPad("pad2","", 0.33, 0.00, 0.66, 1.00);
    TPad* pad3 = new TPad("pad3","", 0.66, 0.00, 1.00, 1.00);

    pad1->SetLeftMargin(0.16);
    pad1->SetRightMargin(0.02);
    pad1->SetBottomMargin(0.14);

    pad2->SetLeftMargin(0.06);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.14);

    pad3->SetLeftMargin(0.06);
    pad3->SetRightMargin(0.10);
    pad3->SetBottomMargin(0.14);

    pad1->Draw();
    pad2->Draw();
    pad3->Draw();

    // 颜色与样式
    Color_t colA  = kBlack;
    Color_t colP  = kRed+1;
    Color_t colN  = kBlue+1;

    for (int i = 0; i < 3; ++i) {
        if (i == 0) pad1->cd();
        if (i == 1) pad2->cd();
        if (i == 2) pad3->cd();


        TFile* f = TFile::Open(fnames[i], "READ");
        if (!f || f->IsZombie()) {
            Error("compare_NcollA_threeModels",
                  "Cannot open file %s", fnames[i]);
            continue;
        }

        auto hA  = (TH1F*) f->Get("h_NcollA");
        auto hAP = (TH1F*) f->Get("h_NcollAP");
        auto hAN = (TH1F*) f->Get("h_NcollAN");

        if (!hA || !hAP || !hAN) {
            Error("compare_NcollA_threeModels",
                  "Histogram missing in %s", fnames[i]);
            f->Close();
            continue;
        }

        // --- 若尚未归一化，可取消注释 ---
        /*
        hA ->Scale(1.0 / hA ->Integral("width"));
        hAP->Scale(1.0 / hAP->Integral("width"));
        hAN->Scale(1.0 / hAN->Integral("width"));
        */

        // 样式
        hA ->SetLineColor(colA);
        hAP->SetLineColor(colP);
        hAN->SetLineColor(colN);

        hA ->SetLineWidth(2);
        hAP->SetLineWidth(2);
        hAN->SetLineWidth(2);

        hA->SetTitle(titles[i]);
        hA->GetXaxis()->SetTitle("N_{coll}");
        hA->GetYaxis()->SetTitle("P(N_{coll})");
        hA->GetYaxis()->SetTitleOffset(1.4);
        hA->GetYaxis()->SetLabelSize(0.045);
        hA->GetXaxis()->SetLabelSize(0.045);
        hA->GetXaxis()->SetTitleSize(0.05);
        hA->GetYaxis()->SetTitleSize(0.05);

        hA->GetXaxis()->SetRangeUser(0,20);
        hA->GetYaxis()->SetRangeUser(0,0.650);

        hA ->Draw("E1");
        hAP->Draw("E1 same");
        hAN->Draw("E1 same");

        // legend
        TLegend* leg = new TLegend(0.60, 0.70, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hA,  "All nucleons", "l");
        leg->AddEntry(hAP, "Protons",      "l");
        leg->AddEntry(hAN, "Neutrons",      "l");
        leg->Draw();

        //f->Close();
    }

    c->Update();
    c->SaveAs("compare_NcollA_threeModels_7080pCen_a3.pdf");
}
