void compare_NcollA_threeModels()
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    // --- 文件名 ---
    const char* fnames[3] = {
        "AuAu_62p4_rapidityloss_05pCen_a3.root",
        "Au2rwAu2rw_62p4_rapidityloss_05pCen_a3.root",
        "Au197pnHFB14Au197pnHFB14_62p4_rapidityloss_05pCen_a3.root"
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
    c->Divide(3,1);

    // 颜色与样式
    Color_t colA  = kBlack;
    Color_t colP  = kRed+1;
    Color_t colN  = kBlue+1;

    for (int i = 0; i < 3; ++i) {
        c->cd(i+1);

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
        hA->GetXaxis()->SetRangeUser(0,10);
        hA->GetYaxis()->SetRangeUser(0,0.65);

        hA ->Draw("HIST");
        hAP->Draw("HIST same");
        hAN->Draw("HIST same");

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
    c->SaveAs("compare_NcollA_threeModels.root");
}
