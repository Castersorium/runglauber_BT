void plot_dNdy_AB() {

    // ===== Open file =====
    TFile *f = TFile::Open("./AA_rapidityloss.root");
    if(!f || f->IsZombie()) {
        printf("Error: cannot open file AA_rapidityloss.root\n");
        return;
    }

    // ===== Get histograms =====
    TH1D *hA = (TH1D*) f->Get("h_dNdyA");
    TH1D *hB = (TH1D*) f->Get("h_dNdyB");
    if(!hA || !hB) {
        printf("Error: cannot find histograms h_dNdyA or h_dNdyB\n");
        return;
    }

    // ===== Create sum =====
    TH1D *hSum = (TH1D*) hA->Clone("hSum");
    hSum->Add(hB);

    // ===== Canvas =====
    TCanvas *c = new TCanvas("c", "dN/dy comparison", 900, 700);
    gStyle->SetOptStat(0);

    Int_t ci = TColor::GetFreeColorIndex();
    auto colorRed = new TColor(ci, 0.894, 0.145, 0.212);
    auto colorBlue = new TColor(ci+1, 0.341, 0.565, 0.988);    


    // ===== Style Settings =====
    hA->SetLineColor(kRed);
    hA->SetFillColor(ci);
    hA->SetLineWidth(2);

    hB->SetLineColor(kBlue);
    hB->SetFillColor(ci+1);
    hB->SetLineWidth(2);



    hSum->SetLineColor(kBlack);
    hSum->SetLineWidth(3);
    hSum->SetLineStyle(2); // dashed for clarity

    // ===== Set labels (optional but recommended) =====
    hA->GetXaxis()->SetTitle("y");
    hA->GetYaxis()->SetTitle("dN/dy");

    // ===== Draw =====
    hSum->Draw("HIST SAME");
    hA->Draw("HIST SAME");
    hB->Draw("HIST SAME");

    // ===== Legend =====
    auto leg = new TLegend(0.65,0.70,0.88,0.88);
    leg->AddEntry(hA,"A side","l");
    leg->AddEntry(hB,"B side","l");
    leg->AddEntry(hSum,"A + B","l");
    leg->Draw();

    // ===== Save =====
    c->SaveAs("dNdy_A_B_and_sum.pdf");
    c->SaveAs("dNdy_A_B_and_sum.root");

    printf("Saved: dNdy_A_B_and_sum.[pdf/root]\n");
}
