TH2D *h2;
TH1D *projX;
TH1D *projY;
TPad *right_pad, *top_pad;

// =============================================================
void plot_DeltaY_Ncoll()
{
    gStyle->SetOptStat(0);

    // ======== Load file ========
    TFile *f = TFile::Open("../pA_rapidityloss.root");
    if(!f || f->IsZombie()) { printf("File error\n"); return; }
    h2 = (TH2D*) f->Get("h2_DeltaY_Ncoll");
    if(!h2) { printf("Histogram not found!\n"); return; }

    // ======== Make projections ========
    projX = h2->ProjectionX("projX");
    projY = h2->ProjectionY("projY");

    // ======== Canvas + Pads ========
    auto c1 = new TCanvas("c1", "DeltaY vs. Ncoll with projections", 1200, 1200);

    // center (2D)
    TPad *center_pad = new TPad("center_pad", "center_pad", 0.0, 0.0, 0.6, 0.6);
    center_pad->Draw();

    // right (Y-projection)
    right_pad = new TPad("right_pad", "right_pad", 0.55, 0.0, 1.0, 0.6);
    right_pad->Draw();

    // top (X-projection)
    top_pad = new TPad("top_pad", "top_pad", 0.0, 0.55, 0.6, 1.0);
    top_pad->Draw();

    // ======== Draw center pad ========
    center_pad->cd();
    gStyle->SetPalette(1);
    h2->Draw("COLZ");

    // ======== Draw X projection ========
    top_pad->cd();
    projX->SetFillColor(kRed+1);
    projX->SetTitle("");
    //projX->GetYaxis()->SetTitle("Counts");
    projX->GetXaxis()->SetTitle("N_{coll}");
    projX->Draw("BAR");

    // ======== Draw Y projection (horizontal bar) ========
    right_pad->cd();
    projY->SetFillColor(kBlue+1);
    projY->SetTitle("");
    //projY->GetXaxis()->SetTitle("Counts");
    //projY->GetYaxis()->SetTitle("#Delta y");
    projY->Draw("HBAR");

    // ======== Add linked zoom ========
    auto ex = new TExec("zoom", "ZoomExec()");
    h2->GetListOfFunctions()->Add(ex);

    // ======== Save canvas ========
    c1->SaveAs("DeltaY_Ncoll_projection.root");
    c1->SaveAs("DeltaY_Ncoll_projection.pdf");
    printf("Saved: DeltaY_Ncoll_projection.pdf\n");
}

// =============================================================
// When user zooms the 2D histogram, update projection ranges
void ZoomExec()
{
    // X range
    int xfirst = h2->GetXaxis()->GetFirst();
    int xlast  = h2->GetXaxis()->GetLast();
    double xmin = h2->GetXaxis()->GetBinLowEdge(xfirst);
    double xmax = h2->GetXaxis()->GetBinUpEdge(xlast);
    projX->GetXaxis()->SetRangeUser(xmin, xmax);
    top_pad->Modified();

    // Y range
    int yfirst = h2->GetYaxis()->GetFirst();
    int ylast  = h2->GetYaxis()->GetLast();
    double ymin = h2->GetYaxis()->GetBinLowEdge(yfirst);
    double ymax = h2->GetYaxis()->GetBinUpEdge(ylast);
    projY->GetXaxis()->SetRangeUser(ymin, ymax);
    right_pad->Modified();
}
