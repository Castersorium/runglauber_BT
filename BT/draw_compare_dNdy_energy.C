#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>

void draw_compare_dNdy_energy() {

    gStyle->SetOptStat(0);

    /* ---------- Experimental data ---------- */
    TFile* fExp = TFile::Open(
       // "./experimentData/ExpData_netP_dNdy.root", "READ"
         "./experimentData/ExpData_netP_dNdy_overHalfNw.root", "READ"
    );
    if (!fExp || fExp->IsZombie()) {
        std::cerr << "Cannot open experiment data file!" << std::endl;
        return;
    }

    /* ---------- Energy setup ---------- */
    const int Nenergy = 3;
    const char* energyLabel[Nenergy] = {
        "17.3 GeV", "62.4 GeV", "200 GeV"
    };

    /* ---------- Alpha ---------- */
    const int Nalpha = 4;
    const double alphaVal[Nalpha] = {3.0, 3.5, 4.0, 5.0};
    const double scale[Nalpha] = { 0.510, 0.510,  0.510, 0.510};
    const int alphaColor[Nalpha] = {kBlue+1, kMagenta, kGreen+2, kRed+1};

    /* ---------- Canvas ---------- */
    TCanvas* c = new TCanvas("c", "net-proton dN/dy over Nw/2", 2400, 700);
    c->Divide(3,1,0.001,0.001);

    /* =========================================================
       Loop over energies
       ========================================================= */
    for (int ie = 0; ie < Nenergy; ++ie) {

        c->cd(ie+1);
        gPad->SetTicks();
        gPad->SetLeftMargin(0.22);
        gPad->SetBottomMargin(0.15);
        if (ie != 0) gPad->SetLeftMargin(0.05);

        TH1F* frame = gPad->DrawFrame(
            -6.0, 0.0,
             6.0, 0.3
        );
        frame->SetTitle(energyLabel[ie]);
        frame->GetXaxis()->SetTitle("y");
        frame->GetXaxis()->CenterTitle();
        //frame->GetYaxis()->SetTitle("dN(p-#bar{p})/dy");
        frame->GetYaxis()->SetTitle("(2/<N_{w}>)dN/dy");

        frame->GetXaxis()->SetTitleSize(0.07);
        frame->GetYaxis()->SetTitleSize(0.07);
        frame->GetXaxis()->SetLabelSize(0.06);
        frame->GetYaxis()->SetLabelSize(0.03);


        TLegend* leg = new TLegend(0.28,0.65,0.88,0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        //leg->SetTextSize(0.055);

        /* ---------- Experimental graphs ---------- */
        if (ie == 0) {
            // NA49 1999
            auto* g99 = (TGraphErrors*)
                fExp->Get("NA49_1999_17p3GeV_netP_dNdy_overNw");
            if (g99) {
                g99->SetMarkerStyle(20);
                g99->SetMarkerSize(1.3);
                g99->SetMarkerColor(kBlack);
                g99->SetLineColor(kBlack);
                g99->Draw("P SAME");
            }

            // NA49 2011
            auto* g11 = (TGraphErrors*)
                fExp->Get("NA49_2011_17p3GeV_netP_dNdy_overNw_C0");
            if (g11) {
                g11->SetMarkerStyle(20);
                g11->SetMarkerSize(1.3);
                g11->SetMarkerColor(kRed+1);
                g11->SetLineColor(kRed+1);
                g11->Draw("P SAME");
            }

            leg->AddEntry(g99,"NA49 1999 0~5%","p");
            leg->AddEntry(g11,"NA49 2011 0~5%","p");
        }

        if (ie == 1) {
            auto* g62B = (TGraphErrors*)
                fExp->Get("BRAHMS_AuAu_62p4GeV_010_netP_dNdy_overNw");
            if (g62B) {
                g62B->SetMarkerStyle(20);
                g62B->SetMarkerSize(1.8);
                g62B->SetMarkerColor(kBlack);
                g62B->Draw("P SAME");
            }
            auto* g62S = (TGraphErrors*)
            fExp->Get("STAR_AuAu_62p4GeV_005_netP_dNdy_overNw");
            if (g62S) {
            g62S->SetMarkerStyle(29);
            g62S->SetMarkerSize(2.0);
            g62S->SetMarkerColor(kRed);
            g62S->Draw("P SAME");
            }

            leg->AddEntry(g62B,"BRAHMS 0~10%","p");
            leg->AddEntry(g62S,"STAR 0~5%","p");
        }

        if (ie == 2) {
            auto* g200B = (TGraphErrors*)
                fExp->Get("BRAHMS_AuAu_200GeV_005_netP_dNdy_overNw");
            if (g200B) {
                g200B->SetMarkerStyle(20);
                g200B->SetMarkerSize(1.8);
                g200B->SetMarkerColor(kBlack);
                g200B->Draw("P SAME");
            }
            auto* g200S = (TGraphErrors*)
                fExp->Get("STAR_AuAu_200GeV_005_netP_dNdy_overNw");
            if (g200S) {
                g200S->SetMarkerStyle(29);
                g200S->SetMarkerSize(2.0);
                g200S->SetMarkerColor(kRed);
                g200S->Draw("P SAME");
        }

        leg->AddEntry(g200B,"BRAHMS 0~5%","p");
        leg->AddEntry(g200S,"STAR   0~5%","p");
        }

        /* ---------- Model ---------- */
        for (int ia = 0; ia < Nalpha; ++ia) {

            TString fname;
            if (ie == 0)
                fname = Form("./rapidity_dNdy/dNdy_PbpnrwPbpnrw_17p3_1999_alpha%.2f_cent0.0_5.0.root", alphaVal[ia]);
            if (ie == 1)
                fname = Form("./rapidity_dNdy/dNdy_Au197pnHFB14Au197pnHFB14_62p4_BRAH_alpha%.2f_cent0.0_10.0.root", alphaVal[ia]);
            if (ie == 2)
                fname = Form("./rapidity_dNdy/dNdy_Au197pnHFB14Au197pnHFB14_200_BRAH_alpha%.2f_cent0.0_5.0.root", alphaVal[ia]);

            TFile* fM = TFile::Open(fname,"READ");
            if (!fM || fM->IsZombie()) continue;

            TH1F* h = (TH1F*) fM->Get("h_dNdy_over_Nw");
            if (!h) { fM->Close(); continue; }

            h->SetDirectory(0);
            h->SetLineColor(alphaColor[ia]);
            // h->Scale(scale[ia]);
            //h->Scale(0.5);
            h->SetLineWidth(2);
            h->SetFillStyle(0);
            h->Draw("HIST C SAME");

            leg->AddEntry(h,Form("#alpha = %.2f",alphaVal[ia]),"l");


            fM->Close();
        }


            // leg->AddEntry((TObject*)0,"Experiment","p");


            // leg->AddEntry((TObject*)0,"Model","l");


            leg->Draw();

    }

    c->SaveAs("compare_netP_dNdy_overHalfNw_energy.pdf");
}
