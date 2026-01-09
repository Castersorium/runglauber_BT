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
    //h->SetFillColor(color);
    h->SetLineColor(color);
    h->SetLineWidth(3);
}

/*------------------------------------------
  主程序
------------------------------------------*/
void comparePbPb() {

    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    /* ---------- 文件与参数 ---------- */
    const int N = 5;
    TString files[N] = {
        "./rapidity_dNdy/dNdy_PbpnrwPbpnrw_17p3_alpha2.00_cent0_5.root",
        "./rapidity_dNdy/dNdy_PbpnrwPbpnrw_17p3_alpha3.00_cent0_5.root",
        "./rapidity_dNdy/dNdy_PbpnrwPbpnrw_17p3_alpha4.00_cent0_5.root",
        "./rapidity_dNdy/dNdy_PbpnrwPbpnrw_17p3_alpha5.00_cent0_5.root",
        "./rapidity_dNdy/dNdy_PbpnrwPbpnrw_17p3_alpha6.00_cent0_5.root"
    };

    //PbPb 17.3 GeV 0~5% ,

    TString labels[N] = {
        "Scaled #alpha = 2.00",
        "Scaled #alpha = 3.00",
        "Scaled #alpha = 4.00",
        "Scaled #alpha = 5.00",
        "Scaled #alpha = 6.00"
    };

    int colors[N]  = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+2};
    int markers[N] = {20, 21, 22, 33, 29};
    double scale[N] = {0.3668, 0.3679, 0.4284, 0.5073, 0.5952};

    //double scale[N] = {0.3668, 0.3679, 0.3906, 0.4284, 0.4640}; // 2 3 3.5 4 4.5

    /* ---------- Canvas ---------- */
    auto* c1 = new TCanvas("c1", "Rapidity loss comparison", 1200, 800);
    c1->SetTicks();
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.12);

    /* ---------- Frame ---------- */
    TH1F* frame = c1->DrawFrame(
        -1.1, 15,  //xmin ymin
         3.0, 80   //xmax ymax
    );

    frame->SetTitle("PbPb 17.3 GeV 0~5%");
    frame->GetXaxis()->SetTitle("y");
    frame->GetYaxis()->SetTitle("dN/dy");

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);

    /* ---------- Legend ---------- */
    auto* leg = new TLegend(0.15, 0.58, 0.48, 0.88);
    //leg->AddEntry((TObject*)0, "Nucleon Scaled #times 0.48","");  
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

        TH1F* h = (TH1F*)f->Get("h_dNdy");
        if (!h) {
            std::cerr << "Histogram h_dNdy not found in "
                      << files[i] << std::endl;
            f->Close();
            continue;
        }

        // 防止文件关闭后 histogram 消失
        h->SetDirectory(0);
        //h->Scale(0.5);
        //h->Scale(scale[i]);

        styleHist(h, colors[i], markers[i]);
        h->Draw("HIST C SAME");
 
        leg->AddEntry(h, labels[i], "l p");
        

        f->Close();
    }
      
    leg->Draw();

    // Add Exp Data comparison
    TFile* fexp = TFile::Open("./experimentData/ExpData_dNdy.root", "READ");
    if (!fexp || fexp->IsZombie()) {
        std::cerr << "Cannot open expData.root" << std::endl;
    } else {
    
        const int Nexp = 5;
        TString expNames[Nexp] = {
            "NA49_17p3GeV_1999",
            "NA49_17p3GeV_2011",
            "NA49_17p3GeV_1999_mirror",
            "NA49_17p3GeV_2011_mirror"
        };
    
        TString expLabels[Nexp] = {
            "NA49(1999) Pb+Pb 17.3 GeV",
            "NA49(2021) Pb+Pb 17.3 GeV Mirror",
            "NA49(1999) Pb+Pb 17.3 GeV",
            "NA49(2021) Pb+Pb 17.3 GeV Mirror"
        };
    
        for (int i = 0; i < Nexp; ++i) {
        
            TGraphErrors* gr =
                (TGraphErrors*)fexp->Get(expNames[i]);
        
            if (!gr) {
                std::cerr << "Cannot find " << expNames[i]
                          << " in expData.root" << std::endl;
                continue;
            }
        
        
            // 直接画，不改任何风格
            gr->Draw("P SAME");
        
            //leg->AddEntry(gr, expLabels[i], "p");
        }
    
        fexp->Close();
    }

    


    //c1->SaveAs("PbPb_17p3_Valpha_005Cen.pdf");
    c1->SaveAs("PbPb_17p3_Valpha_005Cen.root");
}
