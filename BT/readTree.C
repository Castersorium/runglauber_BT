#include <TString.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH3F.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <stdexcept>

#include <stdio.h>
#include <iostream>

#include "include/EXP_dNdy.h"

struct CentNpartEntry {
    const char* energy;
    double centMin;
    double centMax;
    double npartMin;
    double npartMax;
};

static const CentNpartEntry gCentNpartTable[] = {

    // ---------- Au+Au 200 GeV BRAHMS----------
    //{"200",   0,  5, 357-8, 357+8},

    // ---------- Au+Au 200 GeV STAR----------
    {"200",   0,  5, 350.0-2.4, 350.0+2.4},
    {"200",   5, 10, 298.6-4.1, 298.6+4.1},
    {"200",  10, 20, 234.3-4.6, 234.3+4.6},
    {"200",  70, 80,  15.7-2.6,  15.7+2.6},

    // ---------- Au+Au 62.4 GeV BRAHMS----------
    {"62p4",  0,  10, 314-8, 314+8},

    // ---------- Au+Au 62.4 GeV STAR----------
    {"62p4",  0,  5, 346.5-2.8, 346.5+2.8},
    {"62p4",  5, 10, 293.9-4.2, 293.9+4.2},
    {"62p4", 10, 20, 229.8-4.6, 229.8+4.6},
    {"62p4", 70, 80,  15.3-2.4,  15.3+2.4},



    // ---------- Pb+Pb 17.3 GeV ----------
    {"17p3",  0,  5, 352-12, 352+12}//Phys. Rev. Lett. 82, 2471–2475 (1999), arXiv:nucl-ex/9810014.
};

bool LookupNpartRange(
    const TString& energy,
    double centMin,
    double centMax,
    double& npartMin,
    double& npartMax
)
{
    for (const auto& e : gCentNpartTable) {
        if (energy == e.energy &&
            centMin == e.centMin &&
            centMax == e.centMax) {

            npartMin = e.npartMin;
            npartMax = e.npartMax;
            return true;
        }
    }

    std::cout << "Error: centrality bin not found in table"
              << std::endl;
    return false;
}

//root -l -b 'readTree.C("Pbpnrw","Pbpnrw","17p3",3.0,0,5)'
//root -l -b 'readTree.C("Au197pnHFB14","Au197pnHFB14","62p4",3.0,0,5)'
//root -l -b 'readTree.C("Au197pnHFB14","Au197pnHFB14","200",3.0,0,5)'

// for alpha in 1.0 2.0 3.0 4.0 5.0; do
//   root -l -b -q "saveTreeAA.C(\"Au197pnHFB14\",\"Au197pnHFB14\",\"62p4\",${alpha},0,5)"
// done

void readTree(
    const char* projectile,
    const char* target,
    const char* energy,
    double alpha,
    int centMin,
    int centMax
) {

    TString Projectile = projectile;//"Au197pnHFB14"; //Au Au2rw  Au197pnHFB14
    TString Target = target;
    TString Energy = energy;
    TString System = Projectile + Target+ "_" + Energy;

    double npartMin, npartMax;
    if (!LookupNpartRange(Energy, centMin, centMax,
                          npartMin, npartMax)) {
        return;
    }

    // ---------- beam rapidity ----------
    double y_beam = -1;
    if      (Energy == "200")  y_beam = 5.36;
    else if (Energy == "62p4") y_beam = 4.20;
    else if (Energy == "17p3") y_beam = 2.90;
    else {
        std::cout << "Error: Unknown energy" << std::endl;
        return;
    }

    //double alpha = 3.0; 

    // 1. 打开 TGlauber 输出文件
    TFile *fin = TFile::Open("./rapidityTree/tree_" + System + "_alpha4.00_cent0_5.root");
    if (!fin || fin->IsZombie()) { 
        std::cout << "Error: Cannot open file\n"; 
        return; 
    }

    // TString outName;
    // outName.Form(
    //     "tree_%s%s_%s_alpha%.2f_cent%d_%d.root",
    //     projectile,
    //     target,
    //     energy,
    //     alpha,
    //     centMin,
    //     centMax
    // );

    // // TFile* fout = new TFile("PbPb17p3_rapidityloss_05pCen_a3.root", "RECREATE");
    // TFile* fout = new TFile("./rapidityLossFromTree/" + outName, "RECREATE");

    TRandom3 *rnd = new TRandom3(0); 

    // 3. 参数设置
    int Read_TotNevents = 1000000;

    TTree* t = (TTree*)fin->Get("t");
    if (!t) {
        std::cout << "❌ Cannot find TTree 't'" << std::endl;
        return;
    }
    cout << "Reading Tree" <<endl;

    TTree* tout = new TTree("t", "Rapidity loss tree");

    int    evt, Ncoll, Npart;
    bool   isProton,isProjectile;
    double y_init, y_final, dy, b;
    //double alpha;
    
    t->SetBranchAddress("evt", &evt);
    t->SetBranchAddress("Ncoll", &Ncoll);
    t->SetBranchAddress("Npart", &Npart);
    t->SetBranchAddress("b", &b);
    t->SetBranchAddress("isProton", &isProton);
    t->SetBranchAddress("isProjectile", &isProjectile);
    t->SetBranchAddress("y_final", &y_final);
    //t->SetBranchAddress("alpha", &alpha);
    t->SetBranchAddress("dy", &dy);

    TFile *fout = new TFile("qaRapidityTree.root","RECREATE");


    TH1F *h_b =
        new TH1F("h_b",";b [fm]",100,0,20);

    TH1F *h_Npart =
        new  TH1F("h_Npart",";N_{part}",100,0,400);

    TH1F *h_Ncoll =
        new  TH1F("h_Ncoll",";N_{coll}",50,0,50);
    TH1F *h_NcollP =
        new TH1F("h_NcollP",";N_{coll}",50,0,50);
    TH1F *h_NcollN =
        new TH1F("h_NcollN",";N_{coll}",50,0,50);

    TH2F *h_Ncoll_b =
        new TH2F("h_Ncoll_b",
                 ";b;N_{coll}",
                 100,0,20,
                 50,0,50);

    TH1F *h_dNdy =
        new TH1F("h_dNdy",";y;counts",
                 10000,-10,10);

    TH1F *h_dNdyA =
        new TH1F("h_dNdyA",";y;counts",
                 10000,-10,10);

    TH1F *h_dNdyB =
        new TH1F("h_dNdyB",";y;counts",
                 10000,-10,10);

    TH1F *h_dNdyAP =
        new TH1F("h_dNdyAP",";y;counts",
                 10000,-10,10);

    TH1F *h_dNdyAN =
        new TH1F("h_dNdyAN",";y;counts",
                 10000,-10,10);    
                 
    TH1F *h_dy =
        new TH1F("h_dy",";#Delta y",200,-10,10);

    TProfile *p_dy_Ncoll =
        new TProfile("p_dy_Ncoll",
                     ";N_{coll};<#Delta y>",
                     30,0,30);

    Long64_t nentries = t->GetEntries();

    int current_evt = -1;
    int baryon_count_evt = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        // ---- event change ----
        if (evt != current_evt) {
            if (current_evt >= 0) {
                //h_baryon_per_evt->Fill(baryon_count_evt);
            }
            current_evt = evt;
            baryon_count_evt = 0;

            h_b->Fill(b);
            if (isProton) {
                h_NcollP->Fill(Ncoll);
            } else {
                h_NcollN->Fill(Ncoll);
            }

            h_Ncoll->Fill(Ncoll);
            h_Npart->Fill(Npart);
            h_Ncoll_b->Fill(b, Ncoll);
        }

        baryon_count_evt++;


        // ---- Level-2 ----
        h_dNdy->Fill(y_final);
        //h_y_isProjectile->Fill(y_final, isProjectile);

        if (isProjectile) {
            h_dNdyA->Fill(y_final);

            if (isProton) {
                h_dNdyAP->Fill(y_final);
            } else {
                h_dNdyAN->Fill(y_final);
            }
        } else {
            h_dNdyB->Fill(y_final);
        }

        // ---- Level-3 ----
        h_dy->Fill(dy);

        if (isProjectile) {
            //h_dy_proj->Fill(dy);
            p_dy_Ncoll->Fill(Ncoll, dy);
        }
    }

    // last event
    //h_baryon_per_evt->Fill(baryon_count_evt);

    // Normalization

    h_Ncoll->Scale( 1.0 / (h_Ncoll->Integral()));
    h_NcollP->Scale( 1.0 / (h_NcollP->Integral()));
    h_NcollN->Scale( 1.0 / (h_NcollN->Integral()));

    int NeventsSaved = h_Npart->GetEntries();
    double Npart_mid = h_Npart->GetMean();
    h_dNdy  ->Scale(1.0/NeventsSaved, "width");
    h_dNdyA ->Scale(1.0/NeventsSaved, "width");
    h_dNdyB ->Scale(1.0/NeventsSaved, "width");
    h_dNdyAP->Scale(1.0/NeventsSaved, "width");
    h_dNdyAN->Scale(1.0/NeventsSaved, "width");


    double y_hi = 0.1;
    double y_lo = -y_hi;
    double err = 0.0;
    
    int bin_lo = h_dNdy->FindBin(y_lo);
    int bin_hi = h_dNdy->FindBin(y_hi);
    double integral = h_dNdy ->IntegralAndError(bin_lo,    bin_hi, err, "width");
    
    double dy_range = (bin_hi - bin_lo)*h_dNdy->GetBinWidth(0);
    double avg_dNdy = integral / dy_range;
    double avg_err  = err      / dy_range;

    cout << "<dN/dy>_{|y|<" << y_hi << "} = "
         << avg_dNdy << " ± " << avg_err << endl;

double sigma_dNdy   = err; // 你已经算好的 avg_dNdy 误差
double sigma_Npart  = h_Npart->GetMeanError();

double Q = avg_dNdy * 2.0 / Npart_mid;

double sigma_Q = std::sqrt(
    std::pow(2.0 / Npart_mid * sigma_dNdy, 2) +
    std::pow(2.0 * avg_dNdy / (Npart_mid * Npart_mid) * sigma_Npart, 2)
);

    cout << "(2/Npart)<dN/dy>_{|y|<" << y_hi << "} = "
         << Q << " ± " << sigma_Q << endl;

    // ===============================
    // Output
    // ===============================

    fout->Write();
    fout->Close();

}
