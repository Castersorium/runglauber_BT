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
#include <TGraph.h>
#include <TRandom.h>

#include <stdio.h>
#include <iostream>


void readGlauberpA() {

    TFile *f = TFile::Open("pAu200_nucleons_1M.root");
    if (!f || f->IsZombie()) { 
        std::cout << "❌ Cannot open file\n"; 
        return; 
    }

    int Nevents = 5000;
    double y_beam = 5.36; 
    double alpha = 5;  

    TRandom3 *rnd = new TRandom3(0);

    TFile* fout = new TFile("pA_rapidityloss.root", "RECREATE");

    TH1F *h_dNdy = new TH1F("h_dNdy", "Projectile rapidity; y; dN/dy", 200, -15, 15);
    TH1F *h_NcollA = new TH1F("h_NcollA", "Ncoll of projectile; NcollA; ", 30, 0, 30);
    TH1F *h_NcollB = new TH1F("h_NcollB", "Ncoll of target; NcollB; ", 30, 0, 30);

    TH2D *h2_DeltaY_Ncoll = new TH2D(
        "h2_DeltaY_Ncoll",
        "Rapidity loss vs Ncoll; Ncoll; #Delta y",
        30, 0, 30,       // X: Ncoll
        200, 0, 20       // Y: Δy
    );

    TNtuple *nt = (TNtuple*)f->Get("nt_p_Au");
    float Npart, Ncoll, b;
    nt->SetBranchAddress("Npart", &Npart);
    nt->SetBranchAddress("Ncoll", &Ncoll);
    nt->SetBranchAddress("B", &b);

    for (int evt = 0; evt < Nevents; evt++) {

        nt->GetEntry(evt);

        TString arrname = Form("nucleonarray%d", evt);
        TObjArray *arr = (TObjArray*)f->Get(arrname);
        if (!arr) continue;

        // Loop over nucleons
        for (int i=0; i<arr->GetEntriesFast(); i++) {
            TGlauNucleon *nuc = (TGlauNucleon*)arr->At(i);
            if (!nuc) continue;

            int ncoll = nuc->GetNColl();

            // 只处理 projectile nucleon (In pA, A = proton)
            if (nuc->IsInNucleusB()){
                h_NcollB->Fill(ncoll);
                continue;
            }

            h_NcollA->Fill(ncoll);

            if (ncoll == 0) continue; // Ignore spectators

            //计算总 rapidity loss = A independent exponential collisions
            double delta_y_total = 0;

            double y_curr = y_beam;

            // First Collision (p+p)
            double delta_y1 = rnd->Exp(1.0);
            y_curr = y_curr - delta_y1;
            delta_y_total = delta_y_total + delta_y1;

            //2nd~nth Collision (sequential)
            for (int j = 1; j < ncoll; j++) {
                double delta_yi = rnd->Exp(1.0 / alpha);
                y_curr = y_curr - delta_yi;
                delta_y_total = delta_y_total + delta_yi;
            }

            // projectile final rapidity
            double y_final = y_curr;

            // 一个 projectile nucleon 只填一次
            h_dNdy->Fill(y_final);
            h2_DeltaY_Ncoll->Fill(ncoll, delta_y_total);

        }
        
        if (evt % 5000 == 0){
            std::cout << "Processed  " << evt << "#th events" << std::endl;
        }

    }

    fout->Write();
    fout->Close();
}
