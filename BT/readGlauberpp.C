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


void readGlauberpp() {

    TFile *f = TFile::Open("pp200_nucleons_1M.root");
    if (!f || f->IsZombie()) { 
        std::cout << "❌ Cannot open file\n"; 
        return; 
    }

    int Nevents = 10000;
    double y_beam = 5.36; 
    double lambda = 1.0;       

    TFile* fout = new TFile("pp_rapidityloss.root", "RECREATE");

    TH1F *h_dNdy = new TH1F("h_dNdy", "Projectile rapidity; y; dN/dy", 200, -15, 15);
    TH1F *h_NcollA = new TH1F("h_NcollA", "Ncoll of projectile; NcollA; ", 30, 0, 30);
    TH1F *h_NcollB = new TH1F("h_NcollB", "Ncoll of target; NcollB; ", 30, 0, 30);
    TH1F *h_b = new TH1F("h_b","h_b",20,0,1.5);

    TH2D *h2_DeltaY_Ncoll = new TH2D(
        "h2_DeltaY_Ncoll",
        "Rapidity loss vs Ncoll; Ncoll; #Delta y",
        30, 0, 30,       // X: Ncoll
        200, 0, 20       // Y: Δy
    );

    TNtuple *nt = (TNtuple*)f->Get("nt_p_p");
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
            //if (!nuc->IsInNucleusA()) continue;

            h_NcollA->Fill(ncoll);

            // 若没有发生碰撞 → 它是 spectator，不产生 rapidity loss
            if (ncoll == 0) continue;

            //计算总 rapidity loss = A independent exponential collisions
            double delta_y_total = 0;
            for (int j = 0; j < ncoll; j++) {
                delta_y_total += gRandom->Exp(lambda);  // 单次碰撞 rapidity loss
            }

            // projectile final rapidity
            double y_final = y_beam - delta_y_total;

            // 一个 projectile nucleon 只填一次
            h_dNdy->Fill(y_final);
            h2_DeltaY_Ncoll->Fill(ncoll, delta_y_total);

        }
        
        h_b->Fill(b);
        
        if (evt % 1000 == 0){
            std::cout << "Processed  " << evt << "#th events" << std::endl;
            //cout << b << endl;
        }

    }

    fout->Write();
    fout->Close();
}
