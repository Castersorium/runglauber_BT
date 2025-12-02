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

    int Nevents = 100000;
    double y_beam = 5.36; 
    double alpha = 3.0;  

    TRandom3 *rnd = new TRandom3(0);

    TFile* fout = new TFile("pA_rapidityloss_MB.root", "RECREATE");

    TH1F *h_dNdy = new TH1F("h_dNdy", "Projectile rapidity; y; dN/dy", 200, -15, 15);
    TH1F *h_b = new TH1F("h_b", "impact parameter from pA; b; ", 100, 0, 20);

    TH1F *h_NcollA = new TH1F("h_NcollA", "Ncoll of projectile; NcollA; ", 30, 0, 30);
    TH1F *h_NcollB = new TH1F("h_NcollB", "Ncoll of target; NcollB; ", 30, 0, 30);

    TH2D *h2_DeltaY_Ncoll = new TH2D(
        "h2_DeltaY_Ncoll",
        "Rapidity loss vs Ncoll; Ncoll; #Delta y",
        30, 0, 30,       // X: Ncoll
        200, 0, 20       // Y: Δy
    );

    TH2D *h2_dNdy_Ncoll = new TH2D(
        "h2_dNdy_Ncoll",
        "Projectile rapidity vs Ncoll; Ncoll; dN/dy",
         30,  0, 30,       // X: Ncoll
        200,-15, 15        // Y: dNdy
    );
    
    TH2D *h2_b_NcollA = new TH2D(
        "h2_b_NcollA",
        "impact parameter vs NcollA; NcollA; b[fm]",
         30,  0, 30,       // X: Ncoll
         100, 0, 20        // Y: b
    );
    
    
    TNtuple *nt = (TNtuple*)f->Get("nt_p_Au");
    float Npart_tree, Ncoll_tree, b;
    nt->SetBranchAddress("Npart", &Npart_tree);
    nt->SetBranchAddress("Ncoll", &Ncoll_tree);
    nt->SetBranchAddress("B", &b);

    // Loop over event
    for (int evt = 0; evt < Nevents; evt++) {

        nt->GetEntry(evt);
        //if(b>2.0) continue;
        h_b->Fill(b);

        TString arrname = Form("nucleonarray%d", evt);
        TObjArray *arr = (TObjArray*)f->Get(arrname);
        if (!arr) continue;



        // Loop over nucleons
        for (int i=0; i<arr->GetEntriesFast(); i++) {
            TGlauNucleon *nuc = (TGlauNucleon*)arr->At(i);
            if (!nuc) continue;

            int ncoll = nuc->GetNColl();

            
            
            if (ncoll == 0) continue; // Ignore spectators

            // 只处理 projectile nucleon (In pA, A = proton)
            if (nuc->IsInNucleusB()){
                h_NcollB->Fill(ncoll);
                continue;
            }

            h_NcollA->Fill(ncoll);
            h2_b_NcollA->Fill(ncoll,b);


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


            h_dNdy->Fill(y_final);
            h2_DeltaY_Ncoll->Fill(ncoll, delta_y_total);
            h2_dNdy_Ncoll->Fill(ncoll,y_final);



            // XZ:如果要结果compare with 实验的net-proton，那还要把nuclues的变化标出来，不然也是错的。

        }
        
        if (evt % 5000 == 0){
            std::cout << "Processed  " << evt << "#th events" << std::endl;
        }

    }

    fout->Write();
    fout->Close();
}
