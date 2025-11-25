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

#include <stdio.h>
#include <iostream>


void readGlauberAA() {

    // 1. 打开 TGlauber 输出文件
    TFile *f = TFile::Open("AuAu200_nucleons_1M.root");
    if (!f || f->IsZombie()) { 
        std::cout << "❌ Cannot open file\n"; 
        return; 
    }

    // 3. 参数设置
    int Nevents = 50000;
    double y_beam = 5.36; // AuAu 200 GeV beam rapidity
    double lambda = 1.0;  // 指数分布参数: <Δy> = lambda

    TFile* fout = new TFile("AA_rapidityloss.root", "RECREATE");

    TH1F *h_dNdy = new TH1F("h_dNdy", "Particle rapidity; y; dN/dy", 100, -y_beam-12, y_beam+12);

    TH1F *h_dNdyA = new TH1F("h_dNdyA", "Projectile rapidity; y; dN/dy", 100, -y_beam-10, y_beam+10);
    TH1F *h_dNdyB = new TH1F("h_dNdyB", "Target rapidity; y; dN/dy", 100, -y_beam-10, y_beam+10);

    TH1F *h_b = new TH1F("h_b", "impact parameter from AA; b; ", 100, 0, 20);
    TH1F *h_Ncoll = new TH1F("h_Ncoll", "Ncoll from AA; Ncoll; ", 100, 0,1400);
    TH1F *h_NcollA = new TH1F("h_NcollA", "Ncoll of projectile; NcollA; ", 30, 0, 30);
    TH1F *h_NcollB = new TH1F("h_NcollB", "Ncoll of target; NcollB; ", 30, 0, 30);




    TH2D *h2_DeltaY_Ncoll = new TH2D(
        "h2_DeltaY_Ncoll",
        "Rapidity loss vs Ncoll; Ncoll; #Delta y",
        30, 0, 30,       // X: Ncoll
        200, 0, 20       // Y: Δy
    );

    TH1F *h_Npart = new TH1F("h_Npart", "Npart from AA; Npart; ", 100, 0, 400);

    // 读取 ntuple
    TNtuple *nt = (TNtuple*)f->Get("nt_Au_Au");
    if (!nt) {
        printf("❌ Error: cannot find nt_Au_Au\n");
        return;
    }
    float Npart, Ncoll, b;
    nt->SetBranchAddress("Npart", &Npart);
    nt->SetBranchAddress("Ncoll", &Ncoll);
    nt->SetBranchAddress("B", &b);

    // 2. 选择事件
    for (int evt = 0; evt < Nevents; evt++) {

        nt->GetEntry(evt);

        TString arrname = Form("nucleonarray%d", evt);
        TObjArray *arr = (TObjArray*)f->Get(arrname);

        h_b->Fill(b);
        h_Ncoll->Fill(Ncoll);
        h_Npart->Fill(Npart);

        if (!arr) { std::cout << "❌ Cannot find " << arrname << std::endl; return; }       

        // 4. 遍历核子
        for (int i=0;i<arr->GetEntriesFast();i++) {
            TGlauNucleon *nuc = (TGlauNucleon*)arr->At(i);
            if (!nuc) continue;
            
            int ncoll = nuc->GetNColl();
            if (ncoll == 0) continue;
            
            //计算总 rapidity loss = A independent exponential collisions
            double delta_y_total = 0;
            for (int j = 0; j < ncoll; j++) {
                delta_y_total += gRandom->Exp(lambda);  // 单次碰撞 rapidity loss
            }
            double y_final = 0;
            
            if(nuc->IsInNucleusA()){
                h_NcollA->Fill(ncoll);
        
                // projectile final rapidity
                y_final = y_beam - delta_y_total;
                h_dNdyA->Fill(y_final);
                
            }

            if (nuc->IsInNucleusB()){
                h_NcollB->Fill(ncoll);
                
                // target final rapidity
                y_final = -y_beam + delta_y_total;
                h_dNdyB->Fill(y_final);
                
            }

            //// 若没有发生碰撞 → 它是 spectator，不产生 rapidity loss
            //if (ncoll == 0) continue;

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
