#include <TString.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TRandom.h>

#include <stdio.h>
#include <iostream>


void readGlauberAA() {

    // 1. 打开 TGlauber 输出文件
    TFile *f = TFile::Open("AuAu200_nucleons_1k.root");
    if (!f || f->IsZombie()) { 
        std::cout << "❌ Cannot open file\n"; 
        return; 
    }

    // 3. 参数设置
    int Nevents = 1000;
    double y_beam = 5.36; // AuAu 200 GeV beam rapidity
    double lambda = 1.0;  // 指数分布参数: <Δy> = lambda

    TFile* fout = new TFile("AA_rapidityloss.root", "RECREATE");

    TH1F *h_dNdy = new TH1F("h_dNdy", "AA net-baryon rapidity; y; dN/dy", 100, -y_beam-10, y_beam+10);
    TH1F *h_b = new TH1F("h_b", "impact parameter from AA; b; ", 100, 0, 20);
    TH1F *h_Ncoll = new TH1F("h_Ncoll", "Ncoll from AA; Ncoll; ", 100, 0,1400);

    TH1F *h_NcollA = new TH1F("h_NcollA", "NcollA from AA; NcollA; ", 30, 0, 30);
    TH1F *h_NcollB = new TH1F("h_NcollB", "NcollB from AA; NcollB; ", 30, 0, 30);

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

            if (nuc->IsInNucleusA()) {
                h_NcollA->Fill(nuc->GetNColl());
            }
            if (nuc->IsInNucleusB()) {
                h_NcollB->Fill(nuc->GetNColl());
            }

            int ncoll = nuc->GetNColl();  // TGlauberMC 提供
            double delta_y = 0;

            for (int j=0;j<ncoll;j++) {
                delta_y += gRandom->Exp(lambda);
            }

            double y_final;

            if (nuc->IsInNucleusB()) y_final = y_beam - delta_y;
            else if (nuc->IsInNucleusA()) y_final = -y_beam + delta_y;
            else y_final = 0; // 万一既不在 A 也不在 B

            h_dNdy->Fill(y_final);


        }

        if (evt % 5000 == 0){
            std::cout << "Processed  " << evt << "#th events" << std::endl;
        }
    }

    // 5. 绘图
    //TCanvas *c1 = new TCanvas("c1","pA net-baryon dN/dy",800,600);
    //
    //c1->Divide(2,1);
    //c1->cd(1);
    //h_Ncoll->Draw();
    //c1->cd(2);
    //h_Npart->Draw();
    //c1->SetGrid();
    //c1->Print("pA_NcollNpart.root");

    fout->Write();
    fout->Close();


}
