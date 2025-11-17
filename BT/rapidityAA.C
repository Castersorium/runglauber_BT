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


void rapidityAA() {

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

   TH1F *h_dNdy = new TH1F("h_dNdy", "AA net-baryon rapidity; y; dN/dy", 100, -y_beam-10, y_beam+10);


    // 2. 选择事件
    for (int evt = 0; evt < Nevents; evt++) {

        TString arrname = Form("nucleonarray%d", evt);
        TObjArray *arr = (TObjArray*)f->Get(arrname);
        if (!arr) { std::cout << "❌ Cannot find " << arrname << std::endl; return; }       

        // 4. 遍历核子
        for (int i=0;i<arr->GetEntriesFast();i++) {
            TGlauNucleon *nuc = (TGlauNucleon*)arr->At(i);
            if (!nuc) continue;

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
    TCanvas *c1 = new TCanvas("c1","AA net-baryon dN/dy",800,600);
    h_dNdy->Draw();
    c1->SetGrid();
    c1->Print("c1_1k.root");
}
