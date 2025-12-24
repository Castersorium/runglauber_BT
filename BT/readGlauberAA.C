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

#include <stdio.h>
#include <iostream>

#include "EXP_dNdy.h"

TGraphErrors* makeGraphError(
    const std::vector<dNdyPoint>& data,
    const char* title
) {
    int n = data.size();
    auto* gr = new TGraphErrors(n);
    gr->SetName(title);
    gr->SetTitle(title);

    for (int i = 0; i < n; ++i) {
        double y  = data[i].y;
        double dy = 0.5 * (data[i].y_high - data[i].y_low);

        double val = data[i].value;
        double err = std::sqrt(
            data[i].stat_err * data[i].stat_err +
            data[i].sys_err  * data[i].sys_err
        );

        gr->SetPoint(i, y, val);
        gr->SetPointError(i, dy, err);
    }
    return gr;
}

void readGlauberAA() {

    TString Projectile = "Au2rw"; //Au Au2rw  Au197pnHFB14
    TString Target = Projectile;
    TString Energy = "62p4";
    TString System = Projectile+Target+ "_" + Energy;

    double alpha = 3.0; 

    // 1. 打开 TGlauber 输出文件
    TFile *f = TFile::Open(System + "_nucleons_1M.root");
    if (!f || f->IsZombie()) { 
        std::cout << "❌ Cannot open file\n"; 
        return; 
    }
    // TFile* fout = new TFile("PbPb17p3_rapidityloss_05pCen_a3.root", "RECREATE");
    TFile* fout = new TFile(System +"_rapidityloss_05pCen_a3.root", "RECREATE");

    // 3. 参数设置
    int Read_TotNevents = 100000;
    // double y_beam = 5.36; // AuAu 200 GeV beam rapidity
    // double Npart_mid  = 357;
    // double Npart_width = 8;
    // double alpha = 1.0;  
    // double NScale = 1.0;

    double y_beam = 4.2; // AuAu 62.4 GeV beam rapidity
    // double Npart_mid  = 314;//346.5;//+2.8;//314+8;//BRAHMS0~10%
    // double Npart_width = 8;//2.8;//BRAHMS0~10%
    double Npart_mid   = 15.3; //STAR 0~5%
    double Npart_width =  2.4; //STAR 0~5%
    //double Npart_cut_low = 346.5//-2.8;//314-8;
 
    double NScale = 1.0;

    // double y_beam = 2.9; // PbPb 17.3 GeV beam rapidity
    // double Npart_mid  = 352;
    // double Npart_width = 12;
    // double alpha = 3.0;  
    // double NScale = 1.0;


    TRandom3 *rnd = new TRandom3(0); 


    TH1F *h_dNdy  = new TH1F("h_dNdy", "Particle rapidity; y; dN/dy"   , 500, -y_beam-12, y_beam+12);

    TH1F *h_dNdDyA   = new TH1F("h_dNdDyA" ,   "h_dNdDy; #Deltay=y-y_{beam};(2/N_{part})dN^{#it{N-#bar{N}}}/d(y-y_{beam}) "   , 500, -y_beam-12, y_beam+12);
    TH1F *h_dNdDyAP  = new TH1F("h_dNdDyAP", "h_dNdDyAP; #Deltay=y-y_{beam};(2/N_{part})dN^{#it{p-#bar{p}}}/d(y-y_{beam}) "   , 500, -y_beam-12, y_beam+12);
    TH1F *h_dNdDyAN  = new TH1F("h_dNdDyAN", "h_dNdDyAN; #Deltay=y-y_{beam};(2/N_{part})dN^{#it{n-#bar{n}}}/d(y-y_{beam}) "   , 500, -y_beam-12, y_beam+12);

    TH1F *h_dNdyA = new TH1F("h_dNdyA", "Projectile rapidity; y; dN/dy", 500, -y_beam-10, y_beam+10);
    TH1F *h_dNdyB = new TH1F("h_dNdyB", "Target rapidity; y; dN/dy"    , 500, -y_beam-10, y_beam+10);

    TH1F *h_dNdyP  = new TH1F("h_dNdyP", "Proton rapidity P; y; dN/dy"   , 500, -y_beam-12, y_beam+12);
    TH1F *h_dNdyN  = new TH1F("h_dNdyN", "Nucleon rapidity N; y; dN/dy"   , 500, -y_beam-12, y_beam+12);
    TH1F *h_b = new TH1F("h_b", "impact parameter from AA; b; ", 100, 0, 20);
    TH1F *h_Ncoll = new TH1F("h_Ncoll", "Ncoll from AA; Ncoll; ", 100, 0,1400);
    TH1F *h_NcollA = new TH1F("h_NcollA", "Ncoll of projectile; NcollA; ", 30, 0, 30);
    TH1F *h_NcollAP = new TH1F("h_NcollAP", "Ncoll of projectile proton;  NcollAP; ", 30, 0, 30);
    TH1F *h_NcollAN = new TH1F("h_NcollAN", "Ncoll of projectile neutron; NcollAN; ", 30, 0, 30);

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

    TH2D *h2_b_Ncoll = new TH2D(
        "h2_b_Ncoll",
        "impact parameter vs Ncoll; Ncoll; b[fm]",
        100, 0,1400,       // X: Ncoll
        100, 0, 20        // Y: b
    );
    TH2D *h2_b_NcollA = new TH2D(
        "h2_b_NcollA",
        "impact parameter vs NcollA; NcollA; b[fm]",
        30,  0, 30,       // X: NcollA
        100, 0, 20        // Y: b
    );
    

    TH1F *h_Npart = new TH1F("h_Npart", "Npart from AA; Npart; ", 100, 0, 400);

    // 读取 ntuple
    TString TreeName = "nt_" + Projectile + "_" + Target;
    TNtuple *nt = (TNtuple*)f->Get(TreeName);
    // TNtuple *nt = (TNtuple*)f->Get("nt_Au197pnHFB14_Au197pnHFB14");
    if (!nt) {
        printf("❌ Error: cannot find nt tree\n");
        return;
    }
    float Npart_tree, Ncoll_tree, b;
    nt->SetBranchAddress("Npart", &Npart_tree);
    nt->SetBranchAddress("Ncoll", &Ncoll_tree);
    nt->SetBranchAddress("B", &b);

    // Loop over event
    for (int evt = 0; evt < Read_TotNevents; evt++) {

        if (evt % 5000 == 0){
            std::cout << "Processing  " << evt << "#th events" << std::endl;
        }
        nt->GetEntry(evt);
        //if(b>3.31) continue; // 0~5% 10.1103/PhysRevC.79.034909
        if((Npart_tree>(Npart_mid+Npart_width)) || (Npart_tree<(Npart_mid-Npart_width))) continue;
        TString arrname = Form("nucleonarray%d", evt);
        TObjArray *arr = (TObjArray*)f->Get(arrname);

        h_b->Fill(b);
        h_Ncoll->Fill(Ncoll_tree);
        h_Npart->Fill(Npart_tree);

        h2_b_Ncoll->Fill(Ncoll_tree,b);

        if (!arr) { std::cout << "❌ Cannot find " << arrname << std::endl; return; }       

        // Loop over nucleons
        for (int i=0;i<arr->GetEntriesFast();i++) {
            TGlauNucleon *nuc = (TGlauNucleon*)arr->At(i);
            if (!nuc) continue;
            
            int ncoll = nuc->GetNColl();
            h2_b_NcollA->Fill(ncoll,b);
            if (ncoll == 0) continue;

            double delta_y_total = 0;
            double y_curr = y_beam;

            // First Collision (p+p)
            double delta_y1 = rnd->Exp(1.0);
            delta_y_total = delta_y_total + delta_y1;

            if(nuc->IsInNucleusA()){
                h_NcollA->Fill(ncoll);        
                // projectile first rapidity
                y_curr = y_curr - delta_y1; 
                
                if(nuc->IsProton()){        
                    h_NcollAP->Fill(ncoll);             
                }
    
                if (nuc->IsNeutron()){              
                    h_NcollAN->Fill(ncoll);                
                }
            }

            if (nuc->IsInNucleusB()){
                h_NcollB->Fill(ncoll);                
                // target first rapidity
                y_curr = -y_curr + delta_y1;           
            }

            //2nd~nth Collision (sequential)
            for (int j = 1; j < ncoll; j++) {
                double delta_yi = rnd->Exp(1.0 / alpha);
                if(nuc->IsInNucleusA()){        
                    // projectile sequential rapidity
                    y_curr = y_curr - delta_yi;              
                }
    
                if (nuc->IsInNucleusB()){              
                    // target sequential rapidity
                    y_curr = y_curr + delta_yi; // Dont need to flip second time              
                }

                delta_y_total = delta_y_total + delta_yi;
            }

            double y_final = y_curr;

            h_dNdy->Fill(y_final);

            if(nuc->IsInNucleusA()){        
                h_dNdyA->Fill(y_final);   
                h_dNdDyA->Fill(y_final-y_beam);        
                if(nuc->IsProton()){        
                    h_dNdDyAP->Fill(y_final-y_beam);             
                }
    
                if (nuc->IsNeutron()){              
                    h_dNdDyAN->Fill(y_final-y_beam);                
                }    
            }

            if (nuc->IsInNucleusB()){              
                h_dNdyB->Fill(y_final);                
            }

            h2_DeltaY_Ncoll->Fill(ncoll, delta_y_total);
            h2_dNdy_Ncoll->Fill(ncoll,y_final);

            if(nuc->IsProton()){        
                h_dNdyP->Fill(y_final);             
            }

            if (nuc->IsNeutron()){              
                h_dNdyN->Fill(y_final);                
            }
        }
        
    }
    
    int NeventsSaved = h_Npart->GetEntries();

    h_dNdy ->Scale( 1.0 / NeventsSaved);
    h_dNdyA->Scale( 1.0 / NeventsSaved);
    h_dNdyB->Scale( 1.0 / NeventsSaved);

    h_dNdy ->Scale(1.0, "width");
    h_dNdyA->Scale(1.0, "width");
    h_dNdyB->Scale(1.0, "width");

    h_dNdDyA->Scale(NScale * 1.0 / NeventsSaved/(Npart_mid/2));
    h_dNdDyA->Scale(1.0, "width");
    h_dNdDyAP->Scale(NScale * 1.0 / NeventsSaved/(Npart_mid/2));
    h_dNdDyAP->Scale(1.0, "width");
    h_dNdDyAN->Scale(NScale * 1.0 / NeventsSaved/(Npart_mid/2));
    h_dNdDyAN->Scale(1.0, "width");

    h_dNdyP->Scale(1.0 / NeventsSaved);
    h_dNdyP->Scale(1.0, "width");
    h_dNdyN->Scale(1.0 / NeventsSaved);
    h_dNdyN->Scale(1.0, "width");

    h_NcollA->Scale( 1.0 / (h_NcollA->Integral()));
    h_NcollAP->Scale( 1.0 / (h_NcollAP->Integral()));
    h_NcollAN->Scale( 1.0 / (h_NcollAN->Integral()));

    // TGraphErrors* gr_NA49_17p3GeV_netB =  makeGraphError(NA49_PbPb_17p3GeV_005_netB_05Npart_Dy,  "NA49_PbPb_17p3GeV_005_netB_05Npart_Dy");


    // TGraphErrors* gr_BRAHMS_62p4GeV_netB =  makeGraphError(BRAHMS_AuAu_62p4GeV_005_netB_05Npart_Dy,  "BRAHMS_AuAu_62p4GeV_005_netB_05Npart_Dy");
    // TGraphErrors* gr_BRAHMS_200GeV_netB  =  makeGraphError(BRAHMS_AuAu_200GeV_010_netB_05Npart_Dy ,  "BRAHMS_AuAu_200GeV_010_netB_05Npart_Dy");

    // TGraphErrors* gr_STAR_62p4GeV_netB  =  makeGraphError(STAR_AuAu_62p4GeV_005_netP_05Npart_Dy ,  "STAR_AuAu_62p4GeV_005_netP_05Npart_Dy");
    // TGraphErrors* gr_STAR_200GeV_netB  =  makeGraphError(STAR_AuAu_200GeV_005_netP_05Npart_Dy ,  "STAR_AuAu_200GeV_005_netP_05Npart_Dy");


    // gr_NA49_17p3GeV_netB->Write();
    // gr_BRAHMS_62p4GeV_netB->Write();
    // gr_BRAHMS_200GeV_netB->Write();
    // gr_STAR_62p4GeV_netB->Write();
    // gr_STAR_200GeV_netB->Write();


    fout->Write();
    fout->Close();


}
