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

    // ---------- Au+Au 62.4 GeV ----------
    {"62p4",  0,  10, 314+8, 314-8},

    // ---------- Au+Au 62.4 GeV ----------
    {"62p4",  0,  5, 312.9, 317.7},
    {"62p4",  5, 10, 266.0, 271.0},
    {"62p4", 10, 20, 200.0, 210.0},
    {"62p4", 70, 80, 15.3-2.4, 15.3+2.4},

    // ---------- Au+Au 200 GeV ----------
    {"200",   0,  5, 350.0, 360.0},
    {"200",   5, 10, 300.0, 310.0},

    // ---------- Pb+Pb 17.3 GeV ----------
    {"17p3",  0,  5, 340.0, 364.0}
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

//root -l -b -q 'saveTree.C("Au","Au","200",3.0,0,5)'
void saveTreeAA(
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
    TFile *fin = TFile::Open("./nucleonGeneratedData/" + System + "_nucleons_1k.root");
    if (!fin || fin->IsZombie()) { 
        std::cout << "Error: Cannot open file\n"; 
        return; 
    }

    TString outName;
    outName.Form(
        "tree_%s%s_%s_alpha%.2f_cent%d_%d.root",
        projectile,
        target,
        energy,
        alpha,
        centMin,
        centMax
    );

    // TFile* fout = new TFile("PbPb17p3_rapidityloss_05pCen_a3.root", "RECREATE");
    TFile* fout = new TFile(outName, "RECREATE");

    TRandom3 *rnd = new TRandom3(0); 

    // 3. 参数设置
    int Read_TotNevents = 100000;

    // 读取 ntuple
    TString TreeName = "nt_" + Projectile + "_" + Target;
    TNtuple *nt = (TNtuple*)fin->Get(TreeName);

    if (!nt) {
        printf("❌ Error: cannot find nt tree\n");
        return;
    }
    float in_Npart, in_Ncoll, in_b;
    nt->SetBranchAddress("Npart", &in_Npart);
    nt->SetBranchAddress("Ncoll", &in_Ncoll);
    nt->SetBranchAddress("B", &in_b);

    TTree* tout = new TTree("t", "Rapidity loss tree");

    int    out_evt;
    int    out_Npart,out_Ncoll;
    double out_b;
    bool   isProton,isProjectile;
    double y_init, y_final, dy;

    tout->Branch("evt",   &out_evt);
    tout->Branch("Npart", &out_Npart);
    tout->Branch("Ncoll", &out_Ncoll);
    tout->Branch("b",     &out_b);

    tout->Branch("isProton", &isProton);
    tout->Branch("isProjectile", &isProjectile);

    tout->Branch("y_init",  &y_init);
    tout->Branch("y_final", &y_final);
    tout->Branch("dy",      &dy);
    tout->Branch("alpha",   &alpha);

    // Loop over event
    for (int evt = 0; evt < Read_TotNevents; evt++) {

        out_evt = evt; 

        if (evt % 5000 == 0){
            std::cout << "Processing  " << evt << "#th events" << std::endl;
        }
        nt->GetEntry(evt);
        //if(b>3.31) continue; // 0~5% 10.1103/PhysRevC.79.034909
        if((in_Npart>npartMax) || (in_Npart<npartMin)) continue;
        out_Npart = in_Npart; out_b = in_b;

        TString arrname = Form("nucleonarray%d", evt);
        TObjArray *arr = (TObjArray*)fin->Get(arrname);

        if (!arr) { std::cout << "❌ Cannot find " << arrname << std::endl; return; }       

        // Loop over nucleons
        for (int i=0;i<arr->GetEntriesFast();i++) {
            TGlauNucleon *nuc = (TGlauNucleon*)arr->At(i);
            if (!nuc) continue;
            
            int ncoll = nuc->GetNColl();


            if (ncoll == 0) continue;
            out_Ncoll = ncoll;

            double delta_y_total = 0;
            double y_curr = y_beam;
            

            // First Collision (p+p)
            double delta_y1 = rnd->Exp(1.0);
            delta_y_total = delta_y_total + delta_y1;

            //Fill Ncoll
            if(nuc->IsInNucleusA()){      
                // projectile first rapidity
                y_init = y_beam;
                y_curr = y_curr - delta_y1; 
            }

            if (nuc->IsInNucleusB()){          
                // target first rapidity
                y_init = -y_beam;
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
                dy = delta_y_total;
            }

            y_final = y_curr;

            isProton=0; isProjectile=0;
            if(nuc->IsInNucleusA()){        
                isProjectile = 1;
            } 

            if(nuc->IsProton()){        
                isProton = 1;     
            }
            tout->Fill();


        }
        
    }

    fout->Write();
    fout->Close();


}
