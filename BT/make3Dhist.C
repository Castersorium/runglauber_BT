#include <TString.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TROOT.h>


//gSystem->Load("libMathMore");  

void make3Dhist(double target_b = -1.0, double tolerance = 0.1) {

    // 打开文件
    TFile *f = TFile::Open("pAu200_nucleons_1k.root");
    if (!f || f->IsZombie()) {
        printf("❌ Error: cannot open AuAu200_nucleons.root\n");
        return;
    }

    // 读取 ntuple
    TNtuple *nt = (TNtuple*)f->Get("nt_p_Au");
    if (!nt) {
        printf("❌ Error: cannot find nt_Au_Au\n");
        return;
    }

    float Npart, Ncoll, b;
    nt->SetBranchAddress("Npart", &Npart);
    nt->SetBranchAddress("Ncoll", &Ncoll);
    nt->SetBranchAddress("B", &b);

    // 如果用户没指定b，则画第0个事件
    int selectedEvent = -1;
    if (target_b < 0) {
        selectedEvent = 0;
    } else {
        Long64_t nEntries = nt->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            nt->GetEntry(i);
            if (fabs(b - target_b) < tolerance) {
                selectedEvent = i;
                break;
            }
        }
        if (selectedEvent < 0) {
            printf("⚠️  No event found with b ≈ %.2f fm (tolerance=%.2f)\n", target_b, tolerance);
            return;
        }
    }

    printf("✅ Selected event #%d with b = %.2f fm\n", selectedEvent, b);

    // 获取该事件的核子数组
    TString arrname = Form("nucleonarray%d", selectedEvent);
    TObjArray *arr = (TObjArray*)f->Get(arrname);
    if (!arr) {
        printf("❌ Error: cannot find %s\n", arrname.Data());
        return;
    }

    // 定义三组散点
    TGraph2D *gWounded   = new TGraph2D();
    TGraph2D *gSpectator = new TGraph2D();
    //TGraph2D *gOther     = new TGraph2D();

    int iw = 0, is = 0, io = 0;
    for (int i = 0; i < arr->GetEntriesFast(); ++i) {
        TGlauNucleon *nuc = (TGlauNucleon*)arr->At(i);
        if (!nuc) continue;

        double x = nuc->GetX();
        double y = nuc->GetY();
        double z = nuc->GetZ();

        if (nuc->IsWounded()) {
            gWounded->SetPoint(iw++, x, y, z);
        } else if (nuc->IsSpectator()) {
            gSpectator->SetPoint(is++, x, y, z);
        // } else {
        //     gOther->SetPoint(io++, x, y, z);
        }
    }

    // 样式设置
    gWounded->SetMarkerStyle(20);
    gWounded->SetMarkerColor(kRed);

    gSpectator->SetMarkerStyle(24);
    gSpectator->SetMarkerColor(kBlue);

    //gOther->SetMarkerStyle(25);
    //gOther->SetMarkerColor(kGray + 2);

    // 绘图
    TCanvas *c1 = new TCanvas("c1", "Nucleon 3D View", 1000, 800);
    c1->cd();

    // 坐标范围固定为 ±15 fm
    TH3F *frame = new TH3F("frame",
        Form("Au+Au 200 GeV (b=%.2f fm);X [fm];Y [fm];Z [fm]", b),
        10, -15, 15,
        10, -15, 15,
        10, -15, 15);
    frame->SetDirectory(0);
    frame->Draw();

    //gOther->Draw("p0 same");
    gWounded->Draw("p0 same");
    gSpectator->Draw("p0 same");

    auto legend = new TLegend(0.1, 0.85, 0.4, 0.95);
    legend->AddEntry(gWounded, "Wounded", "p");
    legend->AddEntry(gSpectator, "Spectator", "p");
    //legend->AddEntry(gOther, "Other", "p");
    legend->Draw();

    // 保存结果
    TString outname = Form("nucleon_3d_b%.2f.root", b);
    TFile *out = new TFile(outname, "RECREATE");
    gWounded->Write("gWounded");
    gSpectator->Write("gSpectator");
    //gOther->Write("gOther");
    c1->Write("c1");
    out->Close();

    printf("✅ 3D graphs for event #%d (b=%.2f fm) saved to pA200_nucleon_3d.root\n", selectedEvent, b);
}
