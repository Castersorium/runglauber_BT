#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>
#include <map>

void analyzeCumulant() {
    // 1. 打开 ROOT 文件
    TFile* f = TFile::Open("../rapidityTree/tree_Au197pnHFB14Au197pnHFB14_200_STAR_alpha3.00_cent0.0_5.0.root");
    if (!f || f->IsZombie()) {
        std::cout << "❌ Cannot open file!" << std::endl;
        return;
    }

    // 2. 获取 TTree
    TTree* t = (TTree*)f->Get("t");
    if (!t) {
        std::cout << "❌ Cannot find TTree 't'" << std::endl;
        return;
    }
    cout << "Reading Tree" <<endl;

    // 3. 绑定变量
    int    evt, Ncoll;
    bool   isProton;
    double y_final, alpha;
    double w = 5.0;

    t->SetBranchAddress("evt", &evt);
    t->SetBranchAddress("Ncoll", &Ncoll);
    t->SetBranchAddress("isProton", &isProton);
    t->SetBranchAddress("y_final", &y_final);
    t->SetBranchAddress("alpha", &alpha);

    // 4. 按 event 聚合 baryon 数量
    std::map<int,int> baryon_count; // key = evt, value = #baryon

    Long64_t nentries = t->GetEntries();
    std::cout << "Entries= " << nentries << std::endl;
    double ycut_low  = 0 - w;
    double ycut_high = 0 + w;

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);
        baryon_count.try_emplace(evt, 0); // 确保每个 event 都存在
    }

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);
        if ( y_final > ycut_low && y_final < ycut_high) {

            // std::cout << "evt= " << evt << std::endl;
            baryon_count[evt]++;
        }
    }
    std::cout << "Event count: " << baryon_count.size() << std::endl;

    // 5. 计算 cumulants
    double sum1 = 0; // <N>
    double sum2 = 0; // <(N-<N>)^2>
    double sum3 = 0; // <(N-<N>)^3>
    double sum4 = 0; // <(N-<N>)^4>- 3<(N-<N>)^2>^2

    for (auto &kv : baryon_count) {
        double N = kv.second;
        sum1 += N;
    }

    double meanN = sum1 / baryon_count.size();

    for (auto &kv : baryon_count) {
        double x = kv.second;
        sum2 += (x - meanN)*(x - meanN);
        sum3 += (x - meanN)*(x - meanN)*(x - meanN);
        sum4 += (x - meanN)*(x - meanN)*(x - meanN)*(x - meanN);
    }

    double C1_B = meanN;
    double C2_B = sum2 / baryon_count.size();
    double C3_B = sum3 / baryon_count.size();
    double C4_B = sum4 / baryon_count.size()- 3*C2_B*C2_B;

    double C1_P =   (1.0/2.0) * C1_B;
    double C2_P =   (1.0/4.0) * C2_B +    (1.0/4.0) * C1_B;
    double C3_P =   (1.0/8.0) * C3_B +  3*(1.0/8.0) * C2_B;
    double C4_P =  (1.0/16.0) * C4_B +  3*(1.0/8.0) * C3_B + 3.0/16.0 * C2_B - (1.0/8.0) * C1_B;

    
    std::cout << "y_final Baryon range: " << w << "" <<std::endl;
    std::cout << "<N_B>  = "         << C1_B << std::endl;
    std::cout << "C2_B   = "         << C2_B << std::endl;
    std::cout << "C3_B   = "         << C3_B << std::endl;
    std::cout << "C4_B   = "         << C4_B << std::endl;
    std::cout << "B: sigma^2/M   = " << C2_B/C1_B << std::endl;
    std::cout << "B: C3/C1       = " << C3_B/C1_B << std::endl;
    std::cout << "B: S sigma     = " << C3_B/C2_B << std::endl;
    std::cout << "B: K sigma^2   = " << C4_B/C2_B << std::endl;


    std::cout << "<N_P>  = "         << C1_P << std::endl;
    std::cout << "C2_P   = "         << C2_P << std::endl;
    std::cout << "C3_P   = "         << C3_P << std::endl;
    std::cout << "C4_P   = "         << C4_P << std::endl;
    std::cout << "P: sigma^2/M   = " << C2_P/C1_P << std::endl;
    std::cout << "P: C3/C1       = " << C3_P/C1_P << std::endl;
    std::cout << "P: S sigma     = " << C3_P/C2_P << std::endl;
    std::cout << "P: K sigma^2   = " << C4_P/C2_P << std::endl;

    // 6. 可选：生成 event-by-event multiplicity histogram
    TString histName = Form("Baryon Multiplicity for |y| < %.1f", w);
    TH1D* hN = new TH1D("hN", histName, 200, 0, 400);
    for (auto &kv : baryon_count) {
        hN->Fill(kv.second);
        //cout <<kv.second << endl;
    }
    hN->GetXaxis()->SetTitle("N_B");
    hN->GetYaxis()->SetTitle("Number of events");
    hN->Draw("HIST");
}
