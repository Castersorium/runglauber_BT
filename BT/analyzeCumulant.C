#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>
#include <map>

void analyzeCumulant() {
    // 1. 打开 ROOT 文件
    TFile* f = TFile::Open("tree_AuAu_200_alpha2.00_cent0_5.root");
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

    // 3. 绑定变量
    int    evt, Ncoll;
    bool   isProton;
    double y_final, alpha;

    t->SetBranchAddress("evt", &evt);
    t->SetBranchAddress("Ncoll", &Ncoll);
    t->SetBranchAddress("isProton", &isProton);
    t->SetBranchAddress("y_final", &y_final);
    t->SetBranchAddress("alpha", &alpha);

    // 4. 按 event 聚合 proton 数量
    std::map<int,int> proton_count; // key = evt, value = #protons (y_final>0)

    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);
        if (isProton && y_final > 0) {
            proton_count[evt]++;
        }
    }

    // 5. 计算 cumulants
    double sum1 = 0; // <N>
    double sum2 = 0; // <(N-<N>)^2>
    double sum3 = 0; // <(N-<N>)^3>

    for (auto &kv : proton_count) {
        double N = kv.second;
        sum1 += N;
    }

    double meanN = sum1 / proton_count.size();

    for (auto &kv : proton_count) {
        double x = kv.second;
        sum2 += (x - meanN)*(x - meanN);
        sum3 += (x - meanN)*(x - meanN)*(x - meanN);
    }

    double C2 = sum2 / proton_count.size();
    double C3 = sum3 / proton_count.size();

    std::cout << "Event count: " << proton_count.size() << std::endl;
    std::cout << "y_final > 0 proton cumulants:" << std::endl;
    std::cout << "<N>  = " << meanN << std::endl;
    std::cout << "C2   = " << C2 << std::endl;
    std::cout << "C3   = " << C3 << std::endl;

    // 6. 可选：生成 event-by-event multiplicity histogram
    TH1D* hN = new TH1D("hN", "Event-by-event y_final>0 proton multiplicity", 20, 0, 20);
    for (auto &kv : proton_count) {
        hN->Fill(kv.second);
    }
    hN->GetXaxis()->SetTitle("N_proton (y_final>0)");
    hN->GetYaxis()->SetTitle("Number of events");
    hN->Draw();
}
