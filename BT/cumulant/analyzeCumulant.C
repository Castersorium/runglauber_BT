#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>
#include <map>

struct CumulantResult {
    double C1_B, C2_B, C3_B, C4_B;
    double C1_P, C2_P, C3_P, C4_P;
};

CumulantResult ComputeCumulants(const std::map<int,int>& baryon_count) {
    CumulantResult res{};

    int Nevt = baryon_count.size();
    if (Nevt == 0) return res;
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

    res.C1_B = meanN;
    res.C2_B = sum2 / baryon_count.size();
    res.C3_B = sum3 / baryon_count.size();
    res.C4_B = sum4 / baryon_count.size()- 3*res.C2_B*res.C2_B;

    res.C1_P =   (1.0/2.0) * res.C1_B;
    res.C2_P =   (1.0/4.0) * res.C2_B +    (1.0/4.0) * res.C1_B;
    res.C3_P =   (1.0/8.0) * res.C3_B +  3*(1.0/8.0) * res.C2_B;
    res.C4_P =  (1.0/16.0) * res.C4_B +  3*(1.0/8.0) * res.C3_B + 3.0/16.0 * res.C2_B - (1.0/8.0) * res.C1_B;

    return res;
}

void analyzeCumulant() {
    // 1. 打开 ROOT 文件
    TFile* f = TFile::Open("../rapidityTree/tree_Au197pnHFB14Au197pnHFB14_62p4_STAR_alpha3.00_cent0.0_80.0.root");
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
    Int_t  Mult;
    Int_t  evt;
    Int_t  Npart;
    Float_t b;
    std::vector<int>*   Ncoll = nullptr;
    std::vector<int>*   isProton = nullptr;
    std::vector<int>*   isProjectile = nullptr;
    std::vector<float>* y_final = nullptr;
    std::vector<float>* dy = nullptr;
    double w = 0.5;
    // double w = 100;

    // 设定随机数种子
    TRandom3 *rnd = new TRandom3(0); // 0 = 用时间作为种子

    t->SetBranchAddress("evt", &evt);
    t->SetBranchAddress("nMultiplicityTree", &Mult);
    t->SetBranchAddress("Ncoll", &Ncoll);
    t->SetBranchAddress("Npart", &Npart);
    t->SetBranchAddress("isProton", &isProton);
    t->SetBranchAddress("y_final", &y_final);
    // t->SetBranchAddress("alpha", &alpha);

    // 4. 按 event 聚合 baryon 数量
    std::map<int,int> baryon_count; // key = evt, value = #baryon

    Long64_t nentries = t->GetEntries();
    std::cout << "Entries= " << nentries << std::endl;
    double ycut_low  = 0 - w;
    double ycut_high = 0 + w;

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);
         // 确保每个 event 都存在
    }

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);
        baryon_count.try_emplace(evt, 0);
        for (int j=0;j<Mult;j++){
            if ( y_final->at(j) > ycut_low && y_final->at(j) < ycut_high) {

                // std::cout << "evt= " << evt << std::endl;
                baryon_count[evt]++;
            }
        }
    }
    std::cout << "Event count: " << baryon_count.size() << std::endl;

    // 5. 计算 cumulants
    // double sum1 = 0; // <N>
    // double sum2 = 0; // <(N-<N>)^2>
    // double sum3 = 0; // <(N-<N>)^3>
    // double sum4 = 0; // <(N-<N>)^4>- 3<(N-<N>)^2>^2

    // for (auto &kv : baryon_count) {
    //     double N = kv.second;
    //     sum1 += N;
    // }

    // double meanN = sum1 / baryon_count.size();

    // for (auto &kv : baryon_count) {
    //     double x = kv.second;
    //     sum2 += (x - meanN)*(x - meanN);
    //     sum3 += (x - meanN)*(x - meanN)*(x - meanN);
    //     sum4 += (x - meanN)*(x - meanN)*(x - meanN)*(x - meanN);
    // }

    // double C1_B = meanN;
    // double C2_B = sum2 / baryon_count.size();
    // double C3_B = sum3 / baryon_count.size();
    // double C4_B = sum4 / baryon_count.size()- 3*C2_B*C2_B;

    // double C1_P =   (1.0/2.0) * C1_B;
    // double C2_P =   (1.0/4.0) * C2_B +    (1.0/4.0) * C1_B;
    // double C3_P =   (1.0/8.0) * C3_B +  3*(1.0/8.0) * C2_B;
    // double C4_P =  (1.0/16.0) * C4_B +  3*(1.0/8.0) * C3_B + 3.0/16.0 * C2_B - (1.0/8.0) * C1_B;
    CumulantResult Cumulants = ComputeCumulants(baryon_count);

    // 6. 通过抽样方法计算 cumulants 误差
    std::vector<CumulantResult> Cumulant_Array;
    int N_rum = 1000; // 采样次数
    float Sample_Ratio = 0.5; // 每次采样数占总数的比例（概率，并不严格遵守）
    for(int i=0;i<N_rum;i++) {
        std::map<int,int> baryon_count_New;
        for (int j=0;j<baryon_count.size();j++) {
            if(rnd->Rndm()<Sample_Ratio){
                baryon_count_New.try_emplace(j, baryon_count[j]);
                // baryon_count_New[j] = baryon_count[j];
            }
        }
        Cumulant_Array.emplace_back(ComputeCumulants(baryon_count_New));
    }
    // <N_B>
    double Sum = 0, Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C1_B;
        Sum2 += pow(Cumulant_Array[i].C1_B,2);
    }
    double sigma_N_B = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // C2_B 
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C2_B;
        Sum2 += pow(Cumulant_Array[i].C2_B,2);
    }
    double sigma_C2_B = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // C3_B 
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C3_B;
        Sum2 += pow(Cumulant_Array[i].C3_B,2);
    }
    double sigma_C3_B = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // C4_B 
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C4_B;
        Sum2 += pow(Cumulant_Array[i].C4_B,2);
    }
    double sigma_C4_B = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // B: sigma^2/M
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C2_B/Cumulant_Array[i].C1_B;
        Sum2 += pow(Cumulant_Array[i].C2_B/Cumulant_Array[i].C1_B,2);
    }
    double sigma_B_sigma2 = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // B: C3/C1    
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C3_B/Cumulant_Array[i].C1_B;
        Sum2 += pow(Cumulant_Array[i].C3_B/Cumulant_Array[i].C1_B,2);
    }
    double sigma_B_C3C1 = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // B: S sigma  
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C3_B/Cumulant_Array[i].C2_B;
        Sum2 += pow(Cumulant_Array[i].C3_B/Cumulant_Array[i].C2_B,2);
    }
    double sigma_B_S_sigma = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // B: K sigma^2
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C4_B/Cumulant_Array[i].C2_B;
        Sum2 += pow(Cumulant_Array[i].C4_B/Cumulant_Array[i].C2_B,2);
    }
    double sigma_B_K_sigma2 = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // <N_P>
    Sum = 0, Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C1_P;
        Sum2 += pow(Cumulant_Array[i].C1_P,2);
    }
    double sigma_N_P = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // C2_P 
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C2_P;
        Sum2 += pow(Cumulant_Array[i].C2_P,2);
    }
    double sigma_C2_P = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // C3_P 
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C3_P;
        Sum2 += pow(Cumulant_Array[i].C3_P,2);
    }
    double sigma_C3_P = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // C4_P 
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C4_P;
        Sum2 += pow(Cumulant_Array[i].C4_P,2);
    }
    double sigma_C4_P = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // P: sigma^2/M
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C2_P/Cumulant_Array[i].C1_P;
        Sum2 += pow(Cumulant_Array[i].C2_P/Cumulant_Array[i].C1_P,2);
    }
    double sigma_P_sigma2 = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // P: C3/C1    
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C3_P/Cumulant_Array[i].C1_P;
        Sum2 += pow(Cumulant_Array[i].C3_P/Cumulant_Array[i].C1_P,2);
    }
    double sigma_P_C3C1 = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // P: S sigma  
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C3_P/Cumulant_Array[i].C2_P;
        Sum2 += pow(Cumulant_Array[i].C3_P/Cumulant_Array[i].C2_P,2);
    }
    double sigma_P_S_sigma = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    // P: K sigma^2
    Sum = 0; Sum2 = 0;
    for (int i=0;i<N_rum;i++) {
        Sum  +=     Cumulant_Array[i].C4_P/Cumulant_Array[i].C2_P;
        Sum2 += pow(Cumulant_Array[i].C4_P/Cumulant_Array[i].C2_P,2);
    }
    double sigma_P_K_sigma2 = pow(Sum2/N_rum - (Sum*Sum/(N_rum*N_rum)),0.5);
    
    std::cout << "y_final Baryon range: " << w << "" <<std::endl;
    std::cout << "<N_B>  = "         << Cumulants.C1_B << " ± " <<  sigma_N_B << std::endl;
    std::cout << "C2_B   = "         << Cumulants.C2_B << " ± " << sigma_C2_B << std::endl;
    std::cout << "C3_B   = "         << Cumulants.C3_B << " ± " << sigma_C3_B << std::endl;
    std::cout << "C4_B   = "         << Cumulants.C4_B << " ± " << sigma_C4_B << std::endl;
    std::cout << "B: sigma^2/M   = " << Cumulants.C2_B/Cumulants.C1_B << " ± " << sigma_B_sigma2   << std::endl;
    std::cout << "B: C3/C1       = " << Cumulants.C3_B/Cumulants.C1_B << " ± " << sigma_B_C3C1     << std::endl;
    std::cout << "B: S sigma     = " << Cumulants.C3_B/Cumulants.C2_B << " ± " << sigma_B_S_sigma  << std::endl;
    std::cout << "B: K sigma^2   = " << Cumulants.C4_B/Cumulants.C2_B << " ± " << sigma_B_K_sigma2 << std::endl;


    std::cout << "<N_P>  = "         << Cumulants.C1_P << " ± " <<  sigma_N_P << std::endl;
    std::cout << "C2_P   = "         << Cumulants.C2_P << " ± " << sigma_C2_P << std::endl;
    std::cout << "C3_P   = "         << Cumulants.C3_P << " ± " << sigma_C3_P << std::endl;
    std::cout << "C4_P   = "         << Cumulants.C4_P << " ± " << sigma_C4_P << std::endl;
    std::cout << "P: sigma^2/M   = " << Cumulants.C2_P/Cumulants.C1_P << " ± " << sigma_P_sigma2   << std::endl;
    std::cout << "P: C3/C1       = " << Cumulants.C3_P/Cumulants.C1_P << " ± " << sigma_P_C3C1     << std::endl;
    std::cout << "P: S sigma     = " << Cumulants.C3_P/Cumulants.C2_P << " ± " << sigma_P_S_sigma  << std::endl;
    std::cout << "P: K sigma^2   = " << Cumulants.C4_P/Cumulants.C2_P << " ± " << sigma_P_K_sigma2 << std::endl;

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
