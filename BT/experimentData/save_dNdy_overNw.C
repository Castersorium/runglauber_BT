#include <TFile.h>
#include <TGraphErrors.h>
#include <TString.h>

#include <vector>
#include <cmath>
#include <iostream>

#include "../include/EXP_dNdy.h"

/*------------------------------------------
  工具函数：vector<dNdyPoint> -> TGraphErrors
------------------------------------------*/
TGraphErrors* makeGraphError(
    const std::vector<dNdyPoint>& data,
    const char* name
) {
    int n = data.size();
    auto* gr = new TGraphErrors(n);
    gr->SetName(name);
    gr->SetTitle("");

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
/*------------------------------------------
  Graph / <Nw>
------------------------------------------*/
void normalizeByNw(TGraphErrors* gr, double Nw) {
    int n = gr->GetN();
    for (int i = 0; i < n; ++i) {
        double x, y;
        gr->GetPoint(i, x, y);

        double ex = gr->GetErrorX(i);
        double ey = gr->GetErrorY(i);
        //Nw = 1.0;
        //Nw = 0.5 * Nw;
        gr->SetPoint(i, x, y / Nw);
        gr->SetPointError(i, ex, ey / Nw);
    }
}

void normalizeByHalfNw(TGraphErrors* gr, double Nw) {
    int n = gr->GetN();
    for (int i = 0; i < n; ++i) {
        double x, y;
        gr->GetPoint(i, x, y);

        double ex = gr->GetErrorX(i);
        double ey = gr->GetErrorY(i);
        //Nw = 1.0;
        //Nw =  * Nw;
        gr->SetPoint(i, x, 2 * y / Nw);
        gr->SetPointError(i, ex, 2 * ey / Nw);
    }
}
/*------------------------------------------
  主函数
------------------------------------------*/
void save_dNdy_overNw() {

    TFile* fout = new TFile("ExpData_netP_dNdy_overHalfNw.root",
                            "RECREATE");

    /* ---------- <Nw> values ---------- */

    const double Nw_NA49_1999 = 352;

    const double Nw_C0 = 357.0;
    const double Nw_C1 = 288.0;
    const double Nw_C2 = 211.0;
    const double Nw_C3 = 146.0;
    const double Nw_C4 = 85.0;

    const double Nw_BRAHMS_62p4 = 314;
    const double Nw_BRAHMS_200 = 357;

    const double Nw_STAR_62p4 = 346.5;
    const double Nw_STAR_130 = 344.3;
    const double Nw_STAR_200 = 350.6;

    /* ---------- NA49 Centrality C0–C5 ---------- */

    /* ---------- 1999 ---------- */
    auto* gr_NA49_1999_17p3 = makeGraphError(
        NA491999_PbPb_17p3GeV_005_netP_005Cent_dNdy,
        "NA49_1999_17p3GeV_netP_dNdy_overNw"
    );
    normalizeByHalfNw(gr_NA49_1999_17p3, Nw_NA49_1999);

    /* ---------- C0 ---------- */
    auto* gr_C0 = makeGraphError(
        NA492011_PbPb_17p3GeV_C0_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C0"
    );
    normalizeByHalfNw(gr_C0, Nw_C0);

    /* ---------- C1 ---------- */
    auto* gr_C1 = makeGraphError(
        NA492011_PbPb_17p3GeV_C1_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C1"
    );
    normalizeByHalfNw(gr_C1, Nw_C1);

    /* ---------- C2 ---------- */
    auto* gr_C2 = makeGraphError(
        NA492011_PbPb_17p3GeV_C2_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C2"
    );
    normalizeByHalfNw(gr_C2, Nw_C2);

    /* ---------- C3 ---------- */
    auto* gr_C3 = makeGraphError(
        NA492011_PbPb_17p3GeV_C3_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C3"
    );
    normalizeByHalfNw(gr_C3, Nw_C3);

    /* ---------- C4 ---------- */
    auto* gr_C4 = makeGraphError(
        NA492011_PbPb_17p3GeV_C4_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C4"
    );
    normalizeByHalfNw(gr_C4, Nw_C4);


    /* ---------- BRAHMS Au+Au 62.4 GeV, 0–10%, ---------- */
    auto* gr_BRAHMS_62p4 = makeGraphError(
        BRAHMS_AuAu_62p4GeV_010_netP,
        "BRAHMS_AuAu_62p4GeV_010_netP_dNdy_overNw"
    );
    normalizeByHalfNw(gr_BRAHMS_62p4, Nw_BRAHMS_62p4);

    /* ---------- BRAHMS Au+Au 200 GeV, 0–5% ---------- */
    auto* gr_BRAHMS_200 = makeGraphError(
        BRAHMS_AuAu_200GeV_005_netP,
        "BRAHMS_AuAu_200GeV_005_netP_dNdy_overNw"
    );
    normalizeByHalfNw(gr_BRAHMS_200, Nw_BRAHMS_62p4);

    /* ---------- STAR Au+Au 62.4 GeV, 0–5%, ---------- */
    auto* gr_STAR_62p4 = makeGraphError(
        STAR_AuAu_62p4GeV_005_netP,
        "STAR_AuAu_62p4GeV_005_netP_dNdy_overNw"
    );
    normalizeByHalfNw(gr_STAR_62p4, Nw_STAR_62p4);

    /* ---------- STAR Au+Au 130 GeV, 0–6%, ---------- */
    auto* gr_STAR_130 = makeGraphError(
        STAR_AuAu_130GeV_006_netP,
        "STAR_AuAu_130GeV_006_netP_dNdy_overNw"
    );
    normalizeByHalfNw(gr_STAR_130, Nw_STAR_130);

    /* ---------- STAR Au+Au 200 GeV, 0–5% ---------- */
    auto* gr_STAR_200 = makeGraphError(
        STAR_AuAu_200GeV_005_netP,
        "STAR_AuAu_200GeV_005_netP_dNdy_overNw"
    );
    normalizeByHalfNw(gr_STAR_200, Nw_STAR_200);



    /* ---------- 写入文件 ---------- */
    
    gr_NA49_1999_17p3->Write();
    
    gr_C0->Write();
    gr_C1->Write();
    gr_C2->Write();
    gr_C3->Write();
    gr_C4->Write();

    gr_BRAHMS_62p4  ->Write();   
    gr_BRAHMS_200   ->Write();  
    gr_STAR_62p4    ->Write(); 
    gr_STAR_130     ->Write();
    gr_STAR_200     ->Write();

    fout->Write();
    fout->Close();

    std::cout << "Saved net-proton dN/dy over Nw graphs into "
              << fout->GetName() << std::endl;
}
