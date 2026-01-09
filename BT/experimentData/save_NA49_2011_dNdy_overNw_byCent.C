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

        gr->SetPoint(i, x, y / Nw);
        gr->SetPointError(i, ex, ey / Nw);
    }
}

/*------------------------------------------
  主函数
------------------------------------------*/
void save_NA49_2011_dNdy_overNw_byCent() {

    TFile* fout = new TFile("NA49_2011_PbPb_17p3GeV_netP_dNdy.root",
                            "RECREATE");

                                /* ---------- <Nw> values ---------- */
    const double Nw_C0 = 357.0;
    const double Nw_C1 = 288.0;
    const double Nw_C2 = 211.0;
    const double Nw_C3 = 146.0;
    const double Nw_C4 = 85.0;

    /* ---------- Centrality C0–C5 ---------- */

    /* ---------- C0 ---------- */
    auto* gr_C0 = makeGraphError(
        NA492011_PbPb_17p3GeV_C0_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C0"
    );
    normalizeByNw(gr_C0, Nw_C0);

    /* ---------- C1 ---------- */
    auto* gr_C1 = makeGraphError(
        NA492011_PbPb_17p3GeV_C1_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C1"
    );
    normalizeByNw(gr_C1, Nw_C1);

    /* ---------- C2 ---------- */
    auto* gr_C2 = makeGraphError(
        NA492011_PbPb_17p3GeV_C2_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C2"
    );
    normalizeByNw(gr_C2, Nw_C2);

    /* ---------- C3 ---------- */
    auto* gr_C3 = makeGraphError(
        NA492011_PbPb_17p3GeV_C3_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C3"
    );
    normalizeByNw(gr_C3, Nw_C3);

    /* ---------- C4 ---------- */
    auto* gr_C4 = makeGraphError(
        NA492011_PbPb_17p3GeV_C4_netP_dNdy,
        "NA49_2011_17p3GeV_netP_dNdy_overNw_C4"
    );
    normalizeByNw(gr_C4, Nw_C4);

    /* ---------- 写入文件 ---------- */
    gr_C0->Write();
    gr_C1->Write();
    gr_C2->Write();
    gr_C3->Write();
    gr_C4->Write();


    fout->Write();
    fout->Close();

    std::cout << "Saved NA49 2011 net-proton dN/dy graphs (C0–C5) into "
              << fout->GetName() << std::endl;
}
