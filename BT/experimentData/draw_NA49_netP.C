#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#include <vector>
#include <cmath>
#include <iostream>

#include "../include/EXP_dNdy.h"

/*------------------------------------------
  工具函数：生成 TGraphErrors
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
  工具函数：统一美化 Graph
------------------------------------------*/
void styleGraph(TGraphErrors* gr, int color, int marker) {
    gr->SetMarkerStyle(marker);
    gr->SetMarkerSize(2.0);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);
}

/*------------------------------------------
  工具函数：x -> -x（rapidity 左右镜像）
------------------------------------------*/
void mirrorGraphX(TGraphErrors* gr) {
    int n = gr->GetN();
    for (int i = 0; i < n; ++i) {
        double x, y;
        gr->GetPoint(i, x, y);
        gr->SetPoint(i, -x, y);

        // x 误差（dy）不变
        double ex = gr->GetErrorX(i);
        double ey = gr->GetErrorY(i);
        gr->SetPointError(i, ex, ey);
    }
}

/*------------------------------------------
  主函数
------------------------------------------*/
void drawExpData_dNdy() {

    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    TFile* fout = new TFile("ExpData_dNdy.root", "RECREATE");

    /* ---------- 生成 Graph ---------- */
    auto* gr_NA49_17p3_1999 = makeGraphError(
        NA491999_PbPb_17p3GeV_005_netP_005Cent_dNdy,
        "NA49_17p3GeV_1999"
    );

    auto* gr_NA49_17p3_2011 = makeGraphError(
        NA492011_PbPb_17p3GeV_005_netP_005Cent_dNdy,
        "NA49_17p3GeV_2011"
    );

    auto* gr_NA49_17p3_1999_mirror =
    (TGraphErrors*) gr_NA49_17p3_1999->Clone("NA49_17p3GeV_1999_mirror");

    auto* gr_NA49_17p3_2011_mirror =
    (TGraphErrors*) gr_NA49_17p3_2011->Clone("NA49_17p3GeV_2011_mirror");
    mirrorGraphX(gr_NA49_17p3_1999_mirror);
    mirrorGraphX(gr_NA49_17p3_2011_mirror);



    /* ---------- 设置样式（颜色 + 形状） ---------- */
    styleGraph(gr_NA49_17p3_1999,  kBlack,    21); // ●
    styleGraph(gr_NA49_17p3_2011,    kRed,    20); // ■

    styleGraph(gr_NA49_17p3_1999_mirror, kBlack, 25); // open square  
    styleGraph(gr_NA49_17p3_2011_mirror, kRed,   24); // ○  open circle    


    /* ---------- Canvas ---------- */
    auto* c1 = new TCanvas("c1", "net-proton dN/dy", 1200, 800);
    c1->SetTicks();
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.12);

    /* ---------- Frame（控制坐标轴） ---------- */
    TH1F* frame = c1->DrawFrame(
        -3, 0.0,   // x_min, y_min
         3, 60    // x_max, y_max（按你的数据可再调）
    );

    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("y");
    frame->GetYaxis()->SetTitle("dN/dy");

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);

    /* ---------- 画图 ---------- */
    gr_NA49_17p3_1999 ->Draw("P SAME");
    gr_NA49_17p3_2011 ->Draw("P SAME");
    gr_NA49_17p3_1999_mirror ->Draw("P SAME");
    gr_NA49_17p3_2011_mirror ->Draw("P SAME");


    /* ---------- Legend ---------- */
    auto* leg = new TLegend(0.15, 0.70, 0.58, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    leg->AddEntry(gr_NA49_17p3_1999,  "NA49(1999) Pb+Pb 17.3 GeV", "p");
    leg->AddEntry(gr_NA49_17p3_2011,  "NA49(2011) Pb+Pb 17.3 GeV", "p");

    leg->Draw();

    /* ---------- 写文件 ---------- */
    gr_NA49_17p3_1999->Write();
    gr_NA49_17p3_2011->Write();
    gr_NA49_17p3_1999_mirror->Write();
    gr_NA49_17p3_2011_mirror->Write();
    //frame->Write();
    c1->Write();

    fout->Write();
    fout->Close();
}
