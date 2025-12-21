#ifndef EXP_DNDY_H
#define EXP_DNDY_H

#include <vector>

struct dNdyPoint {
    double y;        // y_center
    double y_low;    // bin low
    double y_high;   // bin high
    double value;    // dN/dy
    double stat_err; // statistical uncertainty (symmetric)
    double sys_err;  // systematic uncertainty (symmetric, 0 if not provided)
};

// BRAHMS Au+Au 62.4 GeV, 0–10%, PBAR   10.1016/j.physletb.2009.05.049
const std::vector<dNdyPoint> BRAHMS_AuAu_62p4GeV_010_PBAR = {
    {0.00,  -0.10,  0.10, 10.1, 0.3, 0.0},
    {0.65,   0.40,  0.90,  9.0, 0.1, 0.0},
    {2.30,   2.10,  2.50,  2.2, 0.1, 0.0},
    {3.00,   2.92,  3.08,  0.6, 0.1, 0.0}
};
// BRAHMS Au+Au 62.4 GeV, 0–10%, P - PBAR  10.1016/j.physletb.2009.05.049
const std::vector<dNdyPoint> BRAHMS_AuAu_62p4GeV_010_netP = {
    {0.00,  -0.10,  0.10, 10.4, 0.583, 0.0},
    {0.65,   0.40,  0.90, 12.5, 0.224, 0.0},
    {2.30,   2.10,  2.50, 26.3, 0.608, 0.0},
    {3.00,   2.92,  3.08, 26.1, 0.608, 0.0}
};

// BRAHMS Au+Au 200 GeV, 0–5%, P - PBAR  10.1103/PhysRevLett.93.102301
const std::vector<dNdyPoint> BRAHMS_AuAu_200GeV_005_netP = {
    {-0.05, -0.10,  0.00,  6.32, 0.55, 1.00},
    { 0.05,  0.00,  0.10,  7.08, 0.68, 1.00},
    { 0.45,  0.40,  0.50,  7.09, 0.82, 1.00},
    { 0.55,  0.50,  0.60,  7.33, 0.81, 1.00},
    { 0.80,  0.75,  0.85,  7.67, 0.58, 1.00},
    { 0.90,  0.85,  0.95,  7.65, 0.70, 1.00},
    { 1.85,  1.70,  2.00, 10.50, 0.90, 2.80},
    { 2.20,  2.00,  2.40, 13.20, 1.70, 2.80},
    { 2.90,  2.70,  3.10, 12.40, 0.33, 3.20}
};

// STAR Au+Au, |y|<0.1, net-proton, dN/dy/(0.5*Npart)  10.1103/PhysRevC.96.044904

// const std::vector<dNdyPoint> STAR_AuAu_7p7GeV_005_netP_05Npart = {
//     { 0.0, -0.1, 0.1, 0.32344, 0.00040, 0.0360 }
// };

// const std::vector<dNdyPoint> STAR_AuAu_11p5GeV_005_netP_05Npart = {
//     { 0.0, -0.1, 0.1, 0.2516, 0.00032, 0.0313 }
// };

// const std::vector<dNdyPoint> STAR_AuAu_19p6GeV_005_netP_05Npart = {
//     { 0.0, -0.1, 0.1, 0.1779, 0.00014, 0.0266 }
// };

// const std::vector<dNdyPoint> STAR_AuAu_27GeV_005_netP_05Npart = {
//     { 0.0, -0.1, 0.1, 0.1499, 0.00014, 0.0226 }
// };

// const std::vector<dNdyPoint> STAR_AuAu_39GeV_005_netP_05Npart = {
//     { 0.0, -0.1, 0.1, 0.1054, 0.00011, 0.0181 }
// };

const std::vector<dNdyPoint> STAR_AuAu_7p7GeV_005_netP = {
    {0.0, -0.1, 0.1, 54.51, 6.100, 0.0}
};

const std::vector<dNdyPoint> STAR_AuAu_11p5GeV_005_netP = {
    {0.0, -0.1, 0.1, 42.5, 5.304, 0.0}
};

const std::vector<dNdyPoint> STAR_AuAu_19p6GeV_005_netP = {
    {0.0, -0.1, 0.1, 30.0, 4.528, 0.0}
};

const std::vector<dNdyPoint> STAR_AuAu_27GeV_005_netP = {
    {0.0, -0.1, 0.1, 25.7, 3.864, 0.0}
};

const std::vector<dNdyPoint> STAR_AuAu_39GeV_005_netP = {
    {0.0, -0.1, 0.1, 18.0, 3.068, 0.0}
};


// STAR Au+Au, |y|<0.1, dN/dy, 10.1103/PhysRevC.79.034909
// errors are the quadratic sum of statistical and systematic uncertainties, and are dominated by the latter
const std::vector<dNdyPoint> STAR_AuAu_62p4GeV_005_netP = {
    { 0.0, -0.1, 0.1, 15.4, 0.000, 2.1 }
};

const std::vector<dNdyPoint> STAR_AuAu_130GeV_005_netP = {
    { 0.0, -0.1, 0.1, 8.24, 0.000, 0.93}
};

const std::vector<dNdyPoint> STAR_AuAu_200GeV_005_netP = {
    { 0.0, -0.1, 0.1, 8.0, 0.000, 1.8 }
};



#endif // EXP_DNDY_H