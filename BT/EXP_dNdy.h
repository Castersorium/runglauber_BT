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

// BRAHMS net-baryon DN/D(YRAP-YRAP(BEAM))/(NPART/2) https://doi.org/10.17182/hepdata.89449.v1

const std::vector<dNdyPoint> BRAHMS_AuAu_200GeV_010_netB_05Npart_Dy = {
    { -5.412, -5.462, -5.362, 0.032, 0.014, 0.0 },
    { -5.312, -5.362, -5.262, 0.044, 0.015, 0.0 },
    { -4.912, -4.962, -4.862, 0.053, 0.016, 0.0 },
    { -4.812, -4.862, -4.762, 0.058, 0.016, 0.0 },
    { -4.562, -4.612, -4.512, 0.066, 0.015, 0.0 },
    { -4.462, -4.512, -4.412, 0.067, 0.016, 0.0 },
    { -3.512, -3.662, -3.362, 0.108, 0.036, 0.0 },
    { -3.162, -3.362, -2.962, 0.142, 0.040, 0.0 },
    { -2.462, -2.662, -2.262, 0.135, 0.040, 0.0 }
};

const std::vector<dNdyPoint> BRAHMS_AuAu_62p4GeV_005_netB_05Npart_Dy = {
    { -4.191, -4.291, -4.091, 0.066, 0.021, 0.0 },
    { -3.541, -3.791, -3.291, 0.118, 0.026, 0.0 },
    { -1.891, -2.091, -1.691, 0.338, 0.053, 0.0 },
    { -1.191, -1.271, -1.111, 0.340, 0.053, 0.0 }
};


// STAR Au+Au, |y|<0.1, net-proton, dN/dy/(0.5*Npart)  10.1103/PhysRevC.96.044904

const std::vector<dNdyPoint> STAR_AuAu_7p7GeV_005_netP_05Npart = {
    { 0.0, -0.1, 0.1, 0.32344, 0.00040, 0.0360 }
};

const std::vector<dNdyPoint> STAR_AuAu_11p5GeV_005_netP_05Npart = {
    { 0.0, -0.1, 0.1, 0.2516, 0.00032, 0.0313 }
};

const std::vector<dNdyPoint> STAR_AuAu_19p6GeV_005_netP_05Npart = {
    { 0.0, -0.1, 0.1, 0.1779, 0.00014, 0.0266 }
};

const std::vector<dNdyPoint> STAR_AuAu_27GeV_005_netP_05Npart = {
    { 0.0, -0.1, 0.1, 0.1499, 0.00014, 0.0226 }
};

const std::vector<dNdyPoint> STAR_AuAu_39GeV_005_netP_05Npart = {
    { 0.0, -0.1, 0.1, 0.1054, 0.00011, 0.0181 }
};

const std::vector<dNdyPoint> STAR_AuAu_7p7GeV_005_netP_05Npart_Dy = {
    { -2.091, -2.191, -1.991, 0.32344, 0.00040, 0.0360 }
};

const std::vector<dNdyPoint> STAR_AuAu_11p5GeV_005_netP_05Npart_Dy = {
    { -2.449, -2.549, -2.349, 0.2516, 0.00032, 0.0313 }
};

const std::vector<dNdyPoint> STAR_AuAu_19p6GeV_005_netP_05Npart_Dy = {
    { -3.037, -3.137, -2.937, 0.1779, 0.00014, 0.0266 }
};

const std::vector<dNdyPoint> STAR_AuAu_27GeV_005_netP_05Npart_Dy = {
    { -3.361, -3.461, -3.261, 0.1499, 0.00014, 0.0226 }
};

const std::vector<dNdyPoint> STAR_AuAu_39GeV_005_netP_05Npart_Dy = {
    { -3.734, -3.834, -3.634, 0.1054, 0.00011, 0.0181 }
};


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
    { 0.0, -0.1, 0.1, 15.4, 2.1, 0.000}
}; //npart = 350 \pm  2.4

const std::vector<dNdyPoint> STAR_AuAu_130GeV_006_netP = {
    { 0.0, -0.1, 0.1, 8.24,  0.93,0.000}
}; //npart = 344.3 \pm  3.1

const std::vector<dNdyPoint> STAR_AuAu_200GeV_005_netP = {
    { 0.0, -0.1, 0.1, 8.0, 1.8, 0.000 }
}; //npart = 346.5 \pm  2.8

// Divide 2 for only projectile
const std::vector<dNdyPoint> STAR_AuAu_62p4GeV_005_netP_05Npart_Dy = {
    { -4.198, -4.298, -4.098, 0.0880/2, 0.0120, 0.0 }
};

const std::vector<dNdyPoint> STAR_AuAu_130GeV_006_netP_05Npart_Dy = {
    { -4.932, -5.032, -4.832, 0.0478/2, 0.0054, 0.0}
};

const std::vector<dNdyPoint> STAR_AuAu_200GeV_005_netP_05Npart_Dy = {
    { -5.363, -5.463, -5.263, 0.0462/2, 0.0104, 0.0}
};


//NA49 10.1103/PhysRevC.83.014901
const std::vector<dNdyPoint> NA49_PbPb_17p3GeV_005_netP_05Npart_Dy = {

    // y ∈ [-0.52, -0.32]
    { -3.33, -3.43, -3.23, 0.134, 0.008, 0.134*0.07 },

    // y ∈ [-0.32, -0.12]
    { -3.13, -3.23, -3.03, 0.145, 0.007, 0.145*0.07 },

    // y ∈ [-0.12, 0.08]
    { -2.93, -3.03, -2.83, 0.145, 0.006, 0.145*0.07 },

    // y ∈ [0.08, 0.28]
    { -2.73, -2.83, -2.63, 0.148, 0.006, 0.148*0.07 },

    // y ∈ [0.28, 0.48]
    { -2.53, -2.63, -2.43, 0.147, 0.005, 0.147*0.07 },

    // y ∈ [0.48, 0.68]
    { -2.33, -2.43, -2.23, 0.154, 0.005, 0.154*0.07 },

    // y ∈ [0.68, 0.88]
    { -2.13, -2.23, -2.03, 0.175, 0.005, 0.175*0.07 },

    // y ∈ [0.88, 1.08]
    { -1.93, -2.03, -1.83, 0.164, 0.005, 0.164*0.07 },

    // y ∈ [1.08, 1.28]
    { -1.73, -1.83, -1.63, 0.181, 0.006, 0.181*0.07 },

    // y ∈ [1.28, 1.48]
    { -1.53, -1.63, -1.43, 0.218, 0.019, 0.218*0.07 },

    // y ∈ [1.48, 1.68] (p only, p̄ negl)
    { -1.33, -1.43, -1.23, 0.220, 0.011, 0.220*0.07 }
};

const std::vector<dNdyPoint> NA49_PbPb_17p3GeV_005_netB_05Npart_Dy = {

    { -3.52, -3.62, -3.42, 0.121095, 0.0588912, 0.0 },
    { -3.32, -3.42, -3.22, 0.131217, 0.0478491, 0.0 },
    { -3.12, -3.22, -3.02, 0.160663, 0.0441684, 0.0 },
    { -2.92, -3.02, -2.82, 0.177594, 0.0404877, 0.0 },
    { -2.72, -2.82, -2.62, 0.227283, 0.0358868, 0.0 },
    { -2.52, -2.62, -2.42, 0.271268, 0.0312859, 0.0 },
    { -2.32, -2.42, -2.22, 0.302553, 0.0239245, 0.0 },
    { -2.12, -2.22, -2.02, 0.346630, 0.0211640, 0.0 },
    { -1.92, -2.02, -1.82, 0.360340, 0.0174833, 0.0 },
    { -1.72, -1.82, -1.62, 0.378928, 0.0184035, 0.0 },
    { -1.52, -1.62, -1.42, 0.376167, 0.0220842, 0.0 },
    { -1.32, -1.42, -1.22, 0.350495, 0.0174833, 0.0 },
    { -1.12, -1.22, -1.02, 0.300713, 0.0138026, 0.0 },
    { -0.92, -1.02, -0.82, 0.253232, 0.0165631, 0.0 },
    { -0.72, -0.82, -0.62, 0.189004, 0.0138026, 0.0 },
    { -0.52, -0.62, -0.42, 0.172441, 0.0294456, 0.0 }
};


#endif // EXP_DNDY_H