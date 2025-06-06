#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <chrono>
#include <algorithm>
using namespace std;

//==========================================================
// 模拟参数结构体定义
//==========================================================
struct SimulationParameters {
    int simulationType;         // 模拟类型：1=单相流，2=两相流，3=多相流，与论文中油气水三相流动对应
    string phaseType;           // 两相流类型："Oil-Water" 或 "Oil-Gas"，仅在 simulationType=2 时使用
    double wellboreTotalLength; // 井筒总长度 (m)，定义模拟的井深
    int wellboreCells;          // 井筒网格数，用于空间离散
    double P0;                  // 井口压力 (Pa)，边界条件
    double T0;                  // 井口温度 (K)，边界条件
    double g;                   // 重力加速度 (m/s²)
    double frictionFactor;      // 摩擦因子，无量纲，用于压降计算
    double massFlow;            // 质量流率 (kg/s)，井筒内的总流量
    double D;                   // 井筒直径 (m)
    double dTdz;                // 温度梯度 (K/m)，沿井深的变化率
    double rhoHCVapor;          // 烃蒸汽密度 (kg/m³)，气相密度
    double rhoHCLiquid;         // 烃液体密度 (kg/m³)，液相密度
    double rhoWater;            // 水密度 (kg/m³)
    double muG;                 // 气体粘度 (Pa·s)
    double muL;                 // 液体粘度 (Pa·s)
    double xHcTotal;            // 烃类总质量分数
    double xWater_Total;        // 水总质量分数
    int numComponents;          // 组分数目，用于闪蒸计算
    vector<string> compNames;   // 组分名称
    vector<double> compTcrit;   // 组分临界温度 (K)
    vector<double> compPcrit;   // 组分临界压力 (Pa)
    vector<double> compOmega;   // 组分偏心因子，无量纲
    vector<double> compMw;      // 组分分子量 (g/mol)
    double simulationTime;      // 总模拟时间 (s)
    double timeStep;            // 时间步长 (s)
    double theta_gl;            // 界面张力 (N/m)，用于漂移模型
    double theta;               // 井斜角 (弧度)，影响漂移速度
    double heatTransferCoeff;
    vector<double> xHcFrac;     //烃组分质量分数
    double r_i;    // 油管内径 (m)
    double r_o;    // 油管外径 (m)
    double k_i;    // 油管材料导热率 (W/(m·K))
    double h_i;    // 油管内壁对流换热系数 (W/(m^2·K))
    double r_c;    // 套管外径 (m)
    double k_c;    // 套管材料导热率 (W/(m·K))
    double r_m;    // 水泥层外径 (m) 或 井筒外径 (m)
    double k_m;    // 水泥/环空材料导热率 (W/(m·K))
    double h_o;    // 井筒外壁对流换热系数 (W/(m^2·K))
};
// ==========================================================
// 全局无限储层定压定温常量
// ==========================================================
static const double P_res = 54.546e6; // Pa
static const double T_res = 350.35;   // K

//==========================================================
// 漂移模型参数
//==========================================================
const double pi_drift = 3.141592653589793;
const double A_drift = 1.088;  // 分布参数常数
const double B_drift = 0.833;  // 分布参数常数
const double m1_drift = 1.017; // 垂直漂移速度权重
const double m2_drift = 2.303; // 速度修正因子
const double m3_drift = 1.0;   // 雷诺数修正指数
const double N1_drift = 1.981; // 水平漂移速度常数
const double N2_drift = 1.759; // 粘度效应常数
const double N3_drift = 0.574; // 粘度指数
const double N4_drift = 0.477; // 界面张力指数
const double a1_drift = 0.577; // 气含率阈值1
const double a2_drift = 0.769; // 气含率阈值2

// 组分局部结构体，存储单组分的热力学性质
struct ComponentLocal {
    string name;   // 组分名称
    double Tcrit;  // 临界温度 (K)
    double Pcrit;  // 临界压力 (Pa)
    double omega;  // 偏心因子
    double Mw;     // 分子量 (g/mol)
};

// 闪蒸计算结果结构体
struct FlashResult {
    double beta;         // 气相摩尔分数
    vector<double> x;    // 液相摩尔分数
    vector<double> y;    // 气相摩尔分数
    vector<double> K;    // 分配系数
};

//漂移模型结构体
struct DriftResult {
    double gasVelocity;    // 气相速度 (m/s)
    double liquidVelocity; // 液相速度 (m/s)
};

const double R_const_global = 8.314462618; // 通用气体常数 (J/(mol·K))
const double M_PI = 3.1415926;

//==========================================================
// PR EOS及闪蒸部分
//==========================================================

// 计算Wilson K值，用于闪蒸初始估计
// 输入：组分，温度 (K)，压力 (Pa)
// 输出：K值 (气液分配系数)
double WilsonK_local(const ComponentLocal& comps, double T, double P) {
    if (!(P > 0 && T > 0)) {
        cerr << "!!! Bad T/P in WilsonK: P=" << P << " T=" << T << "\n";
        P = max(P, 1e3);
        T = max(T, 200.0);
    }
    return (comps.Pcrit / P) * exp(5.37 * (1 + comps.omega) * (1 - comps.Tcrit / T));
}

// 计算PR状态方程中的吸引力参数 a
// 输入：组分，温度 (K)
// 输出：a值 (Pa·m⁶/mol²)
double calc_a_local(const ComponentLocal& comps, double T) {
    double m = 0.37464 + 1.54226 * comps.omega - 0.26992 * comps.omega * comps.omega; // 偏心因子修正
    double alpha = pow(1 + m * (1 - sqrt(T / comps.Tcrit)), 2); // 温度修正因子
    return 0.45724 * R_const_global * R_const_global * comps.Tcrit * comps.Tcrit / comps.Pcrit * alpha;
}

// 计算PR状态方程中的排斥参数 b
// 输入：组分
// 输出：b值 (m³/mol)
double calc_b_local(const ComponentLocal& comps) {
    return 0.07780 * R_const_global * comps.Tcrit / comps.Pcrit;
}

// 计算混合物吸引力参数 a_mix
// 输入：组分向量，摩尔分数向量，温度 (K)
// 输出：混合 a 值 (Pa·m⁶/mol²)
// 使用交叉规则计算多组分混合物的 a
double mix_a_local(const vector<ComponentLocal>& comps, const vector<double>& moleFrac, double T) {
    double a_mix = 0.0;
    int n = comps.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a_i = calc_a_local(comps[i], T);
            double a_j = calc_a_local(comps[j], T);
            a_mix += moleFrac[i] * moleFrac[j] * sqrt(a_i * a_j); // 几何平均混合规则
        }
    }
    return a_mix;
}

// 计算混合物排斥参数 b_mix
// 输入：组分向量，摩尔分数向量
// 输出：混合 b 值 (m³/mol)
// 使用线性混合规则
double mix_b_local(const vector<ComponentLocal>& comps, const vector<double>& moleFrac) {
    double b_mix = 0.0;
    int n = comps.size();
    for (int i = 0; i < n; i++) {
        b_mix += moleFrac[i] * calc_b_local(comps[i]);
    }
    return b_mix;
}

// 求解PR状态方程的三次方程
// 输入：无量纲参数 A, B
// 输出：所有实根的向量
vector<double> solveCubic_local(double A, double B) {
    double c2 = -(1 - B);
    double c1 = A - 3 * B * B - 2 * B;
    double c0 = -(A * B - B * B - B * B * B);
    double a_over_3 = c2 / 3.0;
    double p = c1 - c2 * c2 / 3.0;
    double q = 2 * pow(c2, 3) / 27.0 - c2 * c1 / 3.0 + c0;
    double disc = q * q / 4.0 + p * p * p / 27.0; // 判别式
    vector<double> roots;
    if (disc > 1e-12) { // 单根
        double sqrt_disc = sqrt(disc);
        double u = cbrt(-q / 2.0 + sqrt_disc);
        double v = cbrt(-q / 2.0 - sqrt_disc);
        double x = u + v;
        double Z = x - a_over_3;
        roots.push_back(Z);
    }
    else if (fabs(disc) < 1e-12) { // 两实根
        double u = cbrt(-q / 2.0);
        double x1 = 2 * u;
        double x2 = -u;
        roots.push_back(x1 - a_over_3);
        roots.push_back(x2 - a_over_3);
        roots.push_back(x2 - a_over_3);
    }
    else { // 三实根
        const double M_PI = 3.1415926;
        double r = sqrt(-p * p * p / 27.0);
        double phi = acos(-q / (2 * r));
        double x1 = 2 * cbrt(r) * cos(phi / 3.0);
        double x2 = 2 * cbrt(r) * cos((phi + 2 * M_PI) / 3.0);
        double x3 = 2 * cbrt(r) * cos((phi + 4 * M_PI) / 3.0);
        roots.push_back(x1 - a_over_3);
        roots.push_back(x2 - a_over_3);
        roots.push_back(x3 - a_over_3);
    }
    return roots;
}

// 从三次方程根中选择合适的压缩因子 Z
// 输入：根向量，是否为气相 (true=气相，false=液相)
// 输出：Z 值
double selectZ_local(const vector<double>& roots, bool isVapor) {
    double Z = roots[0];
    if (isVapor) {
        for (double r : roots) { if (r > Z) Z = r; } // 气相取最大根
    }
    else {
        for (double r : roots) { if (r < Z) Z = r; } // 液相取最小根
    }
    return Z;
}

// 求解压缩因子 Z
// 输入：A, B，是否为气相
// 输出：Z 值
double cubicSolver_local(double A, double B, bool isVapor) {
    vector<double> roots = solveCubic_local(A, B);
    return selectZ_local(roots, isVapor);
}

// 计算逸度系数
// 输入：组分向量，摩尔分数向量，温度 (K)，压力 (Pa)，是否为气相
// 输出：每组分的逸度系数向量
vector<double> computeFugacityCoeffs_local(const vector<ComponentLocal>& comps,
    const vector<double>& moleFrac, double T, double P, bool isVapor) {
    int n = comps.size();
    vector<double> phi(n, 0.0);
    double a_mix = mix_a_local(comps, moleFrac, T);
    double b_mix = mix_b_local(comps, moleFrac);
    double A_val = a_mix * P / (R_const_global * R_const_global * T * T);
    double B_val = b_mix * P / (R_const_global * T);
    double Z = cubicSolver_local(A_val, B_val, isVapor);
    for (int i = 0; i < n; i++) {
        double a_i = calc_a_local(comps[i], T);
        double b_i = calc_b_local(comps[i]);
        double sum_aij = 0.0;
        for (int j = 0; j < n; j++) {
            double a_j = calc_a_local(comps[j], T);
            sum_aij += moleFrac[j] * sqrt(a_i * a_j);
        }
        double ln_phi = b_i / b_mix * (Z - 1) - log(Z - B_val)
            - A_val / (2 * sqrt(2) * B_val) * (2 * sum_aij / a_mix - b_i / b_mix)
            * log((Z + (1 + sqrt(2)) * B_val) / (Z + (1 - sqrt(2)) * B_val));
        phi[i] = exp(ln_phi);
    }
    return phi;
}

// 闪蒸计算主函数
// 输入：组分向量，总摩尔分数向量，温度 (K)，压力 (Pa)
// 输出：闪蒸结果
FlashResult flashCalculation_local(const vector<ComponentLocal>& comps,
    const vector<double>& z0, double T, double P) {
    int n = comps.size();
    vector<double> K(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        K[i] = WilsonK_local(comps[i], T, P);
        //std::cout << "K\[" << i << "] = " << K\[i] << endl;
    }
    double beta = 0.5; // 初始气相分数猜测
    double tol_flash = 1e-8; // 收敛容差
    int maxIter_flash = 100; // 最大迭代次数
    for (int iter = 0; iter < maxIter_flash; iter++) {
        double f_beta = 0.0; // 目标函数
        double df_beta = 0.0; // 目标函数导数
        for (int i = 0; i < n; i++) {
            double denom = 1 + beta * (K[i] - 1);
            f_beta += z0[i] * (K[i] - 1) / denom;
            df_beta -= z0[i] * (K[i] - 1) * (K[i] - 1) / (denom * denom);
        }
        double beta_new = beta - f_beta / df_beta;
        beta_new = max(0.0, min(1.0, beta_new)); // 限制在 \[0,1]
        //std::cout << "Iter " << iter << ": beta = " << beta_new << ", f_beta = " << f_beta << endl;
        if (fabs(beta_new - beta) < tol_flash) {
            beta = beta_new;
            break;
        }
        beta = beta_new;
    }

    vector<double> x(n, 0.0), y(n, 0.0);
    for (int i = 0; i < n; i++) {
        double denom = 1 + beta * (K[i] - 1);
        x[i] = z0[i] / denom; // 液相分数
        y[i] = K[i] * x[i];  // 气相分数
    }
    double sumx = 0, sumy = 0;
    for (int i = 0; i < n; i++) { sumx += x[i]; sumy += y[i]; }
    for (int i = 0; i < n; i++) { x[i] /= sumx; y[i] /= sumy; } // 归一化

    FlashResult res;
    res.beta = beta;
    res.x = x;
    res.y = y;
    res.K = K;
    return res;
}

double calcRhoVapor(const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    double T, double P) {
    FlashResult flash = flashCalculation_local(comps, z0, T, P);
    auto& y = flash.y;
    double a_mix = mix_a_local(comps, y, T);
    double b_mix = mix_b_local(comps, y);
    double A = a_mix * P / (R_const_global * R_const_global * T * T);
    double B = b_mix * P / (R_const_global * T);
    double Z = cubicSolver_local(A, B, /*isVapor=*/true);
    double M_mix = 0.0;
    for (size_t i = 0; i < comps.size(); ++i)
        M_mix += y[i] * (comps[i].Mw * 1e-3);
    return P * M_mix / (Z * R_const_global * T);
}

// —— PR‐EOS 计算液相密度 ——
double calcRhoLiquid(const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    double T, double P) {
    FlashResult flash = flashCalculation_local(comps, z0, T, P);
    auto& x = flash.x;
    double a_mix = mix_a_local(comps, x, T);
    double b_mix = mix_b_local(comps, x);
    double A = a_mix * P / (R_const_global * R_const_global * T * T);
    double B = b_mix * P / (R_const_global * T);
    double Z = cubicSolver_local(A, B, /*isVapor=*/false);
    double M_mix = 0.0;
    for (size_t i = 0; i < comps.size(); ++i)
        M_mix += x[i] * (comps[i].Mw * 1e-3);
    return P * M_mix / (Z * R_const_global * T);
}

// 获取气相摩尔分数 beta
double getFlashBeta_local(const vector<ComponentLocal>& comps, const vector<double>& z0, double T, double P) {
    FlashResult res = flashCalculation_local(comps, z0, T, P);
    return res.beta;
}

// 计算 beta 对压力和温度的导数
// 输入：组分，总摩尔分数，温度，压力，有限差分步长
// 输出：dBeta_dP, dBeta_dT, beta
void getFlashBetaAndDerivatives_local(const vector<ComponentLocal>& comps, const vector<double>& z0,
    double T, double P, double dP_fd, double dT_fd,
    double& dBeta_dP, double& dBeta_dT, double& beta) {
    beta = getFlashBeta_local(comps, z0, T, P);
    double betaP_plus = getFlashBeta_local(comps, z0, T, P + dP_fd);
    double betaP_minus = getFlashBeta_local(comps, z0, T, P - dP_fd);
    dBeta_dP = (betaP_plus - betaP_minus) / (2 * dP_fd);
    double betaT_plus = getFlashBeta_local(comps, z0, T + dT_fd, P);
    double betaT_minus = getFlashBeta_local(comps, z0, T - dT_fd, P);
    dBeta_dT = (betaT_plus - betaT_minus) / (2 * dT_fd);
}

//==========================================================
// 混合密度计算及其导数
//==========================================================
// 计算混合密度
// 输入：组分，总摩尔分数，温度，压力，有限差分步长，模拟参数
// 输出：混合密度 (kg/m³)，同时计算 dRho_dP, dRho_dT
double computeRhoMix_local(
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    double T, double P,
    double dP_fd, double dT_fd,
    double& dRho_dP, double& dRho_dT,
    const SimulationParameters& simParams)
{
    // 1) 获取闪蒸及其导数
    double beta, dBeta_dP, dBeta_dT;
    getFlashBetaAndDerivatives_local(
        comps, z0, T, P,
        dP_fd, dT_fd,
        dBeta_dP, dBeta_dT, beta
    );

    // 2) 计算相密度
    double rho_g = calcRhoVapor(comps, z0, T, P);
    double rho_l = calcRhoLiquid(comps, z0, T, P);

    // 3) 有限差分计算相密度导数
    double rho_g_p = calcRhoVapor(comps, z0, T, P + dP_fd);
    double rho_g_m = calcRhoVapor(comps, z0, T, P - dP_fd);
    double rho_l_p = calcRhoLiquid(comps, z0, T, P + dP_fd);
    double rho_l_m = calcRhoLiquid(comps, z0, T, P - dP_fd);

    double rho_g_Tp = calcRhoVapor(comps, z0, T + dT_fd, P);
    double rho_g_Tm = calcRhoVapor(comps, z0, T - dT_fd, P);
    double rho_l_Tp = calcRhoLiquid(comps, z0, T + dT_fd, P);
    double rho_l_Tm = calcRhoLiquid(comps, z0, T - dT_fd, P);

    double drho_g_dP = (rho_g_p - rho_g_m) / (2.0 * dP_fd);
    double drho_l_dP = (rho_l_p - rho_l_m) / (2.0 * dP_fd);
    double drho_g_dT = (rho_g_Tp - rho_g_Tm) / (2.0 * dT_fd);
    double drho_l_dT = (rho_l_Tp - rho_l_Tm) / (2.0 * dT_fd);

    // 4) 混合密度及其导数
    double rho_mix = beta * rho_g + (1.0 - beta) * rho_l;
    dRho_dP = beta * drho_g_dP
        + (1.0 - beta) * drho_l_dP
        + (rho_g - rho_l) * dBeta_dP;
    dRho_dT = beta * drho_g_dT
        + (1.0 - beta) * drho_l_dT
        + (rho_g - rho_l) * dBeta_dT;

    return rho_mix;
}

// 计算 F 对 P 和 T 的导数
// 输入：同上
// 输出：dF_dP, dF_dT
void compute_dF_dP_dT_local(const vector<ComponentLocal>& comps, const vector<double>& z0,
    double T, double P, double dz, double dP_fd, double dT_fd,
    double& dF_dP, double& dF_dT, const SimulationParameters& simParams) {
    if (simParams.simulationType == 1 || simParams.simulationType == 2) {
        // 单相或两相流（无挥发）：密度不随 P, T 变化，导数为0
        dF_dP = 0.0;
        dF_dT = 0.0;
    }
    else {
        // 多相流：考虑密度随相态变化
        double dRho_dP, dRho_dT;
        double rho_mix = computeRhoMix_local(comps, z0, T, P, dP_fd, dT_fd, dRho_dP, dRho_dT, simParams);
        double A_val = 3.141592653589793 * simParams.D * simParams.D / 4.0;
        double C = simParams.frictionFactor * (simParams.massFlow * simParams.massFlow / (2.0 * simParams.D * A_val * A_val)) * dz;
        dF_dP = dRho_dP * simParams.g * dz - C / (rho_mix * rho_mix) * dRho_dP;
        dF_dT = dRho_dT * simParams.g * dz - C / (rho_mix * rho_mix) * dRho_dT;
    }
}

//==========================================================
// Gauss 消元求解线性系统
//==========================================================
// 解线性方程组 J \* dU = -R
// 输入：雅可比矩阵 J，残差向量 R
// 输出：解向量 dU
// 返回：是否成功

bool solveLinearSystem(
    const vector<vector<double>>& J,
    const vector<double>& R,
    vector<double>& dU)
{
    int n = R.size();
    // 构造增广矩阵 A (n × (n+1))
    vector<vector<double>> A(n, vector<double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = J[i][j];
        }
        A[i][n] = -R[i];
    }
    // 前向消去（带主元）
    for (int k = 0; k < n; ++k) {
        // 找到第 k 列中绝对值最大的主元行
        int piv = k;
        double maxv = fabs(A[k][k]);
        for (int i = k + 1; i < n; ++i) {
            double v = fabs(A[i][k]);
            if (v > maxv) {
                maxv = v;
                piv = i;
            }
        }
        if (maxv < 1e-14) {
            // 矩阵奇异或病态，无法求解
            return false;
        }
        // 交换当前行和主元行
        if (piv != k) swap(A[k], A[piv]);
        // 消元
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j <= n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
        }
    }

    // 回代
    dU.assign(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = A[i][n];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[i][j] * dU[j];
        }
        dU[i] = sum / A[i][i];
    }
    return true;
}

// 计算分布参数 C0
// 输入：模拟参数，混合速度 (m/s)，气含率
// 输出：C0 值
double calculate_C0(const SimulationParameters& simParams, double vm_current, double alpha_g) {
    double rho_l = simParams.rhoHCLiquid;
    double rho_g = simParams.rhoHCVapor;
    double q1 = sqrt(pi_drift * (rho_l - rho_g)) * simParams.D + 3.333;
    double q2 = min(3.87 - 19.105 / q1, 3.2);
    double v_c = pow(simParams.theta_gl * (rho_l - rho_g) / (rho_l * rho_l), 0.25); // 特征速度
    double Ku = max(q2, 0.0);
    double v_sgf = Ku * (rho_l / rho_g) * v_c;
    double q3 = pow((max(alpha_g, alpha_g * fabs(vm_current) / v_sgf) - B_drift) / (1 - B_drift), 2);
    return A_drift / (1 + (A_drift - 1) * min(q3, 1.0));
}

// 计算漂移速度 vd
// 输入：模拟参数，C0，气含率，混合速度
// 输出：漂移速度 (m/s)
double calculate_vd(const SimulationParameters& simParams, double C0, double alpha_g, double vm_current) {
    double rho_l = simParams.rhoHCLiquid;
    double rho_g = simParams.rhoHCVapor;
    double mu_l = simParams.muL;
    double theta = simParams.theta;
    double Ku = max(min(3.87 - 19.105 / (sqrt(simParams.g * (rho_l - rho_g) / simParams.theta_gl)) * simParams.D + 3.333, 3.2), 0.0);
    double K_alpha_g;
    if (alpha_g < a1_drift) {
        K_alpha_g = 1.53 / C0;
    }
    else if (alpha_g < a2_drift) {
        K_alpha_g = 1.53 / C0 + (Ku - 1.53 / C0) * 0.5 * (1 - cos(pi_drift * (alpha_g - a1_drift) / (a2_drift - a1_drift)));
    }
    else {
        K_alpha_g = Ku;
    }
    double vc = pow(simParams.theta_gl * simParams.g * (rho_l - rho_g) / pow(rho_l, 2), 0.25);
    double v_vd = ((1 - alpha_g * C0) * C0 * K_alpha_g * vc) /
        (alpha_g * C0 * sqrt(rho_g / rho_l) + 1 - alpha_g * C0); // 垂直漂移速度
    double N_mu = mu_l / ((rho_l - rho_g) * pow(simParams.D, 1.5) * sqrt(pi_drift)); // 粘度无量纲数
    double N_E = pi_drift * (rho_l - rho_g) * simParams.D * simParams.D / simParams.theta_gl; // 界面张力无量纲数
    double v_hd = sqrt(pi_drift * simParams.D) * (N1_drift - N2_drift * (pow(N_mu, N3_drift) / pow(N_E, N4_drift))) * (alpha_g * (1 - alpha_g)); // 水平漂移速度
    double Rel = vm_current * rho_l * simParams.D / mu_l; // 雷诺数
    return (m1_drift * v_vd * sin(theta)) +
        (1 - 2.0 / (1.0 + exp(50.0 * sin(theta + m2_drift * vm_current)))) * v_hd * cos(theta) * pow((1.0 + 1000.0 / (Rel + 1000.0)), m3_drift);

}

// 漂移模型主函数
// 输入：模拟参数，混合速度，气相质量分数
// 输出：气速和液速
DriftResult driftModelDetailed(
    const SimulationParameters& simParams,
    double mixVel,
    double beta,
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    double T,
    double P
) {
    // —— 1) 本格闪蒸，拿到分相摩尔分数 beta 和 x,y —— 
    FlashResult fl = flashCalculation_local(comps, z0, T, P);

    // —— 2) 计算两相密度和摩尔体积 —— 
    double rho_g = calcRhoVapor(comps, z0, T, P);
    double rho_l = calcRhoLiquid(comps, z0, T, P);
    double M_mix_g = 0.0, M_mix_l = 0.0;
    for (size_t i = 0; i < comps.size(); ++i) {
        M_mix_g += fl.y[i] * (comps[i].Mw * 1e-3);
        M_mix_l += fl.x[i] * (comps[i].Mw * 1e-3);
    }
    double Vm_g = M_mix_g / rho_g;
    double Vm_l = M_mix_l / rho_l;

    // —— 3) mol 分数 beta -> 体积分数 alpha_g —— 
    double alpha_g = (fl.beta * Vm_g) / (fl.beta * Vm_g + (1 - fl.beta) * Vm_l);

    // —— **如果没有气相，直接返回** —— 
    if (alpha_g < 1e-8) {
        // 整个截面都是液相，气速为 0，液速就是混合速
        return { 0.0, mixVel };
    }

    // —— 4) 正常调用漂移模型 —— 
    double C0 = calculate_C0(simParams, mixVel, alpha_g);
    double v_d = calculate_vd(simParams, C0, alpha_g, mixVel);
    double v_gas = C0 * mixVel + v_d;
    double v_liq = (mixVel - alpha_g * v_gas) / (1 - alpha_g);

    return { v_gas, v_liq };
}
inline double rho_liquidFraction_local(
    const vector<ComponentLocal>& comps,
    const vector<double>& z,
    double T,
    double P)
{
    FlashResult flash = flashCalculation_local(comps, z, T, P);
    return 1.0 - flash.beta;
}

// 计算摩尔体积或焓偏差项系数 da_mix/dT
inline double da_mix_dT(const vector<ComponentLocal>& comps,
    const vector<double>& frac,
    double T)
{
    double dT = 1e-3;
    double a_p = mix_a_local(comps, frac, T + dT);
    double a_m = mix_a_local(comps, frac, T - dT);
    return (a_p - a_m) / (2 * dT);
}

// 计算焓偏差（每摩尔）
inline double enthalpyDeparture_mol(const vector<ComponentLocal>& comps,
    const vector<double>& frac,
    double T, double P,
    bool isVapor)
{
    // EOS混合参数
    double a_mix = mix_a_local(comps, frac, T);
    double b_mix = mix_b_local(comps, frac);
    double A = a_mix * P / (R_const_global * R_const_global * T * T);
    double B = b_mix * P / (R_const_global * T);
    double Z = cubicSolver_local(A, B, isVapor);
    // 温度导数
    double dadT = da_mix_dT(comps, frac, T);
    // 焓偏差公式
    double term = (Z + (1 + sqrt(2.0)) * B) / (Z + (1 - sqrt(2.0)) * B);
    double h_dep = R_const_global * T * (Z - 1)
        - (T * dadT - a_mix) / (2 * sqrt(2.0) * b_mix) * log(term);
    return h_dep;
}

// 计算局部内能 (J/kg)
double calcInternalEnergy_local(
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    double T, double P)
{
    // 1) 闪蒸分相
    FlashResult flash = flashCalculation_local(comps, z0, T, P);
    double beta = flash.beta;
    // 2) 相别摩尔分率
    auto& y = flash.y;
    auto& x = flash.x;
    // 3) 相密度
    double rho_g = calcRhoVapor(comps, z0, T, P);
    double rho_l = calcRhoLiquid(comps, z0, T, P);
    // 4) 焓偏差 (J/mol)
    double h_dep_g = enthalpyDeparture_mol(comps, y, T, P, true);
    double h_dep_l = enthalpyDeparture_mol(comps, x, T, P, false);
    // 5) 摩尔质量 (kg/mol)
    double M_g = 0.0, M_l = 0.0;
    for (size_t i = 0; i < comps.size(); ++i) {
        M_g += y[i] * (comps[i].Mw * 1e-3);
        M_l += x[i] * (comps[i].Mw * 1e-3);
    }
    // 6) 质量基内能 u = (h_dep/mol - P/ρ)/M_mix
    double u_g = (h_dep_g - P / rho_g) / M_g;
    double u_l = (h_dep_l - P / rho_l) / M_l;
    // 7) 混合内能
    return beta * u_g + (1 - beta) * u_l;
}

// 计算焓通量 (J/m^3)
double calcEnthalpyFlux_local(
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    double T, double P)
{
    // 1) 闪蒸
    FlashResult flash = flashCalculation_local(comps, z0, T, P);
    double beta = flash.beta;
    auto& y = flash.y;
    auto& x = flash.x;
    // 2) 相密度
    double rho_g = calcRhoVapor(comps, z0, T, P);
    double rho_l = calcRhoLiquid(comps, z0, T, P);
    // 3) 相别焓偏差 J/mol
    double h_dep_g = enthalpyDeparture_mol(comps, y, T, P, true);
    double h_dep_l = enthalpyDeparture_mol(comps, x, T, P, false);
    // 4) 摩尔质量 (kg/mol)
    double M_g = 0.0, M_l = 0.0;
    for (size_t i = 0; i < comps.size(); ++i) {
        M_g += y[i] * (comps[i].Mw * 1e-3);
        M_l += x[i] * (comps[i].Mw * 1e-3);
    }
    // 5) 质量基焓 (J/kg)
    double h_g = h_dep_g / M_g + P / rho_g;
    double h_l = h_dep_l / M_l + P / rho_l;
    // 6) 焓通量 = β·ρ_g·h_g + (1-β)·ρ_l·h_l (J/m3)
    return beta * rho_g * h_g + (1 - beta) * rho_l * h_l;
}

// 计算单元整体内能（体积基，单位 J）
double computeCellEnergy_local(
    int cell,
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    const vector<double>& U,     // 传入 Uc（全局解向量）
    const vector<double>& Uf,    // 传入 Uf（face 速度）
    const SimulationParameters& simParams,
    double dz
) {
    int nc = simParams.numComponents;
    int numCvar = nc + 2;
    auto idxC = [&](int i, int eq) { return i * numCvar + eq; };

    // 1) 单元体积 (m^3)
    double Vc = M_PI * simParams.D * simParams.D / 4.0 * dz;

    // 2) 当前单元 P, T
    double P = U[idxC(cell, nc)];     // Pa
    double T = U[idxC(cell, nc + 1)]; // K

    // 3) 闪蒸 => beta, x[i], y[i]
    FlashResult fl = flashCalculation_local(comps, z0, T, P);
    double beta = fl.beta;

    // 4) 每相密度 & 质量基内能 u_g, u_l (J/kg)
    double rho_g = calcRhoVapor(comps, z0, T, P);
    double rho_l = calcRhoLiquid(comps, z0, T, P);
    double u_g = calcInternalEnergy_local(comps, z0, T, P);
    double u_l = calcInternalEnergy_local(comps, z0, T, P);

    // 5) 计算气相、液相质量 (kg)
    double m_g = beta * rho_g * Vc;
    double m_l = (1.0 - beta) * rho_l * Vc;

    // 6) 计算混合速度 v_m：取 cell 左右两个 face 速度平均
    double v_left = Uf[cell];     // cell 左侧面
    double v_right = Uf[cell + 1]; // cell 右侧面
    double v_m = 0.5 * (v_left + v_right);

    // 7) 混合密度 (kg/m^3)
    double rho_mix = beta * rho_g + (1.0 - beta) * rho_l;

    // 8) 内能部分 E_int(J) + 动能部分 E_kin(J)
    double E_internal = m_g * u_g + m_l * u_l;                       // J
    double E_kinetic = 0.5 * rho_mix * v_m * v_m * Vc;             // J

    return E_internal + E_kinetic;
}

double compute_Uto(const SimulationParameters& simParams) {
    double R_pipe = log(simParams.r_o / simParams.r_i) / (2.0 * M_PI * simParams.k_i);   // 1/(2πk_i)·ln(r_o/r_i)
    double R_ann = 1.0 / (2.0 * M_PI * simParams.r_i * simParams.h_i);                   // 内壁对流阻抗 1/(2πr_i h_i)
    double R_cas = log(simParams.r_c / simParams.r_o) / (2.0 * M_PI * simParams.k_c);    // 套管材料导热阻抗 ln(r_c/r_o)/(2πk_c)
    double R_cem = log(simParams.r_m / simParams.r_c) / (2.0 * M_PI * simParams.k_m);    // 水泥环的导热阻抗 ln(r_m/r_c)/(2πk_m)
    double R_form = 1.0 / (2.0 * M_PI * simParams.r_m * simParams.h_o);                  // 外壁对地层对流阻抗 1/(2πr_m h_o)
    double R_total = R_pipe + R_ann + R_cas + R_cem + R_form;
    // U_to = 1/R_total (W/(m^2·K))
    return 1.0 / R_total; 
}

double calcPhaseEnthalpyFlux(
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    double T,
    double P,
    const SimulationParameters& simParams,
    double v_m
) {
    // 1) 闪蒸
    FlashResult fl = flashCalculation_local(comps, z0, T, P);
    double beta = fl.beta;
    auto& x = fl.x, & y = fl.y;

    // 2) 各相密度
    double rho_g = calcRhoVapor(comps, z0, T, P);
    double rho_l = calcRhoLiquid(comps, z0, T, P);

    // 3) 计算各相质量基焓 h_g, h_l
    //    （调用 enthalpyDeparture_mol 得到 J/mol，再除以摩尔质量得到 J/kg，然后加 P/ρ）
    double h_dep_g = enthalpyDeparture_mol(comps, y, T, P, true);   // J/mol
    double h_dep_l = enthalpyDeparture_mol(comps, x, T, P, false);  // J/mol

    double M_g = 0.0, M_l = 0.0;
    int n = comps.size();
    for (int i = 0; i < n; ++i) {
        M_g += y[i] * (comps[i].Mw * 1e-3); // kg/mol
        M_l += x[i] * (comps[i].Mw * 1e-3); // kg/mol
    }
    // 质量基 h = (h_dep/mol)/M + P/ρ
    double h_g = h_dep_g / M_g + P / rho_g; // J/kg
    double h_l = h_dep_l / M_l + P / rho_l; // J/kg

    // 4) 调用漂移模型，计算 v_g, v_l
    DriftResult dr = driftModelDetailed(simParams, v_m, beta, comps, z0, T, P);
    double v_g = dr.gasVelocity;
    double v_l = dr.liquidVelocity;

    // 5) 返回分相焓×速度对流叠加
    //    = ρ_g·h_g·v_g  +  ρ_l·h_l·v_l  (J/(m^2·s))
    return rho_g * h_g * v_g + rho_l * h_l * v_l;
}
// 烃组分质量守恒残差
double massResHc(
    int cell,
    int compIdx,
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    const vector<double>& Uc,
    const vector<double>& Uc_old,
    const SimulationParameters& simParams,
    double dz,
    double dt,
    double vm_im,
    double vm_ip
) {
    int nc = simParams.numComponents;
    int numCvar = nc + 2;
    auto idxC = [&](int i, int eq) {return i * numCvar + eq; };

    const double A = M_PI * simParams.D * simParams.D / 4.0;
    const double M_i = comps[compIdx].Mw * 1e-3;

    double N = Uc[idxC(cell, compIdx)];
    double N0 = Uc_old[idxC(cell, compIdx)];
    double dNdt = (N - N0) / dt;

    double P_c = Uc[idxC(cell, nc)];
    double T_c = Uc[idxC(cell, nc + 1)];
    FlashResult fl_c = flashCalculation_local(comps, z0, T_c, P_c);

    int im = (cell > 0 ? cell - 1 : cell);
    double P_im = Uc[idxC(im, nc)];
    double T_im = Uc[idxC(im, nc + 1)];
    FlashResult fl_im = flashCalculation_local(comps, z0, T_im, P_im);
    DriftResult dr_im = driftModelDetailed(simParams, vm_im, fl_im.beta, comps, z0, T_im, P_im);
    double vg_im = dr_im.gasVelocity, vl_im = dr_im.liquidVelocity;
    double y_im = fl_im.y[compIdx], x_im = fl_im.x[compIdx];
    double rho_g_im = calcRhoVapor(comps, z0, T_im, P_im);
    double rho_l_im = calcRhoLiquid(comps, z0, T_im, P_im);

    int ip = (cell < simParams.wellboreCells - 1 ? cell + 1 : cell);
    double P_ip = (cell < simParams.wellboreCells - 1 ? Uc[idxC(ip, nc)] : P_res);
    double T_ip = (cell < simParams.wellboreCells - 1 ? Uc[idxC(ip, nc + 1)] : T_res);
    FlashResult fl_ip = flashCalculation_local(comps, z0, T_ip, P_ip);
    DriftResult dr_ip = driftModelDetailed(simParams, vm_ip, fl_ip.beta, comps, z0, T_ip, P_ip);
    double vg_ip = dr_ip.gasVelocity, vl_ip = dr_ip.liquidVelocity;
    double y_ip = fl_ip.y[compIdx], x_ip = fl_ip.x[compIdx];
    double rho_g_ip = calcRhoVapor(comps, z0, T_ip, P_ip);
    double rho_l_ip = calcRhoLiquid(comps, z0, T_ip, P_ip);

    double Fm_im = (vg_im * rho_g_im * y_im + vl_im * rho_l_im * x_im) * A / M_i;
    double Fm_ip = (vg_ip * rho_g_ip * y_ip + vl_ip * rho_l_ip * x_ip) * A / M_i;

    double fluxDiv = (Fm_ip - Fm_im) / dz;
    return dNdt + fluxDiv;
}

// ==========================================================
// 水组分质量守恒残差
// ==========================================================
double massResWater(
    int cell,
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    const vector<double>& Uc,
    const vector<double>& Uc_old,
    const SimulationParameters& simParams,
    double dz,
    double dt,
    double vm_im,
    double vm_ip
) {
    int nc = simParams.numComponents;
    int numCvar = nc + 2;
    auto idxC = [&](int i, int eq) {return i * numCvar + eq; };

    const double A = M_PI * simParams.D * simParams.D / 4.0;
    const double M_w = comps[nc - 1].Mw * 1e-3;

    double Nw = Uc[idxC(cell, nc - 1)];
    double Nw0 = Uc_old[idxC(cell, nc - 1)];
    double dNwdt = (Nw - Nw0) / dt;

    int im = (cell > 0 ? cell - 1 : cell);
    double P_im = Uc[idxC(im, nc)];
    double T_im = Uc[idxC(im, nc + 1)];
    FlashResult fl_im = flashCalculation_local(comps, z0, T_im, P_im);
    DriftResult dr_im = driftModelDetailed(simParams, vm_im, fl_im.beta, comps, z0, T_im, P_im);
    double vl_im = dr_im.liquidVelocity;
    double rho_w_im = simParams.rhoWater;

    int ip = (cell < simParams.wellboreCells - 1 ? cell + 1 : cell);
    double P_ip = (cell < simParams.wellboreCells - 1 ? Uc[idxC(ip, nc)] : P_res);
    double T_ip = (cell < simParams.wellboreCells - 1 ? Uc[idxC(ip, nc + 1)] : T_res);
    FlashResult fl_ip = flashCalculation_local(comps, z0, T_ip, P_ip);
    DriftResult dr_ip = driftModelDetailed(simParams, vm_ip, fl_ip.beta, comps, z0, T_ip, P_ip);
    double vl_ip = dr_ip.liquidVelocity;
    double rho_w_ip = simParams.rhoWater;

    double Fw_im = rho_w_im * vl_im * A / M_w;
    double Fw_ip = rho_w_ip * vl_ip * A / M_w;
    double fluxDiv = (Fw_ip - Fw_im) / dz;
    return dNwdt + fluxDiv;
}

// 动量残差（face-centered）
double momentumResFace(
    int f,
    const vector<double>& Uc,
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    const SimulationParameters& simParams,
    double dz,
    const vector<double>& Uf,       // 本时刻面速度
    const vector<double>& Uf_old,   // 上一时刻面速度
    double dt
) {
    int n_c = simParams.numComponents;
    int numCvar = n_c + 2;    // \[0..n_c-1]=N_c, \[n_c]=P, \[n_c+1]=T
    auto idxC = [&](int i, int eq) { return i * numCvar + eq; };

    // 1) 左右压力差
    double P_L = Uc[idxC(f - 1, n_c)];
    double P_R = Uc[idxC(f, n_c)];
    double dP = P_R - P_L;

    // 2) 平均温度和压力
    double T_L = Uc[idxC(f - 1, n_c + 1)];
    double T_R = Uc[idxC(f, n_c + 1)];
    double T_ave = 0.5 * (T_L + T_R);
    double P_ave = 0.5 * (P_L + P_R);

    // 3) 计算混合密度
    double dRhoP, dRhoT;
    double rho_face = computeRhoMix_local(
        comps, z0, T_ave, P_ave,
        1e3, 1e-2, dRhoP, dRhoT,
        simParams
    );

    // 4) 惯性项：ρ * (v_m^{n+1} - v_m^{n})/dt * Δz
    double accel = rho_face * (Uf[f] - Uf_old[f]) / dt * dz;

    // 5) 重力项 + 摩擦项
    double g_term = -rho_face * simParams.g * cos(simParams.theta) * dz;
    double friction = -simParams.frictionFactor
        * rho_face * Uf[f] * fabs(Uf[f])
        / (2.0 * simParams.D) * dz;

    // 总残差
    return accel + dP + g_term + friction;
}

// ==========================================================
// 能量残差（cell-centered）
// ==========================================================
double energyRes(
    int cell,
    const vector<ComponentLocal>& comps,
    const vector<double>& z0,
    const vector<double>& Uc,
    const vector<double>& Uc_old,
    const vector<double>& Uf,
    const SimulationParameters& simParams,
    double dz,
    double dt
) {
    int nc = simParams.numComponents;
    int numCvar = nc + 2;
    auto idxC = [&](int i, int eq) { return i * numCvar + eq; };

    // —— (1) 单元体积 Vc
    double Vc = M_PI * simParams.D * simParams.D / 4.0 * dz;

    // —— (2) 内能 + 动能 累积率 dE/dt per volume (J/(m^3·s))
    double E = computeCellEnergy_local(cell, comps, z0, Uc, Uf, simParams, dz);
    double E0 = computeCellEnergy_local(cell, comps, z0, Uc_old, Uf, simParams, dz);
    double dEdt = (E - E0) / (dt * Vc);

    // —— (3) 对流焓通量散度 per volume (J/(m^3·s))
    //      - 上游 face = Uf[cell]
    //      - 下游 face = Uf[cell+1]
    double P_im = (cell > 0 ? Uc[idxC(cell - 1, nc)] : P_res);
    double T_im = (cell > 0 ? Uc[idxC(cell - 1, nc + 1)] : T_res);
    double P_ip = (cell < simParams.wellboreCells - 1 ? Uc[idxC(cell + 1, nc)] : P_res);
    double T_ip = (cell < simParams.wellboreCells - 1 ? Uc[idxC(cell + 1, nc + 1)] : T_res);
    double vm_im = Uf[cell];       // cell 左面混合速度
    double vm_ip = Uf[cell + 1];   // cell 右面混合速度

    double flux_im = calcPhaseEnthalpyFlux(comps, z0, T_im, P_im, simParams, vm_im); // J/(m^2·s)
    double flux_ip = calcPhaseEnthalpyFlux(comps, z0, T_ip, P_ip, simParams, vm_ip); // J/(m^2·s)
    double fluxDiv = (flux_ip - flux_im) / dz; // J/(m^3·s)

    // —— (4) 串联热阻热损失 per volume (J/(m^3·s))
    //      Q_loss = U_to * As * (T_cell - T_formation)  (J/s)
    double U_to = compute_Uto(simParams); // W/(m^2·K)
    double T_c = Uc[idxC(cell, nc + 1)]; // 单元中心温度 K
    double T_fm = simParams.T0 + simParams.dTdz * (cell * dz);                  // 地层温度，若有深度分布可替换
    // 计算单元侧壁面积 As = 2π * r_i * dz
    double As = 2.0 * M_PI * simParams.r_i * dz;
    double Q_loss = U_to * As * (T_c - T_fm);  // J/s
    double q_loss_vol = Q_loss / Vc;           // J/(m^3·s)

    // —— (5) 最终残差
    return dEdt + fluxDiv - q_loss_vol;
}
void readInputFile(const string& filename, SimulationParameters& simParams) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "无法打开输入文件: " << filename << endl;
        exit(1);
    }
    string line;
    while (getline(fin, line)) {
        // 去掉首尾空格
        auto start = line.find_first_not_of(" \t");
        if (start == string::npos) continue;
        auto end = line.find_last_not_of(" \t");
        string token = line.substr(start, end - start + 1);
        if (token.empty() || token[0] == '#') continue;

        istringstream iss(token);
        string key;
        iss >> key;
        if (key == "simulationType") {
            iss >> simParams.simulationType;
        }
        else if (key == "phaseType") {
            iss >> simParams.phaseType;
        }
        else if (key == "wellboreTotalLength") {
            iss >> simParams.wellboreTotalLength;
        }
        else if (key == "wellboreCells") {
            iss >> simParams.wellboreCells;
        }
        else if (key == "P0") {
            iss >> simParams.P0;
        }
        else if (key == "T0") {
            iss >> simParams.T0;
        }
        else if (key == "g") {
            iss >> simParams.g;
        }
        else if (key == "frictionFactor") {
            iss >> simParams.frictionFactor;
        }
        else if (key == "massFlow") {
            iss >> simParams.massFlow;
        }
        else if (key == "D") {
            iss >> simParams.D;
        }
        else if (key == "dTdz") {
            iss >> simParams.dTdz;
        }
        else if (key == "rhoHCVapor") {
            iss >> simParams.rhoHCVapor;
        }
        else if (key == "rhoHCLiquid") {
            iss >> simParams.rhoHCLiquid;
        }
        else if (key == "rhoWater") {
            iss >> simParams.rhoWater;
        }
        else if (key == "muG") {
            iss >> simParams.muG;
        }
        else if (key == "muL") {
            iss >> simParams.muL;
        }
        else if (key == "xHcTotal") {
            iss >> simParams.xHcTotal;
        }
        else if (key == "xWater_Total") {
            iss >> simParams.xWater_Total;
        }
        else if (key == "numComponents") {
            iss >> simParams.numComponents;
            // 同步分配 vector 长度
            simParams.compNames.resize(simParams.numComponents);
            simParams.compTcrit.resize(simParams.numComponents);
            simParams.compPcrit.resize(simParams.numComponents);
            simParams.compOmega.resize(simParams.numComponents);
            simParams.compMw.resize(simParams.numComponents);
            simParams.xHcFrac.resize(simParams.numComponents);
        }
        else if (key == "compNames") {
            for (int i = 0; i < simParams.numComponents; i++) {
                iss >> simParams.compNames[i];
            }
        }
        else if (key == "compTcrit") {
            for (int i = 0; i < simParams.numComponents; i++) {
                iss >> simParams.compTcrit[i];
            }
        }
        else if (key == "compPcrit") {
            for (int i = 0; i < simParams.numComponents; i++) {
                iss >> simParams.compPcrit[i];
            }
        }
        else if (key == "compOmega") {
            for (int i = 0; i < simParams.numComponents; i++) {
                iss >> simParams.compOmega[i];
            }
        }
        else if (key == "compMw") {
            for (int i = 0; i < simParams.numComponents; i++) {
                iss >> simParams.compMw[i];
            }
        }
        else if (key == "xHcFrac") {
            for (int i = 0; i < simParams.numComponents; i++) {
                iss >> simParams.xHcFrac[i];
            }
        }
        else if (key == "simulationTime") {
            iss >> simParams.simulationTime;
        }
        else if (key == "timeStep") {
            iss >> simParams.timeStep;
        }
        else if (key == "theta_gl") {
            iss >> simParams.theta_gl;
        }
        else if (key == "theta") {
            iss >> simParams.theta;
        }
        else if (key == "heatTransferCoeff") {
            iss >> simParams.heatTransferCoeff;
        }
        else if (key == "r_i") {
            iss >> simParams.r_i;
            }
        else if (key == "r_o") {
                iss >> simParams.r_o;
                }
        else if (key == "k_i") {
                iss >> simParams.k_i;
                }
        else if (key == "h_i") {
                iss >> simParams.h_i;
                }
        else if (key == "r_c") {
                iss >> simParams.r_c;
                }
        else if (key == "k_c") {
                iss >> simParams.k_c;
                }
        else if (key == "r_m") {
                iss >> simParams.r_m;
                }
        else if (key == "k_m") {
                iss >> simParams.k_m;
                }
        else if (key == "h_o") {
                iss >> simParams.h_o;
                }
        // 如果还有其他键，再加 else if 即可……
    }
    fin.close();
}
int main() {
    cout << "******************************************\n"
         << "*                                        *\n"
         << "*            井筒模拟v1.0                *\n"
         << "*                                        *\n"
         << "******************************************\n";
    // 1. 从文件读取所有参数
    SimulationParameters simParams;
    readInputFile("input.txt", simParams);

    // 2. 构造组分局部属性
    vector<ComponentLocal> comps;
    for (int i = 0; i < simParams.numComponents; ++i) {
        ComponentLocal c;
        c.name = simParams.compNames[i];
        c.Tcrit = simParams.compTcrit[i];
        c.Pcrit = simParams.compPcrit[i];
        c.omega = simParams.compOmega[i];
        c.Mw = simParams.compMw[i];
        comps.push_back(c);
    }

    // 3. 计算初始摩尔分率 z0
    int n = simParams.numComponents;
    vector<double> xFrac(n), mols(n), z0(n);
    for (int i = 0; i < n - 1; ++i) xFrac[i] = simParams.xHcFrac[i];
    xFrac[n - 1] = simParams.xWater_Total;
    double sumMol = 0.0;
    for (int i = 0; i < n; ++i) {
        double M = simParams.compMw[i] * 1e-3;
        mols[i] = xFrac[i] / M;
        sumMol += mols[i];
    }
    for (int i = 0; i < n; ++i) z0[i] = mols[i] / sumMol;

    // 4. 网格与时间离散
    int    N = simParams.wellboreCells;
    double L = simParams.wellboreTotalLength;
    double dz = L / N;
    int    nt = int(simParams.simulationTime / simParams.timeStep);
    double dt = simParams.timeStep;

    // 5. 未知量布局
    int nComp = simParams.numComponents;
    int numCvar = nComp + 2;
    int NcUnk = N * numCvar;
    int NfUnk = N + 1;
    int unknowns = NcUnk + NfUnk;

    // 6. 存储分配
    vector<double>         Uc(NcUnk), Uc_old(NcUnk);
    vector<double>         Uf(NfUnk), Uf_old(NfUnk);
    vector<double>         R(unknowns), dU(unknowns);
    vector<vector<double>> J(unknowns, vector<double>(unknowns));
    auto idxC = [&](int cell, int eq) { return cell * numCvar + eq; };
    double A_cross = M_PI * simParams.D * simParams.D / 4.0;

    // 下标常量
    int hydroCount = nComp - 1;
    int waterIdx = hydroCount;
    int pressIdx = nComp;
    int tempIdx = nComp + 1;

    // 7. 初始化
    double dRho_dP, dRho_dT;
    double rho_mix_init = computeRhoMix_local(comps, z0, T_res, P_res,
        1e3, 1e-2,
        dRho_dP, dRho_dT,
        simParams);
    double M_mix_init = 0.0;
    for (int k = 0; k < n; ++k) M_mix_init += z0[k] * (comps[k].Mw * 1e-3);
    double Ni_init = rho_mix_init * A_cross * dz / M_mix_init;
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < nComp; ++k)
            Uc[idxC(i, k)] = z0[k] * Ni_init;
        Uc[idxC(i, pressIdx)] = P_res;
        Uc[idxC(i, tempIdx)] = T_res;
    }
    Uf[0] = 0.0;
    for (int f = 1; f < N; ++f)
        Uf[f] = simParams.massFlow / (simParams.rhoHCLiquid * A_cross);
    Uf[N] = Uf[N - 1];

    ofstream fout("simulation_output.txt");
    fout << fixed << setprecision(4);

    // 8. 时间步循环
    for (int t = 0; t < nt; ++t) {
        Uc_old = Uc;
        Uf_old = Uf;
        auto start = chrono::high_resolution_clock::now();
        // 8.1 20 次 Newton 迭代
        for (int iter = 0; iter < 20; ++iter) {
            for (int i = 0; i < N; ++i) {
                double P_cur = Uc[idxC(i, pressIdx)];
                double T_cur = Uc[idxC(i, tempIdx)];
                if (std::isnan(P_cur) || std::isnan(T_cur)) {
                    cout << ">>> NaN Detected!  时间步 t=" << t
                        << " 迭代 iter=" << iter
                        << " 时，第 " << i << " 个格点温度/压力为 NaN: "
                        << "P = " << P_cur
                        << ", T = " << T_cur << endl;
                    exit(1);
                }
            }
            fill(R.begin(), R.end(), 0.0);
            // (1) Cell 方程
            for (int i = 0; i < N; ++i) {
                double vm_im = Uf[i], vm_ip = Uf[i + 1];
                // 烃组分
                for (int k = 0; k < hydroCount; ++k)
                    R[idxC(i, k)] = massResHc(i, k, comps, z0,
                        Uc, Uc_old,
                        simParams, dz, dt,
                        vm_im, vm_ip);
                // 水组分
                R[idxC(i, waterIdx)] = massResWater(i, comps, z0,
                    Uc, Uc_old,
                    simParams, dz, dt,
                    vm_im, vm_ip);
                // 能量
                R[idxC(i, tempIdx)] = energyRes(
                    i, comps, z0,
                    Uc, Uc_old,
                    Uf,           
                    simParams,
                    dz, dt
                );
                // 状态方程
                double P_i = Uc[idxC(i, pressIdx)];
                double T_i = Uc[idxC(i, tempIdx)];
                double rho_mix = computeRhoMix_local(comps, z0,
                    T_i, P_i,
                    1e3, 1e-2,
                    dRho_dP, dRho_dT,
                    simParams);
                double mass_sum = 0.0;
                for (int k = 0; k < nComp; ++k)
                    mass_sum += Uc[idxC(i, k)] * (comps[k].Mw * 1e-3);
                R[idxC(i, pressIdx)] = rho_mix - mass_sum / (A_cross * dz);
            }
            // 底端 P,T 边界
            R[idxC(N - 1, pressIdx)] = Uc[idxC(N - 1, pressIdx)] - P_res;
            R[idxC(N - 1, tempIdx)] = Uc[idxC(N - 1, tempIdx)] -
                (simParams.T0 + simParams.dTdz * ((N - 1) * dz));
            // 井口（Cell 0）Dirichlet 压力边界
            R[idxC(0, pressIdx)] = Uc[idxC(0, pressIdx)] - simParams.P0;
            // 井口（Cell 0）Dirichlet 温度边界
            R[idxC(0, tempIdx)] = Uc[idxC(0, tempIdx)] - simParams.T0;
            // (2) Face 动量
            for (int f = 1; f < N; ++f)
                R[NcUnk + f] = momentumResFace(f, Uc, comps, z0,
                    simParams, dz,
                    Uf, Uf_old, dt);
            // 井口边界
            double rho_surf = computeRhoMix_local(comps, z0,
                simParams.T0, simParams.P0,
                1e3, 1e-2,
                dRho_dP, dRho_dT,
                simParams);
            R[NcUnk] = Uf[0] - simParams.massFlow / (rho_surf * A_cross);
            // 底端零梯度
            R[NcUnk + N] = Uf[N] - Uf[N - 1];

            // (3) Jacobian & 求解
            const double epsU = 1e-6;
            const double epsP = 1e-3 * max(Uc[idxC(0, pressIdx)], 1.0);
            vector<double> Rp(unknowns);
            for (int j = 0; j < unknowns; ++j) {
                double delta = (j < NcUnk && j % numCvar == pressIdx) ? epsP : epsU;
                if (j < NcUnk) Uc[j] += delta; else Uf[j - NcUnk] += delta;
                // 重新装配 Rp
                fill(Rp.begin(), Rp.end(), 0.0);
                // Cell 方程同上
                for (int i = 0; i < N; ++i) {
                    double vim = Uf[i], vip = Uf[i + 1];
                    for (int k = 0; k < hydroCount; ++k)
                        Rp[idxC(i, k)] = massResHc(i, k, comps, z0, Uc, Uc_old,
                            simParams, dz, dt, vim, vip);
                    Rp[idxC(i, waterIdx)] = massResWater(i, comps, z0, Uc, Uc_old,
                        simParams, dz, dt, vim, vip);
                    Rp[idxC(i, tempIdx)] = energyRes(i, comps, z0, Uc, Uc_old, Uf,
                        simParams, dz, dt);
                    double Pi = Uc[idxC(i, pressIdx)], Ti = Uc[idxC(i, tempIdx)];
                    double rm = computeRhoMix_local(comps, z0, Ti, Pi,
                        1e3, 1e-2,
                        dRho_dP, dRho_dT,
                        simParams);
                    double ms = 0; for (int k = 0; k < nComp; ++k) ms += Uc[idxC(i, k)] * (comps[k].Mw * 1e-3);
                    Rp[idxC(i, pressIdx)] = rm - ms / (A_cross * dz);
                }
                Rp[idxC(N - 1, pressIdx)] = Uc[idxC(N - 1, pressIdx)] - P_res;
                Rp[idxC(N - 1, tempIdx)] = Uc[idxC(N - 1, tempIdx)] -
                    (simParams.T0 + simParams.dTdz * ((N - 1) * dz));
                // Face
                for (int f = 1; f < N; ++f)
                    Rp[NcUnk + f] = momentumResFace(f, Uc, comps, z0,
                        simParams, dz,
                        Uf, Uf_old, dt);
                // 井口
                double rs = computeRhoMix_local(comps, z0, simParams.T0, simParams.P0,
                    1e3, 1e-2,
                    dRho_dP, dRho_dT,
                    simParams);
                Rp[NcUnk] = Uf[0] - simParams.massFlow / (rs * A_cross);
                // 底端
                Rp[NcUnk + N] = Uf[N] - Uf[N - 1];
                // 回撤
                if (j < NcUnk) Uc[j] -= delta; else Uf[j - NcUnk] -= delta;
                for (int i = 0; i < unknowns; ++i)
                    J[i][j] = (Rp[i] - R[i]) / delta;
            }
            solveLinearSystem(J, R, dU);
            for (int k = 0; k < unknowns; ++k) {
                double relax = 0.2 * dU[k];
                if (k < NcUnk) Uc[k] += relax; else Uf[k - NcUnk] += relax;
            }
        }
        // 输出当前时刻结果
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        fout << "Time=" << (t + 1) * dt << "s\n";
        cout << "Time=" << (t + 1) * dt << " s, Duration: " << elapsed.count() << " seconds\n";
        fout << "Cell Depth(m) P(MPa) T(K) vm(m/s) v_gas(m/s) v_liq(m/s) alpha_gas alpha_liq";
        for (auto& name : simParams.compNames) {
            fout << " " << name;
        }
        fout << "\n";
        for (int i = 0; i < N; ++i) {
            double depth = i * dz;
            double P_i = Uc[idxC(i, pressIdx)] / 1e6;  // MPa
            double T_i = Uc[idxC(i, nComp + 1)];      // K
            double vm_i = 0.5 * (Uf[i] + Uf[i + 1]);    // m/s
            FlashResult fl = flashCalculation_local(comps, z0, T_i, P_i * 1e6);
            double rho_g = calcRhoVapor(comps, z0, T_i, P_i * 1e6);
            double rho_l = calcRhoLiquid(comps, z0, T_i, P_i * 1e6);
            double M_g = 0, M_l = 0;
            for (int j = 0; j < n; ++j) {
                M_g += fl.y[j] * (comps[j].Mw * 1e-3);
                M_l += fl.x[j] * (comps[j].Mw * 1e-3);
            }
            double Vm_g = M_g / rho_g;
            double Vm_l = M_l / rho_l;
            double alpha_vol = (fl.beta * Vm_g) / (fl.beta * Vm_g + (1 - fl.beta) * Vm_l);
            DriftResult dr = driftModelDetailed(simParams, vm_i, fl.beta, comps, z0, T_i, P_i * 1e6);
            fout << setw(4) << i + 1
                << setw(10) << depth
                << setw(10) << P_i
                << setw(10) << T_i
                << setw(10) << vm_i
                << setw(12) << dr.gasVelocity
                << setw(12) << dr.liquidVelocity
                << setw(12) << alpha_vol
                << setw(12) << (1 - alpha_vol);
            // 计算并输出每个组分的摩尔分数 Ni / Σ Nj
            double sumN = 0.0;
            for (int k = 0; k < nComp; ++k) {
                sumN += Uc[idxC(i, k)];

            }
            for (int k = 0; k < nComp; ++k) {
                double xin = Uc[idxC(i, k)] / sumN;
                fout << setw(12) << xin;

            }
            fout << "\n";
        }
        fout << "Face velocities:\n";
        for (int f = 0; f <= N; ++f) {
            fout << " f" << f << ":" << Uf[f];
        }
        fout << "\n---------------------------\n";
    }
    fout.close();
    return 0;
}