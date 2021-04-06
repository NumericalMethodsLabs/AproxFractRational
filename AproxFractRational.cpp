//
// Created by vadim on 10.03.2021.
//

#include "AproxFractRational.h"

#include "Matrix.h"

AproxFractRationalErr_e AproxFractRational::Init(std::vector<double> x, std::vector<double> f_x, int deg_k, int deg_l) {
    AproxFractRationalErr_e err = AFR_ERR_OK;
    if (x.size() != f_x.size() ||
        deg_k + deg_l + 1 != x.size() ||
        deg_k + deg_l + 1 != f_x.size())
        err = AFR_ERR_SIZE_ERR;
    else {
        _x = x;
        _f_x = f_x;
        _deg_k = deg_k;
        _deg_l = deg_l;
    }
    return err;
}

double AproxFractRational::CalcInPoint(double p, std::vector<double> coef) {
    if (coef.empty())
        coef = _calc_coefs();

    int i = 0;
    double tmp_deg = 1,
           numer = 0,
           denumer = 0;

    for ( ; i <= _deg_k; ++i) {
        numer += (coef[i] * tmp_deg);
        tmp_deg *= p;
    }

    tmp_deg = 1;
    for ( ; i <= _deg_k + _deg_l + 1; ++i) {
        denumer += (coef[i] * tmp_deg);
        tmp_deg *= p;
    }

    return numer / denumer;
}

std::valarray<double> AproxFractRational::CalcInPoints(std::valarray<double> p) {
    auto coef = _calc_coefs();
    std::valarray<double> retv(p.size());

    for (unsigned i = 0; i < p.size(); ++i)
        retv[i] = CalcInPoint(p[i]);

    return retv;
}

std::vector<double> AproxFractRational::_calc_coefs() {
    Matrix m(_x.size());
    Matrix vec(_x.size(), 1U);

    for (int i = 0; i < m.count_row(); ++i) {
        double tmp = 1;
        int j = 0;
        for ( ; j < _deg_k; ++j) {
            m[i][j] = _x[i] * tmp;
            tmp *= _x[i];
        }
        tmp = 1;
        for ( ; j <= _x.size(); ++j) {
            m[i][j] = -tmp * _f_x[i];
            tmp *= _x[i];
        }
        vec[i][0] = -1;
    }

    m = make_identity_matrix(m.concatenation(vec));
    auto retv_vec = new double[m.count_col()];
    m.get_col(m.count_col() - 1, retv_vec);

    std::vector<double> retv(retv_vec, retv_vec + _x.size());

    delete [] retv_vec;

    retv.insert(retv.begin(), 1);

    return retv;
}