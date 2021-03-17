//
// Created by vadim on 10.03.2021.
//

#ifndef INTERPOLATION_APROXFRACTRATIONAL_H
#define INTERPOLATION_APROXFRACTRATIONAL_H

#include <vector>
#include <valarray>

enum AproxFractRationalErr_e{
    AFR_ERR_OK = 0U,
    AFR_ERR_SIZE_ERR
};

class AproxFractRational {
    std::vector<double> _x;
    std::vector<double> _f_x;

    int _deg_k;
    int _deg_l;

    std::vector<double> _calc_coefs();
public:
    AproxFractRational() = default;

    AproxFractRationalErr_e Init(std::vector<double> x, std::vector<double> f_x, int deg_k, int deg_l);

    double CalcInPoint(double p, std::vector<double> coef = {});

    std::valarray<double> CalcInPoints(std::valarray<double> p);
};


#endif //INTERPOLATION_APROXFRACTRATIONAL_H
