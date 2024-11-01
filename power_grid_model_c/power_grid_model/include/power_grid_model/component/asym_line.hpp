#pragma once

#include <numeric>

#include "branch.hpp"

#include "../auxiliary/input.hpp"
#include "../auxiliary/output.hpp"
#include "../auxiliary/update.hpp"
#include "../calculation_parameters.hpp"
#include "../common/common.hpp"
#include "../common/three_phase_tensor.hpp"

namespace power_grid_model {

class AsymLine : public Branch {
  public:
    using InputType = AsymLineInput;
    using UpdateType = BranchUpdate;
    static constexpr char const* name = "asym_line";

    explicit AsymLine(AsymLineInput const& asym_line_input, double system_frequency, double u1, double u2)
        : Branch{asym_line_input}, i_n_{asym_line_input.i_n}, base_i_{base_power_3p / u1 / sqrt3} {
        if (cabs(u1 - u2) > numerical_tolerance) {
            throw ConflictVoltage{id(), from_node(), to_node(), u1, u2};
        }
        ComplexTensor<asymmetric_t> z_series;
        ComplexTensor<asymmetric_t> c_matrix;
        if (is_nan(asym_line_input.r_na)) 
        {
           ComplexTensor<asymmetric_t> r_matrix = ComplexTensor<asymmetric_t>(asym_line_input.r_aa, asym_line_input.r_ba, asym_line_input.r_bb, asym_line_input.r_ca, asym_line_input.r_cb, asym_line_input.r_cc);
           ComplexTensor<asymmetric_t> x_matrix = ComplexTensor<asymmetric_t>(asym_line_input.x_aa, asym_line_input.x_ba, asym_line_input.x_bb, asym_line_input.x_ca, asym_line_input.x_cb, asym_line_input.x_cc);
           z_series = r_matrix + x_matrix * 1.0i;

           c_matrix = ComplexTensor<asymmetric_t>(asym_line_input.c_aa, asym_line_input.c_ba, asym_line_input.c_bb, asym_line_input.c_ca, asym_line_input.c_cb, asym_line_input.c_cc);
        }
        else 
        {
            ComplexTensor4 r_matrix = ComplexTensor4(asym_line_input.r_aa, asym_line_input.r_ba, asym_line_input.r_bb, asym_line_input.r_ca, asym_line_input.r_cb, asym_line_input.r_cc, asym_line_input.r_na, asym_line_input.r_nb, asym_line_input.r_nc, asym_line_input.r_nn);
            ComplexTensor4 x_matrix = ComplexTensor4(asym_line_input.x_aa, asym_line_input.x_ba, asym_line_input.x_bb, asym_line_input.x_ca, asym_line_input.x_cb, asym_line_input.x_cc, asym_line_input.x_na, asym_line_input.x_nb, asym_line_input.x_nc, asym_line_input.x_nn);
            
            ComplexTensor4 y = r_matrix + 1.0i * x_matrix;
            
            ComplexTensor<asymmetric_t> Y_reduced = this->kron_reduction(y);

            // Continue with a calculation
            DoubleComplex a = std::pow(e, 1.0*(2/3)* pi);
            ComplexTensor<asymmetric_t> a_matrix = ComplexTensor<asymmetric_t>(1, 1, pow(a, 2), 1, a, pow(a,2));
            ComplexTensor<asymmetric_t> a_matrix_inv = a_matrix.inverse();
            z_series = a_matrix_inv * Y_reduced * a_matrix;
            ComplexTensor4 c_matrix_neutral = ComplexTensor4(asym_line_input.c_aa, asym_line_input.c_ba, asym_line_input.c_bb, asym_line_input.c_ca, asym_line_input.c_cb, asym_line_input.c_cc, asym_line_input.c_na, asym_line_input.c_nb, asym_line_input.c_nc, asym_line_input.c_nn);
            c_matrix = this->kron_reduction(c_matrix_neutral);
        }

        y_series = inv(z_series);
        y_shunt = 2 * pi * system_frequency * c_matrix * 1.0i;
    }

    // override getter
    double base_i_from() const override { return base_i_; }
    double base_i_to() const override { return base_i_; }
    double loading(double /* max_s */, double max_i) const override { return max_i / i_n_; };
    double phase_shift() const override { return 0.0; }
    bool is_param_mutable() const override { return false; }

  private:
    double i_n_;
    double base_i_;
    ComplexTensor<asymmetric_t> y_series;
    ComplexTensor<asymmetric_t> y_shunt;

    // You might want to move these functions to some kind of helper function class
    ComplexTensor<asymmetric_t> kron_reduction(const ComplexTensor4& matrix_to_reduce) const {
        ComplexTensor4 Y = matrix_to_reduce;
        ComplexTensor<asymmetric_t> Y_aa = ComplexTensor<asymmetric_t>(Y(0,0), Y(1,0), Y(1, 1), Y(2,0), Y(2,1), Y(2,2)); // 3 * 3
        ComplexValue<asymmetric_t> Y_ab = ComplexValue<asymmetric_t>(Y(0,3), Y(1,3), Y(2,3)); // 3 * 1
        Eigen::Matrix<DoubleComplex, 1, 3> Y_ba; // 1 * 3
        Eigen::Vector3d columnVector;
        columnVector << 1, 2, 3;
        Y_ba << Y(3,0), Y(3,1), Y(3,2);
        DoubleComplex Y_bb = Y(3,3);
        DoubleComplex Y_bb_inv = 1.0 / Y_bb; // 1 * 1
        return Y_aa - ((Y_ab * Y_bb_inv) * Y_ba);
    }

    DoubleComplex average_of_diagonal_of_matrix(const ComplexTensor<asymmetric_t> &matrix) const {
        return (matrix(0,0) + matrix(1,1) + matrix(2,2)) / 3.0;
    }

    DoubleComplex average_of_off_diagonal_of_matrix(const ComplexTensor<asymmetric_t> &matrix) const {
        return (matrix(0,2) + matrix(1,1) + matrix(2,0)) / 3.0;
    }

    BranchCalcParam<symmetric_t> sym_calc_param() const override {
        DoubleComplex y1_series_ = average_of_diagonal_of_matrix(y_series) - average_of_off_diagonal_of_matrix(y_series);
        DoubleComplex y1_shunt_ = average_of_diagonal_of_matrix(y_shunt) - average_of_off_diagonal_of_matrix(y_shunt);

        // return calc_param_y_sym(y1_series_, y1_shunt_, 1.0);
        return calc_param_y_sym(0.0, 0.0, 1.0);
    }

    BranchCalcParam<asymmetric_t> asym_calc_param() const override {
        BranchCalcParam<asymmetric_t> param{};
        // not both connected
        if (!branch_status()) {
            // single connected
            if (from_status() || to_status()) {
                // branch_shunt = 0.5 * y_shunt + 1.0 / (1.0 / y_series + 2.0 / y_shunt);
                ComplexTensor<asymmetric_t> branch_shunt = 0.5 * inv(y_shunt) + inv(inv(y_series) + 2.0 * inv(y_shunt));
                // from or to connected
                param.yff() = from_status() ? branch_shunt : ComplexTensor<asymmetric_t>();
                param.ytt() = to_status() ? branch_shunt :  ComplexTensor<asymmetric_t>();
            }
        }
        // both connected
        else {
            param.ytt() = y_series + 0.5 * y_shunt;
            param.yff() = param.ytt();
            param.yft() = y_series;
            param.ytf() = y_series;
        }
        return param;
    }
};
}