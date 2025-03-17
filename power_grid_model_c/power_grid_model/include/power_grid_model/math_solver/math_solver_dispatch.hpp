// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

// runtime dispatch for math solver
// this can separate math solver into a different translation unit

#pragma once

#include "../calculation_parameters.hpp"
#include "../common/common.hpp"
#include "../common/exception.hpp"
#include "../common/three_phase_tensor.hpp"
#include "../common/timer.hpp"

#include <memory>

namespace power_grid_model {

namespace math_solver {

// forward declare YBus
template <symmetry_tag sym> class YBus;

// abstract base class
template <symmetry_tag sym> class MathSolverBase {
  public:
    virtual ~MathSolverBase() = default;

    virtual SolverOutput<sym> run_power_flow(PowerFlowInput<sym> const& input, double err_tol, Idx max_iter,
                                             CalculationInfo& calculation_info, CalculationMethod calculation_method,
                                             YBus<sym> const& y_bus) = 0;
    virtual SolverOutput<sym> run_state_estimation(StateEstimationInput<sym> const& input, double err_tol, Idx max_iter,
                                                   CalculationInfo& calculation_info,
                                                   CalculationMethod calculation_method, YBus<sym> const& y_bus) = 0;
    virtual ShortCircuitSolverOutput<sym> run_short_circuit(ShortCircuitInput const& input,
                                                            CalculationInfo& calculation_info,
                                                            CalculationMethod calculation_method,
                                                            YBus<sym> const& y_bus) = 0;
    virtual void clear_solver() = 0;
    virtual void parameters_changed(bool changed) = 0;

  protected:
    MathSolverBase(MathSolverBase const&) = default;
    MathSolverBase& operator=(MathSolverBase const&) = default;
    MathSolverBase(MathSolverBase&&) noexcept = default;
    MathSolverBase& operator=(MathSolverBase&&) noexcept = default;
};

// tag of math solver concrete types
template <template <symmetry_tag> class MathSolverType> struct math_solver_tag {};

class MathSolverDispatcher {
  public:
    template <symmetry_tag sym> struct Config {
        template <template <symmetry_tag> class MathSolverType>
        constexpr Config(math_solver_tag<MathSolverType> /* unused */)
            : create{[](std::shared_ptr<MathModelTopology const> const& topo_ptr) -> void* {
                  return new MathSolverType<sym>{topo_ptr};
              }},
              destroy{[](void const* solver) { delete reinterpret_cast<MathSolverType<sym> const*>(solver); }},
              copy{[](void const* solver) -> void* {
                  return new MathSolverType<sym>{*reinterpret_cast<MathSolverType<sym> const*>(solver)};
              }},
              run_power_flow{[](void* solver, PowerFlowInput<sym> const& input, double err_tol, Idx max_iter,
                                CalculationInfo& calculation_info, CalculationMethod calculation_method,
                                YBus<sym> const& y_bus) {
                  return reinterpret_cast<MathSolverType<sym>*>(solver)->run_power_flow(
                      input, err_tol, max_iter, calculation_info, calculation_method, y_bus);
              }},
              run_state_estimation{[](void* solver, StateEstimationInput<sym> const& input, double err_tol,
                                      Idx max_iter, CalculationInfo& calculation_info,
                                      CalculationMethod calculation_method, YBus<sym> const& y_bus) {
                  return reinterpret_cast<MathSolverType<sym>*>(solver)->run_state_estimation(
                      input, err_tol, max_iter, calculation_info, calculation_method, y_bus);
              }},
              run_short_circuit{[](void* solver, ShortCircuitInput const& input, CalculationInfo& calculation_info,
                                   CalculationMethod calculation_method, YBus<sym> const& y_bus) {
                  return reinterpret_cast<MathSolverType<sym>*>(solver)->run_short_circuit(input, calculation_info,
                                                                                           calculation_method, y_bus);
              }},
              clear_solver{[](void* solver) { reinterpret_cast<MathSolverType<sym>*>(solver)->clear_solver(); }},
              parameters_changed{[](void* solver, bool changed) {
                  reinterpret_cast<MathSolverType<sym>*>(solver)->parameters_changed(changed);
              }} {}

        std::add_pointer_t<void*(std::shared_ptr<MathModelTopology const> const&)> create;
        std::add_pointer_t<void(void const*)> destroy;
        std::add_pointer_t<void*(void const*)> copy;
        std::add_pointer_t<SolverOutput<sym>(void*, PowerFlowInput<sym> const&, double, Idx, CalculationInfo&,
                                             CalculationMethod, YBus<sym> const&)>
            run_power_flow;
        std::add_pointer_t<SolverOutput<sym>(void*, StateEstimationInput<sym> const&, double, Idx, CalculationInfo&,
                                             CalculationMethod, YBus<sym> const&)>
            run_state_estimation;
        std::add_pointer_t<ShortCircuitSolverOutput<sym>(void*, ShortCircuitInput const&, CalculationInfo&,
                                                         CalculationMethod, YBus<sym> const&)>
            run_short_circuit;
        std::add_pointer_t<void(void*)> clear_solver;
        std::add_pointer_t<void(void*, bool)> parameters_changed;
    };

    template <template <symmetry_tag> class MathSolverType>
    constexpr MathSolverDispatcher(math_solver_tag<MathSolverType> /* unused */)
        : sym_config_{math_solver_tag<MathSolverType>{}}, asym_config_{math_solver_tag<MathSolverType>{}} {}

    template <symmetry_tag sym> Config<sym> const& get_dispatcher_config() const {
        if constexpr (is_symmetric_v<sym>) {
            return sym_config_;
        } else {
            return asym_config_;
        }
    }

  private:
    Config<symmetric_t> sym_config_;
    Config<asymmetric_t> asym_config_;
};

template <symmetry_tag sym> class MathSolverProxy {
  public:
    explicit MathSolverProxy(MathSolverDispatcher const* dispatcher,
                             std::shared_ptr<MathModelTopology const> const& topo_ptr)
        : dispatcher_{dispatcher},
          solver_{dispatcher_->get_dispatcher_config<sym>().create(topo_ptr),
                  dispatcher_->get_dispatcher_config<sym>().destroy} {}
    MathSolverProxy(MathSolverProxy const& other)
        : dispatcher_{other.dispatcher_},
          solver_{dispatcher_->get_dispatcher_config<sym>().copy(other.get_ptr()),
                  dispatcher_->get_dispatcher_config<sym>().destroy} {}
    MathSolverProxy& operator=(MathSolverProxy const& other) {
        if (this != &other) {
            solver_.reset();
            dispatcher_ = other.dispatcher_;
            solver_ = std::unique_ptr<void, std::add_pointer_t<void(void const*)>>{
                dispatcher_->get_dispatcher_config<sym>().copy(other.get_ptr()),
                dispatcher_->get_dispatcher_config<sym>().destroy};
        }
        return *this;
    }
    MathSolverProxy(MathSolverProxy&& other) noexcept = default;
    MathSolverProxy& operator=(MathSolverProxy&& other) noexcept = default;
    ~MathSolverProxy() = default;

    SolverOutput<sym> run_power_flow(PowerFlowInput<sym> const& input, double err_tol, Idx max_iter,
                                     CalculationInfo& calculation_info, CalculationMethod calculation_method,
                                     YBus<sym> const& y_bus) {
        return dispatcher_->get_dispatcher_config<sym>().run_power_flow(get_ptr(), input, err_tol, max_iter,
                                                                        calculation_info, calculation_method, y_bus);
    }

    SolverOutput<sym> run_state_estimation(StateEstimationInput<sym> const& input, double err_tol, Idx max_iter,
                                           CalculationInfo& calculation_info, CalculationMethod calculation_method,
                                           YBus<sym> const& y_bus) {
        return dispatcher_->get_dispatcher_config<sym>().run_state_estimation(
            get_ptr(), input, err_tol, max_iter, calculation_info, calculation_method, y_bus);
    }

    ShortCircuitSolverOutput<sym> run_short_circuit(ShortCircuitInput const& input, CalculationInfo& calculation_info,
                                                    CalculationMethod calculation_method, YBus<sym> const& y_bus) {
        return dispatcher_->get_dispatcher_config<sym>().run_short_circuit(get_ptr(), input, calculation_info,
                                                                           calculation_method, y_bus);
    }

    void clear_solver() { dispatcher_->get_dispatcher_config<sym>().clear_solver(get_ptr()); }

    void parameters_changed(bool changed) {
        dispatcher_->get_dispatcher_config<sym>().parameters_changed(get_ptr(), changed);
    }

  private:
    MathSolverDispatcher const* dispatcher_{};
    std::unique_ptr<void, std::add_pointer_t<void(void const*)>> solver_;

    void* get_ptr() { return solver_.get(); }
    void const* get_ptr() const { return solver_.get(); }
};

} // namespace math_solver

template <symmetry_tag sym> using MathSolverProxy = math_solver::MathSolverProxy<sym>;

using MathSolverDispatcher = math_solver::MathSolverDispatcher;

} // namespace power_grid_model
