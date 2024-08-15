#ifndef EXECUTOR_H
#define EXECUTOR_H

#include <chrono>
#include <iterator>
#include <type_traits>
#include <utility> // std::move

#if USE_BENCHMARK
#include <benchmark/benchmark.h>
#endif /* USE_BENCHMARK */

#include "grackle.h"

#include "FieldData.h"
#include "operation.h"
#include "utils.h"

// this is basically a dummy object to let us operate without installing
// google benchmark
namespace grackle {

class BenchState {

  int max_iter_;
  double total_elapsed_seconds_;

public:
  BenchState(int max_iter)
    : max_iter_(max_iter), total_elapsed_seconds_(0)
  {
    GRCLI_REQUIRE(max_iter > 0, "max_iter must be positive");
  }

  void PauseTiming() { }
  void ResumeTiming() { }

  void SetIterationTime(double num_seconds) {
    total_elapsed_seconds_ += num_seconds;
  }

  double GetTotalElapsedSeconds() {
    return total_elapsed_seconds_;
  }

  class iterator {
    int count_;

  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = int;
    using difference_type = int;
    using pointer = int*;
    using reference = int&;

    explicit iterator(int count) : count_(count) { }
    iterator& operator++() { *this = iterator(count_+1); return *this;}
    iterator operator++(int) { iterator out = *this; ++(*this); return out; }
    bool operator==(iterator other) const { return count_ == other.count_; }
    bool operator!=(iterator other) const { return count_ != other.count_; }
    reference operator*() { return count_; } // should never be nullptr
  };

  iterator begin() {return iterator(0);}
  iterator end() { return iterator(max_iter_); }
};

} // namespace grackle

class GrackleDriver {

  chemistry_data* my_chem_;
  chemistry_data_storage* my_rates_;
  code_units my_units_;
  FieldData wrapped_my_field_;
  OperationSpec operation_;

  
  /// This template function does the actual work
  template<typename T, OperationKind op>
  void helper_(T& state) {

    double dt = this->operation_.dt;
    // do any setup!
    chemistry_data* my_chem = this->my_chem_;
    chemistry_data_storage* my_rates = this->my_rates_;
    code_units* my_units = &my_units_;

    // set the following up!
    grackle_field_data* copied_fields = impl::allocate_and_init_gr_field_data(
        *my_chem,
        wrapped_my_field_.get_ptr()->grid_rank,
        wrapped_my_field_.get_ptr()->grid_dimension);


    gr_float* out = new gr_float[wrapped_my_field_.grid_size()];

    using time_point = 
      std::chrono::time_point<std::chrono::high_resolution_clock>;

    for (auto _ : state) {
      // declare 2 variables that are sometimes used
      [[maybe_unused]] time_point t_start;
      [[maybe_unused]] time_point t_end;

      // because we modify field values in place in all operations (even when
      // simply doing unit conversions), the input values will vary between
      // iterations (the amount of drift obviously depends upon the operation)
      // -> for consistency, we create a new copy of the fields at the start of
      //    each loop
      // -> right now, we are using PauseTiming and ResumeTiming, but it may be
      //    better to manually time these things
      state.PauseTiming();
      clone_field_data(copied_fields, wrapped_my_field_.get_ptr());
      state.ResumeTiming();

      if constexpr (std::is_same_v<T, grackle::BenchState>) {
        t_start = std::chrono::high_resolution_clock::now();
      }

      // here is where we do the actual work
      if constexpr (op == OperationKind::calc_cooling_time) {
        local_calculate_cooling_time(my_chem, my_rates, my_units,
                                     copied_fields, out);
      } else if constexpr (op == OperationKind::calc_dust_temperature) {
        local_calculate_dust_temperature(my_chem, my_rates, my_units,
                                         copied_fields, out);
      } else if constexpr (op == OperationKind::calc_pressure) {
        local_calculate_pressure(my_chem, my_rates, my_units,
                                 copied_fields, out);
      } else if constexpr (op == OperationKind::calc_temperature) {
        local_calculate_temperature(my_chem, my_rates, my_units,
                                    copied_fields, out);
      } else if constexpr (op == OperationKind::solve_chemistry) {
        local_solve_chemistry(my_chem, my_rates, my_units, copied_fields, dt);
      }

      if constexpr (std::is_same_v<T, grackle::BenchState>) {
        t_end = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds elapsed_nanoseconds =
          std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start);
        state.SetIterationTime(elapsed_nanoseconds.count() / 1e9);
      }
    }

    // this is where we cleanup
    impl::destroy_selfcontained_field_data(copied_fields);
    delete[] out;
  }

public:

  GrackleDriver(chemistry_data* my_chem, chemistry_data_storage* my_rates,
                code_units my_units,
                FieldData&& wrapped_fields, OperationSpec operation)
    : my_chem_(my_chem),
      my_rates_(my_rates),
      my_units_(my_units),
      wrapped_my_field_(std::move(wrapped_fields)),
      operation_(operation)
  {
  }

  /// function that actually dispatches to the correct specialization of the
  /// `helper_` function
  template <typename T>
  void operator()(T& state) {
    switch(this->operation_.kind) {
      case OperationKind::calc_cooling_time:
        return helper_<T, OperationKind::calc_cooling_time>(state);
      case OperationKind::calc_dust_temperature:
        return helper_<T, OperationKind::calc_dust_temperature>(state);
      case OperationKind::calc_pressure:
        return helper_<T, OperationKind::calc_pressure>(state);
      case OperationKind::calc_temperature:
        return helper_<T, OperationKind::calc_temperature>(state);
      case OperationKind::solve_chemistry:
        return helper_<T, OperationKind::solve_chemistry>(state);
    }
  }

};

#endif /* EXECUTOR_H */