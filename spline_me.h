/****************************************************************
  @file spline_me.h

  Define helper functions for generating radial matrix elements
  with spline integration.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 08/08/18 (pjf): Created.
****************************************************************/

#ifndef SPLINE_ME_H_
#define SPLINE_ME_H_

#include <cmath>

#include "spline/wavefunction_class.h"

namespace spline
{
enum class BasisType : char {
  kOscillator = 'o',
  kLaguerre = 'l',
};

enum class OperatorType : char {
  kR = 'r',
  kK = 'k',
};

inline double RadialMatrixElement(
    int bra_n, int bra_l, double bra_b, BasisType bra_basis,
    int ket_n, int ket_l, double ket_b, BasisType ket_basis,
    OperatorType operator_type, int operator_order,
    int num_steps = 3000)
// Calculate matrix elements for r and ik operator using spline integration.
//
// Arguments:
//   bra_N, ket_N (input): bra and ket N(=2n+l) quantum numbers
//   bra_l, ket_l (input): bra and ket angular momentum quantum numbers
//   operator_sign (int): sign selecting coordinate or momentum
//     representation (+1 for "r", -1 for "k")
//
// Returns:
//   radial matrix element
{
  // override operator type for overlaps
  if (operator_order == 0) operator_type = OperatorType::kR;

  // select basis type
  Basis bra_basis_type, ket_basis_type;
  if (operator_type == OperatorType::kR) {
    if (bra_basis == BasisType::kOscillator) {
      bra_basis_type = Basis::HC;
    } else if (bra_basis == BasisType::kLaguerre) {
      bra_basis_type = Basis::LC;
    }
    if (ket_basis == BasisType::kOscillator) {
      ket_basis_type = Basis::HC;
    } else if (ket_basis == BasisType::kLaguerre) {
      ket_basis_type = Basis::LC;
    }
  } else if (operator_type == OperatorType::kK) {
    if (bra_basis == BasisType::kOscillator) {
      bra_basis_type = Basis::HM;
    } else if (bra_basis == BasisType::kLaguerre) {
      bra_basis_type = Basis::LM;
    }
    if (ket_basis == BasisType::kOscillator) {
      ket_basis_type = Basis::HM;
    } else if (ket_basis == BasisType::kLaguerre) {
      ket_basis_type = Basis::LM;
    }
  }
  // get bra and ket states
  spline::WaveFunction bra_wavefunction(bra_n, bra_l, bra_b, bra_basis_type);
  spline::WaveFunction ket_wavefunction(ket_n, ket_l, ket_b, ket_basis_type);

  double matrix_element =
      bra_wavefunction.MatrixElement(num_steps, ket_wavefunction, operator_order);

  return matrix_element;
}

};      // namespace spline
#endif  // SPLINE_ME_H_
