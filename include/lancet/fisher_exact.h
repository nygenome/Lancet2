/*!
 * \brief Fisher exact test implementation mostly taken from htslib.
 *
 *  Copyright (C) 2010, 2013 Genome Research Ltd.
 *  Copyright (C) 2011 Attractive Chaos <attractor@live.co.uk>
 *  https://github.com/samtools/htslib/blob/b39e724b6764d12169c1b22154df4759b4cf2a56/htslib/kfunc.h
 *  https://github.com/samtools/htslib/blob/b39e724b6764d12169c1b22154df4759b4cf2a56/kfunc.c
 */

#pragma once

namespace lancet {
class FisherExact {
 public:
  FisherExact() = default;

  struct Result {
    double leftPVal = 1.0;
    double rightPVal = 1.0;
    double twoSidedPVal = 1.0;
    double probability = 0.0;
  };

  /*!
   *  n11  n12  | n1_
   *  n21  n22  | n2_
   * -----------+----
   *  n_1  n_2  | n
   */
  [[nodiscard]] static auto Test(int n_11, int n_12, int n_21, int n_22) -> Result;

 private:
  struct HgaccT {
    int n11{0};
    int n1_{0};
    int n_1{0};
    int n{0};
    double p{0.0};
  };

  // log binom{n}{k}
  [[nodiscard]] static auto LBinom(int n, int k) -> double;

  /// n11  n12  | n1_
  /// n21  n22  | n2_
  ///-----------+----
  /// n_1  n_2  | n
  /// hypergeometric distribution
  [[nodiscard]] static auto Hypergeo(int n_11, int nrt, int nct, int n) -> double;

  /// incremental version of hypergeometric distribution
  [[nodiscard]] static auto HypergeoAcc(int n_11, int nrt, int nct, int n, HgaccT* aux) -> double;
};

[[nodiscard]] auto PhredScaledProbability(const FisherExact::Result& result) -> double;
}  // namespace lancet
