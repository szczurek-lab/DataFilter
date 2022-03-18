//
// Created by senbaikang on 21.05.21.
//

#ifndef DATAFILTER_PROBABILITIES_H
#define DATAFILTER_PROBABILITIES_H

#include <dlib/matrix.h>
#include <numeric>


double logBetaBinCountsTerm(
    double sup,
    double cov
    );

double logBetaBinMixedTerm(
    double sup,
    double cov,
    double mean,
    double overDis
    );

double logBetaBinParamsTerm(
    double mean,
    double overDis
    );

double logBetaBinPDF(
    double,
    double,
    double,
    double
    );

double computeRawWildLogScore(
    Config const &,
    double,
    double
    );

double computeRawHeteroMutLogScore(
    Config const &,
    double,
    double
    );

double computeRawHomoMutLogScore(
    Config const &,
    double,
    double
);

double addLogProb(double, double);

double logNChoose2(u_int32_t);

double logNChooseK(u_int32_t, u_int32_t, double);

double logNChooseK(u_int32_t n, u_int32_t k);

/*
 * This function is used to optimize mean and overdispersion
 * for a given loci
 */
struct OptimizeBetaBinMeanOverDis
{
  const std::vector<std::pair<u_int32_t , u_int32_t>> &counts;

  explicit OptimizeBetaBinMeanOverDis(
      const std::vector<std::pair<u_int32_t, u_int32_t>> &
      );

  double operator()(const dlib::matrix<double,0,1> &) const;

};

struct OptimizeBetaBinMeanOverDisDerivates
{
  const std::vector<std::pair<u_int32_t , u_int32_t>> &counts;

  explicit OptimizeBetaBinMeanOverDisDerivates(
      const std::vector<std::pair<u_int32_t , u_int32_t>> &counts
  );

  dlib::matrix<double> operator()(const dlib::matrix<double,0,1> &) const;

};

struct OptimizeBetaBinOverDis
{
  const std::vector<std::pair<u_int32_t, u_int32_t>> &counts;
  double meanFilter;

  OptimizeBetaBinOverDis(
      const std::vector<std::pair<u_int32_t, u_int32_t>> &,
      double
      );

  double operator()(const dlib::matrix<double,0,1> &) const;

};

struct OptimizeBetaBinOverDisDerivates
{
  const std::vector<std::pair<u_int32_t, u_int32_t>> &counts;
  double meanFilter;

  OptimizeBetaBinOverDisDerivates(
      const std::vector<std::pair<u_int32_t, u_int32_t>> &counts,
      double meanFilter
      );

  dlib::matrix<double> operator()(const dlib::matrix<double,0,1> &) const;

};

template <class T>
decltype(auto) logSumExp(
    T &begin,
    T &end
    )
{
  using R = typename std::iterator_traits<T>::value_type;

  if (begin == end)
    return R{0};

  auto _max = *std::max_element(begin, end);
  auto _sum = std::accumulate(
      begin,
      end,
      R{0},
      [_max](R a, R b)
      {
        return a + std::exp(b - _max);
      }
      );

  return _max + std::log(_sum);
}

template <class T>
T logSumExp(T loga, T logb)
{
  const T _max = std::max(loga, logb);
  const T _sum = std::exp(loga - _max) + std::exp(logb - _max);

  return _max + std::log(_sum);
}

#endif // DATAFILTER_PROBABILITIES_H
