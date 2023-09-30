//
// Created by senbaikang on 21.05.21.
//

#include <cmath>
#include <config.h>

#include <boost/math/special_functions/digamma.hpp>

#include "probabilities.h"

double logBetaBinCountsTerm(
    double sup,
    double cov
    )
{
  return std::lgamma(cov + 1.0) -  std::lgamma(sup + 1.0) - std::lgamma(cov - sup + 1.0);
}

double logBetaBinMixedTerm(
    double sup,
    double cov,
    double mean,
    double overDis
    )
{
  return std::lgamma(sup + mean * overDis) + std::lgamma(cov - sup + overDis * (1.0 - mean)) - std::lgamma(cov + overDis);
}

double logBetaBinParamsTerm(
    double mean,
    double overDis
    )
{
  return std::lgamma(overDis) - std::lgamma(mean * overDis) - std::lgamma(overDis * (1.0 - mean));
}

double logBetaBinPDF(
    double sup,
    double cov,
    double mean,
    double overDis
)
{
  if (cov == 0)
    return 0;

  return logBetaBinCountsTerm(sup, cov) +
         logBetaBinMixedTerm(sup, cov, mean, overDis) +
         logBetaBinParamsTerm(mean, overDis);
}

double computeRawWildLogScore(
    Config const &config,
    double altCount,
    double coverage
    )
{
  return logBetaBinPDF(
      altCount,
      coverage,
      config.getEffectiveSeqErrRate() / 3.0,
      config.getWildOverdispersion()
      );
//  return logBetaBinPDF(
//      altCount,
//      coverage,
//      config.getEffectiveSeqErrRate(),
//      config.getWildOverdispersion()
//  );
}

double computeRawHeteroMutLogScore(
    Config const &config,
    double altCount,
    double coverage
    )
{
  return logBetaBinPDF(
      altCount,
      coverage,
      0.5 - config.getEffectiveSeqErrRate() / 3.0,
      config.getMuOverdispersion()
      );
//  return logBetaBinPDF(
//      altCount,
//      coverage,
//      0.5 - (2.0 / 3.0 * config.getEffectiveSeqErrRate()),
//      config.getMuOverdispersion()
//  );
}

double computeRawHomoMutLogScore(
    Config const &config,
    double covMinusSup,
    double coverage
    )
{
  return logBetaBinPDF(
      covMinusSup,
      coverage,
      config.getEffectiveSeqErrRate(),
      config.getWildOverdispersion()
      );
}

// Add two values in real space by first exponentiating
double addLogProb(double x, double y)
{
  double maxScore;
  double minScore;

  if (x > y)
  {
    maxScore = x;
    minScore = y;
  }
  else
  {
    maxScore = y;
    minScore = x;
  }

  return std::log(1.0 + std::exp(minScore - maxScore)) + maxScore;
}

double logNChoose2(u_int32_t numMut)
{
  return log(static_cast<double>(numMut)) + log((static_cast<double>(numMut) - 1.0) / 2.0);
}

double logNChooseK(u_int32_t n, u_int32_t k, double logNChoosekMinusOne)
{
  if (k == 0)
    return 0;

  return logNChoosekMinusOne + log(static_cast<double>(n + 1 - k) / static_cast<double>(k));
}

double logNChooseK(u_int32_t n, u_int32_t k)
{
  if (k == 0)
    return 0;

  return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}

OptimizeBetaBinMeanOverDis::OptimizeBetaBinMeanOverDis(
    const std::vector<std::pair<u_int32_t, u_int32_t>> &counts
    ):
        counts(counts)
{}

double OptimizeBetaBinMeanOverDis::operator()(const dlib::matrix<double,0,1> &x) const
{
  double result = 0;
  for (const auto & count : this->counts)
    result += logBetaBinPDF(count.first, count.second, x(0), x(1));

  return result;
}

OptimizeBetaBinMeanOverDisDerivates::OptimizeBetaBinMeanOverDisDerivates(
    const std::vector<std::pair<u_int32_t , u_int32_t>> &counts
    ):
    counts(counts)
{}

dlib::matrix<double> OptimizeBetaBinMeanOverDisDerivates::operator()(const dlib::matrix<double,0,1> &x) const
{
  double mean = x(0);
  double overDis = x(1);
  dlib::matrix<double,0,1> res = {0,0};

  double temp = 0;
  unsigned counter = 0;
  for (const auto & count : this->counts)
  {
    unsigned k = count.first;
    unsigned n = count.second;
    temp += overDis *  boost::math::digamma(k + mean * overDis)
            - overDis *  boost::math::digamma(n - k + overDis - overDis * mean);
    ++counter;
  }
  res(0) = counter * (-overDis *  boost::math::digamma(mean * overDis) + overDis *  boost::math::digamma(overDis - overDis * mean)) + temp;

  temp = 0;
  for (const auto & count : this->counts)
  {
    unsigned k = count.first;
    unsigned n = count.second;
    temp += mean * boost::math::digamma(k + mean * overDis) +
            (1.0 - mean) * boost::math::digamma(n - k + overDis - overDis * mean) -
            boost::math::digamma(n + overDis);
  }
  res(1) = counter * (boost::math::digamma(overDis) - mean * boost::math::digamma(mean * overDis) - (1.0 - mean) * boost::math::digamma(overDis - overDis * mean)) + temp;

  return res;
}

OptimizeBetaBinOverDis::OptimizeBetaBinOverDis(
    const std::vector<std::pair<u_int32_t, u_int32_t>> &counts,
    double meanFilter
):
    counts(counts),
    meanFilter(meanFilter)
{}

double OptimizeBetaBinOverDis::operator()(const dlib::matrix<double,0,1> &x) const
{
  double result = 0;
  for (const auto & count : this->counts)
    result += logBetaBinPDF(count.first, count.second, meanFilter, x(0));

  return result;
}

OptimizeBetaBinOverDisDerivates::OptimizeBetaBinOverDisDerivates(
    const std::vector<std::pair<u_int32_t, u_int32_t>> &counts,
    double meanFilter
    ):
    counts(counts),
    meanFilter(meanFilter)
{}

dlib::matrix<double> OptimizeBetaBinOverDisDerivates::operator()(const dlib::matrix<double,0,1> &x) const
{
  double mean = this->meanFilter;
  double overDis = x(0);
  dlib::matrix<double,0,1> res = {0};

  double temp = 0;
  u_int32_t counter = 0;
  for (const auto & count : this->counts)
  {
    unsigned k = count.first;
    unsigned n = count.second;
    temp += mean * boost::math::digamma(k + mean * overDis) +
            (1.0 - mean) * boost::math::digamma(n - k + overDis - overDis * mean) -
            boost::math::digamma(n + overDis);
    ++counter;
  }
  res(0) = counter * (boost::math::digamma(overDis) - mean * boost::math::digamma(mean * overDis) - (1.0 - mean) * boost::math::digamma(overDis - overDis * mean)) + temp;

  return res;
}
