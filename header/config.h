//
// Created by senbaikang on 13.05.21.
//

#ifndef DATAFILTER_CONFIG_H
#define DATAFILTER_CONFIG_H

#include <set>
#include <string>
#include <vector>
#include <filesystem>

#include "output_data.h"

// Non-greedy pattern that only the very last `.` and everything on the right till the end will be found.
#define SUFFIX_PATTERN "\\..*?$"

class Config
{

private:

  u_int16_t threads = 1;

  // "in"
  // Name of the BAM files used to create the mpileup, i.e., the names of single cells.
  std::string cellNamesFileName;

  // Suffix of cell names in `cellNamesFileName`. [`\\..*?$`]
  std::string cellNameSuffix = SUFFIX_PATTERN;

  // "ex"
  // File name of exclusion list (VCF format), containing loci which should be ignored.
  std::string excludedSitesFileName;

  // "me"
  // File name of mutations to exclude during the sequencing error rate estimation (VCF format).
  std::string excludedMutatedSitesFileName;

  // "inc"
  // File name of inclusion list containing Variants (CHROM, POS) that should be included.
  std::string includedSitesFileName;

  // "mi"
  // File name of inclusion list (VCF format) containing Variants (CHROM, POS, REF, ALT) that should be included.
  std::string includedMutatedSitesFileName;

  // Compressed mpileup files in `*.gz` format or uncompressed normal mpileup files.
  std::vector<std::string> inputFileNames;

  // If the input files are compressed in `*.gz` format, reading length in bytes (2^n) should be set. [524288]
  u_int64_t readLength = 524288;

  // "o"
  // Output file name.
  std::filesystem::path outputFileName;

  // "pr"
  // Prior mutation rate. [0.0001]
  double muRatePrior = 0.0001;

  // Prior germline rate. [0.001]
  double germlineRatePrior = 0.001;

  // "wildMean"
  // Initial Effective sequencing error rate (including amplification error rate and sequencing error rate). [0.001]
  double effectiveSeqErrRate = 0.001;

  // "wildOverDis"
  // Initial overdispersion for wild type. [100.0]
  double wildOverdispersion = 100.0;

  // "muOverDis"
  // Initial overdispersion for mutant type. [2.0]
  double muOverdispersion = 2.0;

  // Prior allelic dropout rate.
  double adoRatePrior = 0.1;

  // "ese"
  // Estimate the effective sequencing error rate (true) or not (false). [true]
  bool estimateErrRate = true;

  // Max number of sites per input file for estimating effective sequencing error rate. [100000]
  u_int32_t errRateEstLoops = 100000;

  // "cwm"
  // Number of tumor cells required to have a mutation in order to be called. [2]
  u_int32_t minTumorCellsToCallMu = 2;

  // "mnp"
  // Number of cells which need to pass the filters described below. [2]
  u_int32_t minTumorCellsToPass = 2;

  // Minimum coverage required per normal cell. [5]
  u_int32_t minCovNormalCell = 5;

  // "mc"
  // Minimum coverage required per cell. [1]
  u_int32_t minCov = 1;

  // "ms"
  // Minimum number of reads required to support the alternative. [3]
  u_int32_t minSup = 3;

  // "mf"
  // Minimum required frequency of reads supporting the alternative per cell. [0]
  double minFreq = 0.0;

  // "mff"
  // Mean of acceptable variant allele frequency across all cells for a specific locus.
  // Mapping artifacts may result in low allele frequencies across cells.
  // In order to filter these events out we apply a log-likelihood ratio test where the sequencing
  // error model has a mean of this value. [0.25]
  double meanFreqSite = 0.25;

  // "bns"
  // Loci with up to this number of alternative supporting reads in the bulk control sample will be skipped as germline. [2]
  u_int32_t maxSupInBulk = 2;

  // "bnc"
  // Minimum required coverage of reads in the bulk control sample. [6]
  u_int32_t minCovInBulk = 6;

  // "ncf"
  // Normal cell filter. Currently there are three options:
  // (0) Do not use the normal cells for filtering;
  // (1) use a simple filtering scheme excluding mutations if the probability of being mutated is higher than not
  // being mutated for any cell independently;
  // (2) filter mutations where the probability that at least one cell is mutated is higher than no cell is mutated.
  // Note that in contrast to (1) the cells are not independent and cells with no alternative support need to be
  // explained via dropout events. [1]
  u_int16_t normalCellFilterMode = 1;

  // "mnc"
  // Maximum number of control cells allowed to be mutated. [0]
  u_int32_t maxNormalCellsToMutate = 0;

  u_int32_t tumorCellNum{};

public:

  Config() = default;

  explicit Config(
      u_int16_t,
      std::string,
      std::vector<std::string>,
      std::string
      );

  explicit Config(
      u_int16_t,
      std::string,
      std::string,
      std::vector<std::string>,
      u_int32_t,
      std::string,
      std::string,
      std::string,
      std::string,
      double,
      double,
      double,
      double,
      double,
      bool,
      u_int32_t,
      u_int32_t,
      u_int32_t,
      u_int32_t,
      u_int32_t,
      double,
      double,
      u_int32_t,
      u_int32_t,
      u_int16_t,
      u_int32_t
      );

  [[nodiscard]] u_int16_t getThreads() const;
  [[nodiscard]] const std::string &getCellNamesFileName() const;
  [[nodiscard]] const std::string &getCellNameSuffix() const;
  [[nodiscard]] const std::string &getExcludedSitesFileName() const;
  [[nodiscard]] const std::string &getExcludedMutatedSitesFileName() const;
  [[nodiscard]] const std::string &getIncludedSitesFileName() const;
  [[nodiscard]] const std::string &getIncludedMutatedSitesFileName() const;
  [[nodiscard]] const std::vector<std::string> &getInputFileNames() const;
  [[nodiscard]] u_int64_t getReadLength() const;
  [[nodiscard]] const std::filesystem::path &getOutputFileName() const;
  [[nodiscard]] double getMuRatePrior() const;
  [[nodiscard]] double getGermlineRatePrior() const;
  [[nodiscard]] double getEffectiveSeqErrRate() const;
  void setEffectiveSeqErrRate(double effectiveSeqErrRate);
  [[nodiscard]] double getWildOverdispersion() const;
  [[nodiscard]] double getMuOverdispersion() const;
  [[nodiscard]] double getAdoRatePrior() const;
  [[nodiscard]] bool isEstimateErrRate() const;
  [[nodiscard]] u_int32_t getErrRateEstLoops() const;
  [[nodiscard]] u_int32_t getMinTumorCellsToCallMu() const;
  [[nodiscard]] u_int32_t getMinTumorCellsToPass() const;
  [[nodiscard]] u_int32_t getMinCovNormalCell() const;
  [[nodiscard]] u_int32_t getMinCov() const;
  [[nodiscard]] u_int32_t getMinSup() const;
  [[nodiscard]] double getMinFreq() const;
  [[nodiscard]] double getMeanFreqSite() const;
  [[nodiscard]] u_int32_t getMaxSupInBulk() const;
  [[nodiscard]] u_int32_t getMinCovInBulk() const;
  [[nodiscard]] u_int16_t getNormalCellFilterMode() const;
  [[nodiscard]] u_int32_t getMaxNormalCellsToMutate() const;
  [[nodiscard]] u_int32_t getTumorCellNum() const;
  void setTumorCellNum(u_int32_t tumorCellNum);

  friend int parseCommandLineArgs(
      Config &,
      int,
      char* []
  );
};

#endif // DATAFILTER_CONFIG_H
