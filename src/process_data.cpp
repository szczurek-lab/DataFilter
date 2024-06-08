//
// Created by senbaikang on 14.05.21.
//

#include <cfloat>
#include <iostream>
#include <regex>

#include <dlib/global_optimization.h>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <datafilter/config.h>
#include <datafilter/thread_pool.h>
#include <datafilter/output_data.h>
#include <datafilter/process_data.h>
#include <datafilter/probabilities.h>

void process_data(Config &config) {
  // Build a thread pool.
  ThreadPool threadPool(config.getThreads());

  // Get pointers to data files.
  auto filePtrs(std::move(getReadCountsFiles(config)));

  // Object storing output data.
  Data data;

  // Mutex for data synchronization.
  std::mutex dataMutex;

  // Read the cellular information.
  std::vector<u_int32_t> tumorCellPos;
  std::vector<u_int32_t> normalCellPos;
  u_int32_t tumorBulkPos;
  u_int32_t normalBulkPos;
  readCellInformation(
          tumorCellPos,
          normalCellPos,
          tumorBulkPos,
          normalBulkPos,
          config,
          data
  );

  // Get the number of tumor cells.
  const u_int32_t tumorCellNum = tumorCellPos.size();
  config.setTumorCellNum(tumorCellNum);

  // Get sites to skip.
  std::set<std::pair<ChromosomeLabel, u_int32_t>> exMap;
  if (!config.getExcludedSitesFileName().empty()) {
    std::cout << "> Reading sites to be excluded... ";
    readSpecialSites(config.getExcludedSitesFileName(), exMap);
    std::cout << "Done. " << exMap.size() << " sites will be excluded." << std::endl;
  }

  // Get mutated sites to skip during the error rate estimation.
  std::set<std::pair<ChromosomeLabel, u_int32_t>> exMuMap;
  if (!config.getExcludedMutatedSitesFileName().empty()) {
    std::cout << "> Reading sites to be excluded... ";
    readSpecialSites(config.getExcludedMutatedSitesFileName(), exMuMap);
    std::cout << "Done. " << exMuMap.size() << " sites will be excluded." << std::endl;
  }

  // Get sites to keep.
  std::set<std::pair<ChromosomeLabel, u_int32_t>> incMap;
  if (!config.getIncludedSitesFileName().empty()) {
    std::cout << "> Reading sites to be included... ";
    readSpecialSites(config.getIncludedSitesFileName(), incMap);
    std::cout << "Done. " << incMap.size() << " sites will be included." << std::endl;
  }

  // Get mutated sites to keep.
  std::set<std::tuple<ChromosomeLabel, u_int32_t, char, char>> incMuMap;
  if (!config.getIncludedMutatedSitesFileName().empty())
    readInclusionVCF(config.getIncludedMutatedSitesFileName(), incMuMap);

  // Estimate effective error rate if necessary using multithreading.
  if (config.isEstimateErrRate()) {
    // <error, cov>
    u_int32_t _error = 0;
    u_int32_t _cov = 0;
    std::vector<std::future<std::pair<u_int32_t, u_int32_t>>> collection{};

    for (auto &i: filePtrs)
      collection.emplace_back(
              threadPool.addTask(
                      estimateSeqErrorRate,
                      std::ref(config),
                      std::ref(i),
                      std::ref(exMap),
                      std::ref(exMuMap),
                      std::ref(tumorCellPos),
                      std::ref(normalCellPos)
              )
      );

    for (auto &&i: collection) {
      const auto &j = i.get();

      _error += std::get<0>(j);
      _cov += std::get<1>(j);
    }

    config.setEffectiveSeqErrRate(
            static_cast<double>(_error) / static_cast<double>(_cov)
    );
  }

  // Multithreading.
  for (auto &i: filePtrs)
    threadPool.addTask(
            process_single_file,
            std::ref(config),
            std::ref(i),
            std::ref(tumorCellNum),
            std::ref(tumorCellPos),
            std::ref(normalCellPos),
            std::ref(tumorBulkPos),
            std::ref(normalBulkPos),
            std::ref(exMap),
            std::ref(incMap),
            std::ref(incMuMap),
            std::ref(dataMutex),
            std::ref(data)
    );

  std::cout << threadPool.waitUntilFinished() << " files processed." << std::endl;

  // Write to output file.
  data.createNoiseCounts();
  data.sortEntries();

  try {
    writeAltNucInfo(data, config);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
  }

  closeReadCountsFile(filePtrs);
}

/**
 * Get indices of different cells / cell groups.
 *
 * @param config        configuration object
 * @param data          object storing output data
 */
void readCellInformation(
        std::vector<u_int32_t> &tumorCellPos,
        std::vector<u_int32_t> &normalCellPos,
        u_int32_t &tumorBulkPos,
        u_int32_t &normalBulkPos,
        const Config &config,
        Data &data
) {
  std::cout << "> Reading cell names... ";

  readCellNames(config, data);

  tumorBulkPos = UINT_MAX;
  normalBulkPos = UINT_MAX;

  std::ifstream inputStream(config.getCellNamesFileName());

  std::vector<std::string> splitVec;
  u_int32_t counter = 0;
  std::string currLine;

  while (getline(inputStream, currLine)) {
    if (!currLine.empty()) {
      boost::split(splitVec, currLine, boost::is_any_of("\t"));

      if (splitVec[1] == "BN") {
        if (normalBulkPos != UINT_MAX)
          std::cout << "WARNING: Multiple bulk control files provided. Only using the first one!";
        else
          normalBulkPos = counter;
      }

      // TODO:
      if (splitVec[1] == "BT") {
      }

      if (splitVec[1] == "CN")
        normalCellPos.emplace_back(counter);

      if (splitVec[1] == "CT")
        tumorCellPos.emplace_back(counter);

      ++counter;
    }
  }

  std::cout << "Done." << std::endl;
}

/**
 * Get names of cancer cells.
 *
 * @param config configuration object
 * @param data   object storing output data
 */
void readCellNames(const Config &config, Data &data) {
  std::ifstream inputStream(config.getCellNamesFileName());
  std::vector<std::string> splitVec;
  std::vector<std::string> splitVecEntry;

  std::string currLine;

  while (getline(inputStream, currLine)) {
    if (!currLine.empty()) {
      boost::split(splitVec, currLine, boost::is_any_of("\t"));

      if (splitVec[1] == "CT") {
        boost::split(splitVecEntry, splitVec[0], boost::is_any_of("/"));

        std::string name = splitVecEntry.back();
        std::regex suf(config.getCellNameSuffix());

        if (std::regex_search(name, suf))
          name = std::regex_replace(name, suf, "");

        data.addCellName(std::move(name));
      }
    }
  }
}

/**
 * Get sites which will be treated specially (excluded or included).
 */
void readSpecialSites(
        const std::string &v,
        std::set<std::pair<ChromosomeLabel, u_int32_t>> &_map
) {
  std::ifstream inputStream(v);
  if (!inputStream)
    throw std::runtime_error("Error: Could not open " + v);

  std::string currLine;
  std::vector<std::string> splitVec;

  while (std::getline(inputStream, currLine)) {
    if (currLine[0] != '#') {
      boost::split(splitVec, currLine, boost::is_any_of("\t"));

      _map.emplace(
              std::make_pair(
                      std::move(ChromosomeLabel(std::move(splitVec[0]))),
                      std::stoi(splitVec[1])
              )
      );
    }
  }
}

/**
 * Get sites which will be included.
 */
void readInclusionVCF(
        const std::string &v,
        std::set<std::tuple<ChromosomeLabel, u_int32_t, char, char>> &_map
) {
  std::cout << "> Reading sites to be included... ";

  std::ifstream inputStream(v);
  if (!inputStream)
    throw std::runtime_error("Error: Could not open " + v);

  std::string currLine;
  std::vector<std::string> splitVec;

  while (std::getline(inputStream, currLine)) {
    if (currLine[0] != '#') {
      boost::split(splitVec, currLine, boost::is_any_of("\t"));

      if (splitVec[3].size() != 1)
        throw std::runtime_error("Error: Illegal format of reference nucleotide: " + splitVec[3]);

      for (auto i: splitVec[4]) {
        _map.emplace(
                std::make_tuple(
                        std::move(ChromosomeLabel(std::move(splitVec[0]))),
                        std::stoi(splitVec[1]),
                        splitVec[3][0],
                        i
                )
        );
      }
    }
  }

  std::cout << "Done. " << _map.size() << " sites will be included." << std::endl;
}

std::vector<std::unique_ptr<ReadCountsFile>> getReadCountsFiles(const Config &config) {
  const std::regex name(GZ_FILE_PATTERN);

  std::vector<std::unique_ptr<ReadCountsFile>> ret;

  for (const std::string &i: config.getInputFileNames()) {
    if (std::regex_match(i, name))
      ret.emplace_back(std::unique_ptr<ReadCountsFile>(new GZFile(i, config.getReadLength())));
    else
      ret.emplace_back(std::unique_ptr<ReadCountsFile>(new MPileupFile(i)));
  }

  return ret;
}

void closeReadCountsFile(std::vector<std::unique_ptr<ReadCountsFile>> &files) {
  for (auto &i: files)
    i.reset();
}

// This function is used to collect the sequencing information of all cells.
void extractSeqInformation(
        std::vector<std::array<u_int32_t, 5>> &counts,
        const std::vector<std::string> &splitVec,
        const std::vector<u_int32_t> &positions
) {
  // loop over the cells
  for (size_t i = 0; i < positions.size(); ++i)
    extractSeqInformation(counts[i], splitVec, positions[i]);
}

// This functions is used to extract information about
// sequencing errors. Only if a single cell shows a mutation
// the error rates will be collected.
void updateSeqErrorStats(
        u_int32_t &seqErrors,
        u_int32_t &seqErrorsCombCov,
        const std::vector<std::array<u_int32_t, 5>> &counts
) {
  for (size_t j = 0; j < 4; ++j) // loop over the nucleotides
    for (const auto &count: counts) // loop over all cells
      seqErrors += count[j];

  for (const auto &count: counts) // update the coverage
    seqErrorsCombCov += count[4];
}

std::pair<u_int32_t, u_int32_t> estimateSeqErrorRate(
        const Config &config,
        std::unique_ptr<ReadCountsFile> &filePtr,
        const std::set<std::pair<ChromosomeLabel, u_int32_t>> &exMap,
        const std::set<std::pair<ChromosomeLabel, u_int32_t>> &errExMap,
        const std::vector<u_int32_t> &tumorCellPos,
        const std::vector<u_int32_t> &normalCellPos
) {
  u_int32_t maxEstLoop = config.getErrRateEstLoops();

  std::string currLine;
  std::vector<std::string> splitVec;
  std::vector<std::array<u_int32_t, 5>> tumor_counts(tumorCellPos.size(), {{0, 0, 0, 0, 0}});
  std::vector<std::array<u_int32_t, 5>> normal_counts(normalCellPos.size(), {{0, 0, 0, 0, 0}});

  u_int32_t seqErrors = 0;
  u_int32_t seqErrorsCombCov = 0;

  for (size_t lineNumber = 0; lineNumber < maxEstLoop && filePtr->getLine(currLine); lineNumber++) {
    if (currLine.empty()) {
      lineNumber--;
      continue;
    }

    boost::split(splitVec, currLine, boost::is_any_of("\t"));

    if (!isRefKnown(splitVec[2])) {
      lineNumber--;
      continue;
    }

    const auto &tmp = std::make_pair(
            std::move(ChromosomeLabel(splitVec[0])),
            std::stoi(splitVec[1])
    );
    auto itEx = exMap.find(tmp);
    auto itErrEx = errExMap.find(tmp);

    if (itEx == exMap.end() && itErrEx == errExMap.end()) {
      extractSeqInformation(tumor_counts, splitVec, tumorCellPos);
      updateSeqErrorStats(seqErrors, seqErrorsCombCov, tumor_counts);

      extractSeqInformation(normal_counts, splitVec, normalCellPos);
      updateSeqErrorStats(seqErrors, seqErrorsCombCov, normal_counts);
    } else
      lineNumber--;
  }

  // add pseudo counts
  seqErrors++;

  filePtr->rewind();

  return std::make_pair(seqErrors, seqErrorsCombCov);
}

// Check if the reference base is known.
bool isRefKnown(const std::string &n) {
  if (n[0] == 'n' || n[0] == 'N')
    return false;

  return true;
}

// This function is used to collect the sequencing information of all cells.
void extractSeqInformation(
        std::array<u_int32_t, 5> &count,
        const std::vector<std::string> &splitVec,
        u_int32_t position
) {
  count = {{0, 0, 0, 0, 0}}; // (a,c,g,t, coverage)
  count[4] = std::stoi(splitVec[position * 3 + 3]);
  extractCellNucCountInformation(count, splitVec[position * 3 + 4]);
}

void extractCellNucCountInformation(
        std::array<u_int32_t, 5> &counts,
        const std::string &nucleotides
) {
  //counts = {{0, 0, 0, 0, 0}}; // (a,c,g,t, coverage)
  for (size_t j = 0; j < nucleotides.size(); ++j) // loop over the nucleotides of a cell in the pileup
  {
    u_int16_t index = charToIndex(nucleotides[j]);
    if (index < 4)
      ++counts[index];
    else if (index == 9)
      --counts[4];
    else if (index == 4)
      j = skipIndels(nucleotides, j);
    else if (index == 5)
      ++j;
    else if (index == 6)
      continue;
  }
}

u_int16_t charToIndex(char c) {
  switch (std::toupper(c)) {
    case ('A'):
      return 0;
    case ('C'):
      return 1;
    case ('G'):
      return 2;
    case ('T'):
      return 3;
    case ('-'):
    case ('+'):
      return 4;
    case ('^'):
      return 5;
    case ('$'):
      return 6;
    case ('.'):
    case (','):
      return 7;
    case ('*'):
      return 8;
  }

  if (('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z')) {
    return 9;
  }

  std::cerr << "Unknown character \"" << c << "\" in pileup sequence!" << std::endl;
  std::exit(EXIT_FAILURE);

  return 10;
}

// This function is used to skip indels in the pileup.
u_int32_t skipIndels(
        const std::string &nucleotides,
        u_int32_t currentPos
) {
  if (nucleotides[currentPos] != '-' && nucleotides[currentPos] != '+')
    return currentPos;

  u_int32_t numIndels = 0;
  u_int32_t i = currentPos + 1; // skip the '-' or '+'

  for (; i < nucleotides.size(); ++i) {
    if (nucleotides[i] >= '0' && nucleotides[i] <= '9') {
      numIndels *= 10;
      numIndels += static_cast<int32_t>(nucleotides[i]) - 48;
    } else
      break;
  }
  return i + (numIndels - 1);
}

bool passNormalFilter(
        const std::array<u_int32_t, 5> &normalCounts,
        const Config &config
) {
  if (normalCounts[4] >= config.getMinCovInBulk()) {
    for (size_t i = 0; i < 4; ++i)
      if (normalCounts[i] >= config.getMaxSupInBulk())
        return false;

    return true;
  }

  return false;
}

bool passNormalCellCoverage(
        const std::vector<std::array<u_int32_t, 5>> &normalCellCounts,
        const Config &config
) {
  u_int32_t maxCov = 0;

  for (const auto &normalCellCount: normalCellCounts)
    if (normalCellCount[4] > maxCov)
      maxCov = normalCellCount[4];

  if (maxCov < config.getMinCovNormalCell())
    return false;

  return true;
}

void applyCoverageFilterPerCell(
        std::vector<std::array<u_int32_t, 5>> &counts,
        const Config &config
) {
  for (auto &count: counts)
    if (count[4] < config.getMinCov())
      for (size_t j = 0; j < 5; ++j)
        count[j] = 0;
}

bool applyFilterAcrossCells(
        const std::vector<std::array<u_int32_t, 5>> &counts,
        const Config &config,
        u_int32_t nucId
) {
  u_int32_t numCellsAboveThresh = 0;
  for (auto &count: counts) {
    if (passCovFilter(count[4], config.getMinCov()) &&
        passFreqFilter(count[nucId], count[4], config.getMinFreq()) &&
        passSuppFilter(count[nucId], config.getMinSup())) {
      ++numCellsAboveThresh;

      if (numCellsAboveThresh >= config.getMinTumorCellsToPass())
        return true;
    }
  }

  return false;
}

bool passSuppFilter(
        u_int32_t altCount,
        u_int32_t minSupport
) {
  if (altCount >= minSupport)
    return true;

  return false;
}

bool passFreqFilter(
        double altCount,
        double coverage,
        double minFreq
) {
  if (altCount / coverage >= minFreq)
    return true;

  return false;
}

bool passCovFilter(
        u_int32_t coverage,
        u_int32_t minCov
) {
  if (coverage >= minCov)
    return true;

  return false;
}

bool passNormalCellFilter(
        const std::vector<std::array<u_int32_t, 5>> &normalCellCounts,
        u_int16_t j,
        const Config &config
) {
  u_int32_t numMutated = 0;

  for (const auto &normalCellCount: normalCellCounts) {
    const auto homoRef = computeRawWildLogScore(config, normalCellCount[j], normalCellCount[4]);
    const auto heteroMu = computeRawHeteroMutLogScore(config, normalCellCount[j], normalCellCount[4]);
    const auto homoMu = computeRawHomoMutLogScore(config, normalCellCount[4] - normalCellCount[j], normalCellCount[4]);

    if (homoRef < heteroMu || homoRef < homoMu)
      ++numMutated;
  }

  if (numMutated > config.getMaxNormalCellsToMutate())
    return false;

  return true;
}

double updateLogH1Temp(
        double logH1Temp,
        u_int32_t numCells,
        u_int32_t numMuts,
        bool tumorCells,
        double dropOut
) {
  if (tumorCells)
    return logH1Temp + 2 * logNChooseK(numCells, numMuts) - std::log(2 * numMuts - 1) -
           logNChooseK(2 * numCells, 2 * numMuts);

  return logH1Temp + logNChooseK(numCells, numMuts) + std::log(std::pow(1.0 - dropOut, numMuts)) +
         std::log(std::pow(dropOut, numCells - numMuts)) - (1.0 - std::pow(dropOut, numCells));
}

bool mustH0Win(
        double &logH1Max,
        double logH1Temp,
        double logNumCells,
        double logH0
) {
  if (logH1Temp >= logH1Max) {
    logH1Max = logH1Temp;
    return false;
  } else {
    if (logH1Max + logNumCells < logH0)
      return true;
  }

  return false;
}

bool computeProbCellsAreMutated(
        const Config &config,
        std::vector<long double> &logProbs,
        std::vector<long double> &tempLogProbs,
        std::vector<double> &logProbTempValues,
        std::vector<std::array<u_int32_t, 5>> &filteredCounts,
        std::vector<double> &cellsNotMutated,
        std::vector<double> &cellsMutated,
        u_int32_t currentChar,
        bool tumorCells
) {
  u_int32_t numCells = filteredCounts.size();
  double logNumCells = log(numCells);

  tempLogProbs[0] = 0.0;
  for (size_t i = 1; i < filteredCounts.size() + 1; ++i) {
    cellsNotMutated[i - 1] = computeRawWildLogScore(config, filteredCounts[i - 1][currentChar],
                                                    filteredCounts[i - 1][4]);
    cellsMutated[i - 1] = logSumExp(
            computeRawHeteroMutLogScore(config, filteredCounts[i - 1][currentChar], filteredCounts[i - 1][4]),
            computeRawHomoMutLogScore(config, filteredCounts[i - 1][4] - filteredCounts[i - 1][currentChar],
                                      filteredCounts[i - 1][4]));
//    cellsMutated[i - 1] = computeRawHeteroMutLogScore(config, filteredCounts[i - 1][currentChar], filteredCounts[i - 1][4]);
    tempLogProbs[i] = tempLogProbs[i - 1] + cellsNotMutated[i - 1];
  }

  swap(logProbs, tempLogProbs);

  const double prior = tumorCells ? config.getMuRatePrior() : config.getGermlineRatePrior();

  const double logH0 = logProbs.back() + log(1.0 - prior);

  logProbTempValues[0] = logH0;

  // compute the probabilitues of observing 1, 2, 3, ... mutations
  double logH1Max = -DBL_MAX; // the current best alternative score
  long double logNOverK = 0;  // helper to efficiently compute nChooseK
  size_t numMut = 1;          // number of mutations currently computet

  for (; numMut <= numCells; ++numMut) {
    double logProbAllPrevCellsMutated = logProbs[numMut - 1];
    double currentCellMutated = cellsMutated[numMut - 1];
    tempLogProbs[numMut] = logProbAllPrevCellsMutated + currentCellMutated;
    for (size_t i = numMut + 1; i < filteredCounts.size() + 1; ++i) {
      double previousCellNotMutated = logProbs[i - 1];
      currentCellMutated = cellsMutated[i - 1];
      double previousCellMutated = tempLogProbs[i - 1];
      double currentCellNotMutated = cellsNotMutated[i - 1];
      tempLogProbs[i] = addLogProb(
              previousCellNotMutated + currentCellMutated,
              previousCellMutated + currentCellNotMutated
      );

    }
    swap(logProbs, tempLogProbs);
    logNOverK = logNChooseK(numCells, numMut, logNOverK);
    double logH1Temp = logProbs.back() + log(prior) - logNOverK;
    logH1Temp = updateLogH1Temp(logH1Temp, numCells, numMut, tumorCells, config.getAdoRatePrior() / 2.0);

    // check whether the alternative hypothesis can win
    bool h0Wins = mustH0Win(logH1Max, logH1Temp, logNumCells, logH0);
    logProbTempValues[numMut] = logH1Temp;
    if (h0Wins)
      return true;
  }

  return false;
}

double sumValuesInLogSpace(
        std::vector<double>::const_iterator itBegin,
        std::vector<double>::const_iterator itEnd
) {
  auto it = itBegin;
  double maxLogValue = *it;
  ++it;

  for (; it != itEnd; ++it)
    if (maxLogValue < *it)
      maxLogValue = *it;

  it = itBegin;
  double h1 = 0;

  for (; it != itEnd; ++it)
    h1 += exp(*it - maxLogValue);

  return log(h1) + maxLogValue;
}

char indexToChar(u_int16_t index) {
  switch (index) {
    case (0):
      return 'A';
    case (1):
      return 'C';
    case (2):
      return 'G';
    case (3):
      return 'T';
    default:
      return 'N';
  }
}

void writeAltNucInfo(
        const Data &data,
        const Config &config
) {
  if (config.getOutputFileName().has_parent_path()) {
    std::string mkdir =
            "mkdir -p " +
            static_cast<std::string>(config.getOutputFileName().parent_path());
    std::system(mkdir.c_str());
  }

  std::ofstream out;
  out.open(config.getOutputFileName());

  out << "=numSamples=" << "\n";
  out << config.getTumorCellNum() << "\n";

  out << "=numCandidateMutatedSites=" << "\n";
  out << data.getCandidateMutatedSitesNum() << "\n";

  if (data.getBackgroundSitesNum() % static_cast<u_int64_t>(config.getTumorCellNum()) != 0)
    std::cerr
            << "WARNING! Read counts data of some cells for background sites are missing. Make sure your mpileups are generated correctly.\n";

  out << "=numBackgroundSites=" << "\n";
  out << std::fixed;
  out.precision(0);
  out << std::floor(data.getBackgroundSitesNum() / static_cast<u_int64_t>(config.getTumorCellNum())) << "\n";

  out << "=mutations=" << "\n";
  for (const auto &i: data.getCandidateMutatedReadCounts()) {
    // chrom
    out << std::get<0>(std::get<0>(i)) << "\t";

    // pos
    out << std::get<1>(std::get<0>(i)) << "\t";

    // ref nuc
    out << std::get<2>(std::get<0>(i)) << "\t";

    // significant alt nucs
    if (std::get<3>(std::get<0>(i)).empty())
      out << "N";
    else {
      for (size_t j = 0; j < std::get<3>(std::get<0>(i)).size(); j++) {
        out << std::get<3>(std::get<0>(i))[j];

        if (j < std::get<3>(std::get<0>(i)).size() - 1)
          out << ',';
      }
    }

    // each cell
    for (size_t j = 0; j < config.getTumorCellNum(); j++)
      out << "\t" << std::get<1>(i)[j];

    out << "\n";
  }

  out << "=background=" << "\n";
  out << data.getBackgroundSitesReadCounts();

  out.close();
}

void process_single_file(
        const Config &config,
        std::unique_ptr<ReadCountsFile> &filePtr,
        const u_int32_t &tumorCellNum,
        const std::vector<u_int32_t> &tumorCellPos,
        const std::vector<u_int32_t> &normalCellPos,
        const u_int32_t &tumorBulkPos,
        const u_int32_t &normalBulkPos,
        const std::set<std::pair<ChromosomeLabel, u_int32_t>> &exMap,
        const std::set<std::pair<ChromosomeLabel, u_int32_t>> &incMap,
        const std::set<std::tuple<ChromosomeLabel, u_int32_t, char, char>> &incMuMap,
        std::mutex &_mutex,
        Data &data
) {
  const std::thread::id _id(std::this_thread::get_id());
  std::cout << "> [" << _id << "] Processing " + filePtr->getFileName() + "... " << std::endl;

  TData _data;
  std::vector<CellReadCounts> cellReadCounts{};
  cellReadCounts.reserve(tumorCellNum);
  SignificantAltNucs significantAltNucs{};
  std::array<u_int16_t, 3> altNucs{};
  u_int16_t altNucIdx{};
  ContinuousNoiseCounts continuousNoiseCounts{};
  bool positionKept, positionMutated, hasSignificantAltNucs;

  std::array<u_int32_t, 5> normalBulkCounts{};

  // vector to hold the nucleotide information {a,c,g,t,coverage}
  std::vector<std::array<u_int32_t, 5>> counts(tumorCellNum, {{0, 0, 0, 0, 0}});

  // vector to hold the nucleotide information {a,c,g,t,coverage}
  std::vector<std::array<u_int32_t, 5>> countsNormal(normalCellPos.size(), {{0, 0, 0, 0, 0}});

  // vector to hold the filtered nucleotide information {a,c,g,t,coverage}
  std::vector<std::array<u_int32_t, 5>> filteredCounts(tumorCellNum, {{0, 0, 0, 0, 0}});

  // probabilities of observing 0, 1, 2, 3, 4 ... mutations
  std::vector<long double> logProbsNormal(normalCellPos.size() + 1, 0);

  // helper array for probabilities of observing 0, 1, 2, 3, 4 ... mutations
  std::vector<long double> tempLogProbsNormal(normalCellPos.size() + 1, 0);

  std::vector<double> logProbTempValuesNormal(normalCellPos.size() + 1);

  // probabilities of observing 0, 1, 2, 3, 4 ... mutations
  std::vector<long double> logProbs(tumorCellNum + 1, 0);

  std::vector<double> logProbTempValues(tumorCellNum + 1);

  // helper array for probabilities of observing 0, 1, 2, 3, 4 ... mutations
  std::vector<long double> tempLogProbs(tumorCellNum + 1, 0);

  std::vector<double> cellsNotMutated(tumorCellNum);
  std::vector<double> cellsNotMutatedNormal(normalCellPos.size());
  std::vector<double> cellsMutated(tumorCellNum);
  std::vector<double> cellsMutatedNormal(normalCellPos.size());

  std::string currentLine;
  std::vector<std::string> splitVec;

  filePtr->setClearCache(true);

  while (filePtr->getLine(currentLine)) {
    if (currentLine.empty())
      continue;

    significantAltNucs.resetSigAltNucs();
    altNucIdx = 0;
    positionKept = false;
    positionMutated = false;
    hasSignificantAltNucs = false;
    cellReadCounts.clear();

    // split the current line into easily accessible chunks
    boost::split(splitVec, currentLine, boost::is_any_of("\t"));

    if (!isRefKnown(splitVec[2]))
      continue;

    // check if the current pos is to be excluded
    auto exIt = exMap.find(
            std::make_pair(
                    std::move(ChromosomeLabel(splitVec[0])),
                    std::stoi(splitVec[1])
            )
    );

    // check if the current pos is to be included
    auto incIt = incMap.find(
            std::make_pair(
                    std::move(ChromosomeLabel(splitVec[0])),
                    std::stoi(splitVec[1])
            )
    );

    if (exIt != exMap.end() && incIt != incMap.end()) {
      std::cerr << splitVec[0] << "'s " << splitVec[1] << " appears both in " <<
                config.getExcludedSitesFileName() << " and " <<
                config.getIncludedSitesFileName() << "." << std::endl;
      exit(1);
    } else if (exIt != exMap.end()) {
      continue;
    } else {
      if (incIt != incMap.end())
        positionKept = true;

      if (normalBulkPos != UINT_MAX) {
        extractSeqInformation(normalBulkCounts, splitVec, normalBulkPos);

        if (!passNormalFilter(normalBulkCounts, config))
          continue;
      }

      if (!normalCellPos.empty()) {
        extractSeqInformation(countsNormal, splitVec, normalCellPos);
        if (!passNormalCellCoverage(countsNormal, config))
          continue;
      }

      extractSeqInformation(counts, splitVec, tumorCellPos);

      filteredCounts = counts;
      applyCoverageFilterPerCell(filteredCounts, config);

      for (size_t altAlleleIdx = 0; altAlleleIdx < 4; ++altAlleleIdx) {
        if (altAlleleIdx == charToIndex(splitVec[2][0]))
          continue;

        altNucs[altNucIdx++] = altAlleleIdx;
        bool h0Wins = !applyFilterAcrossCells(filteredCounts, config, altAlleleIdx);

        if (!normalCellPos.empty()) {
          // use the normal cell filter
          if (config.getNormalCellFilterMode() == 1) {
            if (!passNormalCellFilter(countsNormal, altAlleleIdx, config)) {
              positionMutated = true;
              h0Wins = true;
              continue;
            }
          } else if (config.getNormalCellFilterMode() == 2) {
            bool h0WinsNormal = computeProbCellsAreMutated(
                    config,
                    logProbsNormal,
                    tempLogProbsNormal,
                    logProbTempValuesNormal,
                    countsNormal,
                    cellsNotMutatedNormal,
                    cellsMutatedNormal,
                    altAlleleIdx,
                    false
            );

            if (!h0WinsNormal) {
              double logH0Normal = sumValuesInLogSpace(
                      logProbTempValuesNormal.begin(),
                      logProbTempValuesNormal.begin() +
                      config.getMaxNormalCellsToMutate() + 1);

              double logH1Normal = sumValuesInLogSpace(
                      logProbTempValuesNormal.begin() +
                      config.getMaxNormalCellsToMutate() + 1,
                      logProbTempValuesNormal.end());

              if (logH0Normal < logH1Normal) {
                positionMutated = true;
                h0Wins = true;
                continue;
              }
            }
          }
        }

        double logH0 = -DBL_MAX;
        double logH1 = -DBL_MAX;
        if (!h0Wins)
          h0Wins = computeProbCellsAreMutated(
                  config, logProbs, tempLogProbs, logProbTempValues,
                  filteredCounts, cellsNotMutated, cellsMutated, altAlleleIdx,
                  true);

        if (h0Wins) {
          logH1 = -DBL_MAX;
          logH0 = DBL_MAX;
        } else {
          logH0 =
                  sumValuesInLogSpace(
                          logProbTempValues.begin(),
                          logProbTempValues.begin() + config.getMinTumorCellsToCallMu()
                  );
          logH1 = sumValuesInLogSpace(
                  logProbTempValues.begin() + config.getMinTumorCellsToCallMu(),
                  logProbTempValues.end()
          );
        }

        if (logH1 > logH0 ||
            incMuMap.count(
                    std::make_tuple(
                            std::move(ChromosomeLabel(splitVec[0])),
                            std::stoi(splitVec[1]),
                            splitVec[2][0],
                            indexToChar(altAlleleIdx)
                    )
            ) != 0
                ) {
          positionMutated = true;

          std::vector<std::pair<u_int32_t, u_int32_t>> testCounts;
          for (size_t cell = 0; cell < counts.size(); ++cell)
            if (cellsNotMutated[cell] < cellsMutated[cell])
              testCounts.emplace_back(counts[cell][altAlleleIdx], counts[cell][4]);

          dlib::matrix<double, 0, 1> startingPointMeanOverDis = {0.05, 5.0};
          dlib::matrix<double, 0, 1> dLibMinMeanOverDis = {0.001, 0.1};
          dlib::matrix<double, 0, 1> dLibMaxMeanOverDis = {0.999, 10000.0};
          OptimizeBetaBinMeanOverDis optBetaBinMeanOverDis(testCounts);
          OptimizeBetaBinMeanOverDisDerivates optBetaBinDerMeanOverDis(
                  testCounts);
          double resultMeanOverDis = dlib::find_max_box_constrained(
                  dlib::bfgs_search_strategy(), // Use BFGS search algorithm
                  dlib::objective_delta_stop_strategy(
                          1e-7), // Stop when the change in rosen() is less than 1e-7
                  optBetaBinMeanOverDis, optBetaBinDerMeanOverDis,
                  startingPointMeanOverDis, dLibMinMeanOverDis,
                  dLibMaxMeanOverDis);

          dlib::matrix<double, 0, 1> startingPointOverDis = {2.0};
          dlib::matrix<double, 0, 1> dLibMinOverDis = {0.1};
          dlib::matrix<double, 0, 1> dLibMaxOverDis = {10000.0};
          OptimizeBetaBinOverDis optBetaBinOverDis(
                  testCounts, config.getMeanFreqSite());
          OptimizeBetaBinOverDisDerivates optBetaBinDerOverDis(
                  testCounts, config.getMeanFreqSite());
          double resultOverDis = dlib::find_max_box_constrained(
                  dlib::bfgs_search_strategy(), // Use BFGS search algorithm
                  dlib::objective_delta_stop_strategy(
                          1e-7), // Stop when the change in rosen() is less than 1e-7
                  optBetaBinOverDis, optBetaBinDerOverDis, startingPointOverDis,
                  dLibMinOverDis, dLibMaxOverDis);

          double pValue;
          if (resultOverDis >
              resultMeanOverDis) // if true the results are within the optimization limit
            pValue = 1;
          else {
            double testStat = -2 * (resultOverDis - resultMeanOverDis);
            boost::math::chi_squared mydist(1);
            pValue = 1 - boost::math::cdf(mydist, testStat);
          }

          if (pValue > config.getThreshold() ||
              startingPointMeanOverDis(0) >= config.getMeanFreqSite() ||
              incMuMap.count(
                      std::make_tuple(std::move(ChromosomeLabel(splitVec[0])),
                                      std::stoi(splitVec[1]), splitVec[2][0],
                                      indexToChar(altAlleleIdx))
              ) != 0) {
            hasSignificantAltNucs = true;
            significantAltNucs.addSigAltNuc(SignificantAltNuc(
                    altAlleleIdx, startingPointMeanOverDis(0), pValue));
          }
        }
      }

      if (!positionMutated)
        for (size_t cell = 0; cell < tumorCellNum; cell++)
          continuousNoiseCounts.add(std::array<u_int32_t, 4>{
                  counts[cell][altNucs[0]], counts[cell][altNucs[1]],
                  counts[cell][altNucs[2]], counts[cell][4]});

      // only treat the site as a candidate snv if it contains significant alternative nucleotides, or if it must be kept
      if (hasSignificantAltNucs || positionKept) {
        // sort significant alternative nucleotides
        significantAltNucs.getOrderedSigAltNucs();

        // collect read counts for each alternative nucleotide in each cell
        for (size_t cell = 0; cell < tumorCellNum; cell++) {
          cellReadCounts.emplace_back(
                  AltNucReadCount(altNucs[0], indexToChar(altNucs[0]),
                                  counts[cell][altNucs[0]], significantAltNucs),
                  AltNucReadCount(altNucs[1], indexToChar(altNucs[1]),
                                  counts[cell][altNucs[1]], significantAltNucs),
                  AltNucReadCount(altNucs[2], indexToChar(altNucs[2]),
                                  counts[cell][altNucs[2]], significantAltNucs),
                  counts[cell][4]);

          cellReadCounts[cell].sortAltNucReadCounts();
        }

        // save the data
        _data.emplace_back(std::move(TDataEntry(
                std::move(
                        TPosInfo(std::move(ChromosomeLabel(std::move(splitVec[0]))),
                                 std::stoi(splitVec[1]), splitVec[2][0],
                                 std::move(significantAltNucs.convertAltNucsType(
                                         indexToChar)))),
                std::move(cellReadCounts))));
      }
    }
  }

  {
    std::lock_guard<std::mutex> lock(_mutex);
    data.addReadCounts(std::move(_data));
    data.addContinuousNoiseCounts(continuousNoiseCounts);

    std::cout << "> [" << _id << "] Done." << std::endl;
  }
}
