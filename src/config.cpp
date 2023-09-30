//
// Created by senbaikang on 15.05.21.
//

#include <boost/program_options.hpp>
#include <fstream>

#include <datafilter/config.h>

Config::Config(
    u_int16_t threads,
    std::string cellNamesFileName,
    std::vector<std::string> inputFileNames,
    std::string outputFileName
    ):
        threads(threads),
        cellNamesFileName(std::move(cellNamesFileName)),
        inputFileNames(std::move(inputFileNames)),
        outputFileName(std::move(static_cast<std::filesystem::path>(std::move(outputFileName))))
{}

Config::Config(
    u_int16_t threads,
    std::string cellNamesFileName,
    std::string cellNameSuffix,
    std::vector<std::string> inputFileNames,
    u_int32_t readLength,
    std::string outputFileName,
    std::string excludedSitesFileName,
    std::string excludedMutatedSitesFileName,
    std::string includedSitesFileName,
    double muRatePrior,
    double germlineRatePrior,
    double effectiveSeqErrRate,
    double wildOverdispersion,
    double muOverdispersion,
    bool estimateErrRate,
    u_int32_t errRateEstLoops,
    u_int32_t minTumorCellsToCallMu,
    u_int32_t minTumorCellsToPass,
    u_int32_t minCov,
    u_int32_t minSup,
    double minFreq,
    double meanFreqSite,
    u_int32_t maxSupInBulk,
    u_int32_t minCovInBulk,
    u_int16_t normalCellFilterMode,
    u_int32_t maxNormalCellsToMutate
    ):
        threads(threads),
        cellNamesFileName(std::move(cellNamesFileName)),
        cellNameSuffix(std::move(cellNameSuffix)),
        inputFileNames(std::move(inputFileNames)),
        readLength(readLength * 16),
        outputFileName(std::move(static_cast<std::filesystem::path>(std::move(outputFileName)))),
        excludedSitesFileName(std::move(excludedSitesFileName)),
        excludedMutatedSitesFileName(std::move(excludedMutatedSitesFileName)),
        includedSitesFileName(std::move(includedSitesFileName)),
        muRatePrior(muRatePrior),
        germlineRatePrior(germlineRatePrior),
        effectiveSeqErrRate(effectiveSeqErrRate),
        wildOverdispersion(wildOverdispersion),
        muOverdispersion(muOverdispersion),
        estimateErrRate(estimateErrRate),
        errRateEstLoops(errRateEstLoops),
        minTumorCellsToCallMu(minTumorCellsToCallMu),
        minTumorCellsToPass(minTumorCellsToPass),
        minCov(minCov),
        minSup(minSup),
        minFreq(minFreq),
        meanFreqSite(meanFreqSite),
        maxSupInBulk(maxSupInBulk),
        minCovInBulk(minCovInBulk),
        normalCellFilterMode(normalCellFilterMode),
        maxNormalCellsToMutate(maxNormalCellsToMutate)
{}

u_int16_t Config::getThreads() const { return threads; }

const std::string &Config::getCellNamesFileName() const { return cellNamesFileName; }

const std::string &Config::getCellNameSuffix() const { return cellNameSuffix; }

const std::string &Config::getExcludedSitesFileName() const { return excludedSitesFileName; }

const std::string &Config::getExcludedMutatedSitesFileName() const { return excludedMutatedSitesFileName; }

const std::string &Config::getIncludedSitesFileName() const { return includedSitesFileName; }

const std::string &Config::getIncludedMutatedSitesFileName() const { return includedMutatedSitesFileName; }

const std::vector<std::string> &Config::getInputFileNames() const { return inputFileNames; }

u_int64_t Config::getReadLength() const { return readLength; }

const std::filesystem::path &Config::getOutputFileName() const { return outputFileName; }

double Config::getMuRatePrior() const { return muRatePrior; }

double Config::getGermlineRatePrior() const { return germlineRatePrior; }

double Config::getEffectiveSeqErrRate() const { return effectiveSeqErrRate; }

void Config::setEffectiveSeqErrRate(double v) { effectiveSeqErrRate = v; }

double Config::getWildOverdispersion() const { return wildOverdispersion; }

double Config::getMuOverdispersion() const { return muOverdispersion; }

double Config::getAdoRatePrior() const { return adoRatePrior; }

bool Config::isEstimateErrRate() const { return estimateErrRate; }

u_int32_t Config::getErrRateEstLoops() const { return errRateEstLoops; }

u_int32_t Config::getMinTumorCellsToCallMu() const { return minTumorCellsToCallMu; }

u_int32_t Config::getMinTumorCellsToPass() const { return minTumorCellsToPass; }

u_int32_t Config::getMinCovNormalCell() const { return minCovNormalCell; }

u_int32_t Config::getMinCov() const { return minCov; }

u_int32_t Config::getMinSup() const { return minSup; }

double Config::getMinFreq() const { return minFreq; }

double Config::getMeanFreqSite() const { return meanFreqSite; }

u_int32_t Config::getMaxSupInBulk() const { return maxSupInBulk; }

u_int32_t Config::getMinCovInBulk() const { return minCovInBulk; }

u_int16_t Config::getNormalCellFilterMode() const { return normalCellFilterMode; }

u_int32_t Config::getMaxNormalCellsToMutate() const { return maxNormalCellsToMutate; }

u_int32_t Config::getTumorCellNum() const { return tumorCellNum; }

void Config::setTumorCellNum(u_int32_t v) { tumorCellNum = v; }

int parseCommandLineArgs(
    Config &config,
    int argc,
    char* argv[]
)
{
  boost::program_options::options_description generic("Generic options");
  generic.add_options()("help,h", "Print help message.");

  boost::program_options::options_description required("Required options");
  required.add_options()
      (
          "cellNames",
          boost::program_options::value<decltype(config.cellNamesFileName)>(&config.cellNamesFileName)->required(),
          "File name of names of the BAM files used to create the mpileup."
      )
      (
          "out,o",
          boost::program_options::value<decltype(config.outputFileName)>(&config.outputFileName)->required(),
          "Output file name."
      )
      ;

  boost::program_options::options_description mutex("Required but mutually exclusive options");
  mutex.add_options()
      (
          "inFile",
          boost::program_options::value<std::string>(),
          "File name of new line separated paths to input mpileup files."
      )
      (
          "in",
          boost::program_options::value<decltype(config.inputFileNames)>(&config.inputFileNames)->multitoken(),
          "Multiple input mpileup files."
      )
      ;

  boost::program_options::options_description optional("Optional options");
  optional.add_options()
      (
          "threads,t",
          boost::program_options::value<decltype(config.threads)>(&config.threads),
          "Number of threads. [1]"
      )
      (
          "cellNameSuf",
          boost::program_options::value<decltype(config.cellNameSuffix)>(&config.cellNameSuffix),
          "Suffix in regex form of cell names to remove. [\\\\..*?$]"
      )
      (
          "ex",
          boost::program_options::value<decltype(config.excludedSitesFileName)>(&config.excludedSitesFileName),
          "File name of exclusion list (VCF format), containing loci which should be ignored."
      )
      (
          "me",
          boost::program_options::value<decltype(config.excludedMutatedSitesFileName)>(&config.excludedMutatedSitesFileName),
          "File name of mutations to exclude during the sequencing error rate estimation (VCF format)."
      )
      (
          "inc",
          boost::program_options::value<decltype(config.includedSitesFileName)>(&config.includedSitesFileName),
          "File name of inclusion list (VCF format) containing Variants (CHROM, POS) that should be included."
      )
      (
          "mi",
          boost::program_options::value<decltype(config.includedMutatedSitesFileName)>(&config.includedMutatedSitesFileName),
          "File name of inclusion list (VCF format) containing Variants (CHROM, POS, REF, ALT) that should be included."
      )
      (
          "readLength",
          boost::program_options::value<decltype(config.readLength)>(&config.readLength),
          "Read length of gunzipped input mpileup. [524288]"
      )
      (
          "mupr",
          boost::program_options::value<decltype(config.muRatePrior)>(&config.muRatePrior),
          "Prior mutation rate. [0.0001]"
      )
      (
          "germpr",
          boost::program_options::value<decltype(config.germlineRatePrior)>(&config.germlineRatePrior),
          "Prior germline rate. [0.001]"
      )
      (
          "adopr",
          boost::program_options::value<decltype(config.adoRatePrior)>(&config.adoRatePrior),
          "Prior allelic drop-out rate. [0.1]"
      )
      (
          "wildMean",
          boost::program_options::value<decltype(config.effectiveSeqErrRate)>(&config.effectiveSeqErrRate),
          "Mean error rate. If the effective sequencing error rate should not be learned \"--ese 0\" one can specify it. [0.001]"
      )
      (
          "wildOverDis",
          boost::program_options::value<decltype(config.wildOverdispersion)>(&config.wildOverdispersion),
          "Initial overdispersion for wild type. [100.0]"
      )
      (
          "muOverDis",
          boost::program_options::value<decltype(config.muOverdispersion)>(&config.muOverdispersion),
          "Initial overdispersion for mutant type. [2.0]"
      )
      (
          "ese",
          boost::program_options::value<decltype(config.estimateErrRate)>(&config.estimateErrRate),
          "Estimate the effective sequencing error rate. [on]"
      )
      (
          "maxese",
          boost::program_options::value<decltype(config.errRateEstLoops)>(&config.errRateEstLoops),
          "Max number of sites per input file for estimating effective sequencing error rate. [100000]"
      )
      (
          "cwm",
          boost::program_options::value<decltype(config.minTumorCellsToCallMu)>(&config.minTumorCellsToCallMu),
          "Number of tumor cells required to have a mutation in order to be called. [2]"
      )
      (
          "mnp",
          boost::program_options::value<decltype(config.minTumorCellsToPass)>(&config.minTumorCellsToPass),
          "Number of cells which need to pass the filters described below. [2]"
      )
      (
          "mcn",
          boost::program_options::value<decltype(config.minCovNormalCell)>(&config.minCovNormalCell),
          "Minimum coverage required per normal cell. [5]"
      )
      (
          "mc",
          boost::program_options::value<decltype(config.minCov)>(&config.minCov),
          "Minimum coverage required per cell. [1]"
      )
      (
          "ms",
          boost::program_options::value<decltype(config.minSup)>(&config.minSup),
          "Minimum number of reads required to support the alternative. [3]"
      )
      (
          "mf",
          boost::program_options::value<decltype(config.minFreq)>(&config.minFreq),
          "Minimum required frequency of reads supporting the alternative per cell. [0.0]"
      )
      (
          "mff",
          boost::program_options::value<decltype(config.meanFreqSite)>(&config.meanFreqSite),
          "Mean of acceptable variant allele frequency across all cells for a specific locus. Mapping artifacts may result in low allele frequencies across cells. In order to filter these events out we apply a log-likelihood ratio test where the sequencing error model has a mean of this value. [0.25]"
      )
      (
          "bns",
          boost::program_options::value<decltype(config.maxSupInBulk)>(&config.maxSupInBulk),
          "Loci with up to this number of alternative supporting reads in the bulk control sample will be skipped as germline. [2]"
      )
      (
          "bnc",
          boost::program_options::value<decltype(config.minCovInBulk)>(&config.minCovInBulk),
          "Minimum required coverage of reads in the bulk control sample. [6]"
      )
      (
          "ncf",
          boost::program_options::value<decltype(config.normalCellFilterMode)>(&config.normalCellFilterMode),
          "Normal cell filter. Currently there are three options: (0) Do not use the normal cells for filtering; (1) use a simple filtering scheme excluding mutations if the probability of being mutated is higher than not being mutated for any cell independently; (2) filter mutations where the probability that at least one cell is mutated is higher than no cell is mutated. Note that in contrast to (1) the cells are not independent and cells with no alternative support need to be explained via dropout events. [1]"
      )
      (
          "mnc",
          boost::program_options::value<decltype(config.maxNormalCellsToMutate)>(&config.maxNormalCellsToMutate),
          "Maximum number of control cells allowed to be mutated. [0]"
      )
      ;

  boost::program_options::options_description allOptions;
  allOptions.add(generic).add(required).add(mutex).add(optional);

  boost::program_options::variables_map vm;

  try
  {
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(allOptions).run(), vm);

    // show help options
    if (vm.count("help"))
    {
      std::cout << allOptions << std::endl;
      exit(EXIT_SUCCESS);
    }

    boost::program_options::notify(vm);
  }
  catch (boost::program_options::required_option& e) {
    if (e.get_option_name() == "--cellNames")
      std::cerr << "Error! Missing option: --cellNames." << std::endl;
    else if (e.get_option_name() == "--out" || e.get_option_name() == "-o")
      std::cerr << "Error! Missing option: --out or -o." << std::endl;
    else
      std::cerr << "Error! " << e.what() << std::endl;

    exit(EXIT_FAILURE);
  }
  catch (boost::program_options::error& e)
  {
    std::cerr << "Error! " << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((vm.count("inFile") && vm.count("in")) || (!vm.count("inFile") && !vm.count("in")))
  {
    std::cerr << "Error! only one of --inFile and --in can be defined." << std::endl;
    exit(EXIT_FAILURE);
  }
  else if (vm.count("inFile"))
  {
    std::string line;
    std::ifstream input(vm["inFile"].as<std::string>(), std::ifstream::in);
    while (getline(input, line))
    {
      if (std::find(config.inputFileNames.begin(), config.inputFileNames.end(), line) == config.inputFileNames.end())
        config.inputFileNames.emplace_back(line);
    }
  }

  return 0;
}
