//
// Created by senbaikang on 14.05.21.
//

#ifndef DATAFILTER_PROCESS_DATA_H
#define DATAFILTER_PROCESS_DATA_H

#include <set>
#include <vector>

#include <datafilter/config.h>
#include <datafilter/read_counts_file.h>

#define GZ_FILE_PATTERN ".*?\\.gz$"

void process_data(Config &);

void readCellInformation(
        std::vector<u_int32_t> &,
        std::vector<u_int32_t> &,
        u_int32_t &,
        u_int32_t &,
        const Config &,
        Data &
);

void readCellNames(const Config &, Data &);

void readSpecialSites(
        const std::string &,
        std::set<std::pair<ChromosomeLabel, u_int32_t>> &
);

void readInclusionVCF(
        const std::string &,
        std::set<std::tuple<ChromosomeLabel, u_int32_t, char, char>> &
);

std::vector<std::unique_ptr<ReadCountsFile>> getReadCountsFiles(const Config &);

void closeReadCountsFile(std::vector<std::unique_ptr<ReadCountsFile>> &);

void extractSeqInformation(
        std::vector<std::array<u_int32_t, 5>> &,
        const std::vector<std::string> &,
        const std::vector<u_int32_t> &
);

void updateSeqErrorStats(
        u_int32_t &,
        u_int32_t &,
        const std::vector<std::array<u_int32_t, 5>> &
);

std::pair<u_int32_t, u_int32_t> estimateSeqErrorRate(
        const Config &,
        std::unique_ptr<ReadCountsFile> &,
        const std::set<std::pair<ChromosomeLabel, u_int32_t>> &,
        const std::set<std::pair<ChromosomeLabel, u_int32_t>> &,
        const std::vector<u_int32_t> &,
        const std::vector<u_int32_t> &
);

bool isRefKnown(const std::string &);

void extractSeqInformation(
        std::array<u_int32_t, 5> &,
        const std::vector<std::string> &,
        u_int32_t position
);

void extractCellNucCountInformation(
        std::array<u_int32_t, 5> &,
        const std::string &
);

u_int16_t charToIndex(char);

u_int32_t skipIndels(
        const std::string &,
        u_int32_t
);

bool passNormalFilter(
        const std::array<u_int32_t, 5> &,
        const Config &
);

bool passNormalCellCoverage(
        const std::vector<std::array<u_int32_t, 5>> &,
        const Config &
);

void applyCoverageFilterPerCell(
        std::vector<std::array<u_int32_t, 5>> &,
        const Config &
);

bool applyFilterAcrossCells(
        const std::vector<std::array<u_int32_t, 5>> &,
        const Config &,
        u_int32_t
);

bool passSuppFilter(
        u_int32_t,
        u_int32_t
);

bool passFreqFilter(
        double,
        double,
        double
);

bool passCovFilter(
        u_int32_t,
        u_int32_t
);

bool passNormalCellFilter(
        std::vector<std::array<u_int32_t, 5>> const &,
        u_int16_t,
        Config const &
);

double updateLogH1Temp(
        double,
        u_int32_t,
        u_int32_t,
        bool,
        double
);

bool mustH0Win(
        double &,
        double,
        double,
        double
);

bool computeProbCellsAreMutated(
        const Config &,
        std::vector<long double> &,
        std::vector<long double> &,
        std::vector<double> &,
        std::vector<std::array<u_int32_t, 5>> &,
        std::vector<double> &,
        std::vector<double> &,
        u_int32_t,
        bool
);

double sumValuesInLogSpace(
        std::vector<double>::const_iterator,
        std::vector<double>::const_iterator
);

char indexToChar(u_int16_t);

void writeAltNucInfo(
        const Data &,
        const Config &
);

void process_single_file(
        const Config &,
        std::unique_ptr<ReadCountsFile> &,
        const u_int32_t &,
        const std::vector<u_int32_t> &,
        const std::vector<u_int32_t> &,
        const u_int32_t &,
        const u_int32_t &,
        const std::set<std::pair<ChromosomeLabel, u_int32_t>> &,
        const std::set<std::pair<ChromosomeLabel, u_int32_t>> &,
        const std::set<std::tuple<ChromosomeLabel, u_int32_t, char, char>> &,
        std::mutex &,
        Data &
);

#endif // DATAFILTER_PROCESS_DATA_H
