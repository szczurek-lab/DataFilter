//
// Created by senbaikang on 14.02.21.
//

#ifndef DATAFILTER_OUTPUT_DATA_H
#define DATAFILTER_OUTPUT_DATA_H

#include <algorithm>
#include <array>
#include <iostream>
#include <regex>

#define CHROMOSOME_PATTERN "^([a-zA-Z]*)([1-9]|1\\d|2[012]|[XYxy]|MT|mt)$"
#define LABEL_POSITION 2


struct SignificantAltNuc {

private:
    u_int16_t sigAltNuc{0};
    double mean{0.0};
    double pValue{0.0};

public:
    SignificantAltNuc() = default;

    SignificantAltNuc(const SignificantAltNuc &) = default;

    SignificantAltNuc(SignificantAltNuc &&) noexcept = default;

    explicit SignificantAltNuc(
            u_int16_t,
            double,
            double
    );

    SignificantAltNuc &operator=(const SignificantAltNuc &) = default;

    SignificantAltNuc &operator=(SignificantAltNuc &&) noexcept = default;

    bool operator>(const SignificantAltNuc &v) const;

    [[nodiscard]] u_int16_t getSigAltNuc() const;

    [[nodiscard]] double getMean() const;

    [[nodiscard]] double getPValue() const;

};


struct SignificantAltNucs {

private:
    std::vector<SignificantAltNuc> sigAltNucs{};
    std::vector<u_int16_t> orderedSigAltNucs{};

public:
    SignificantAltNucs() = default;

    SignificantAltNucs(const SignificantAltNucs &) = default;

    SignificantAltNucs(SignificantAltNucs &&) noexcept;

    SignificantAltNucs &operator=(const SignificantAltNucs &) = default;

    SignificantAltNucs &operator=(SignificantAltNucs &&) noexcept;

    void addSigAltNuc(const SignificantAltNuc &);

    void addSigAltNuc(SignificantAltNuc &&);

    void getOrderedSigAltNucs();

    [[nodiscard]] int16_t getAltNucOrder(const u_int16_t &) const;

    void resetSigAltNucs();

    template<class T, class P>
    std::vector<T> convertAltNucsType(T (*)(P));

};

template<class T, class P>
std::vector<T> SignificantAltNucs::convertAltNucsType(T (*func)(P)) {
  std::vector<T> results(orderedSigAltNucs.size());

  for (size_t i = 0; i != orderedSigAltNucs.size(); i++)
    results[i] = (*func)(orderedSigAltNucs[i]);

  return results;
}

struct AltNucReadCount {

private:
    u_int16_t altNuc{0};
    char altNucChar{0};
    u_int32_t readCount{0};
    SignificantAltNucs sigAltNucs{};

public:
    AltNucReadCount() = default;

    AltNucReadCount(const AltNucReadCount &) = default;

    AltNucReadCount(AltNucReadCount &&) noexcept;

    explicit AltNucReadCount(
            u_int16_t,
            char,
            u_int32_t,
            SignificantAltNucs
    );

    AltNucReadCount &operator=(const AltNucReadCount &) = default;

    AltNucReadCount &operator=(AltNucReadCount &&) noexcept;

    bool operator>(const AltNucReadCount &v) const;

    [[nodiscard]] u_int16_t getAltNuc() const;

    [[nodiscard]] char getAltNucChar() const;

    [[nodiscard]] u_int32_t getReadCount() const;

    [[nodiscard]] const SignificantAltNucs &getSigAltNucs() const;
};


struct CellReadCounts {

private:
    std::array<AltNucReadCount, 3> altNucReadCounts{};
    u_int32_t cov{0};

public:
    CellReadCounts() = default;

    CellReadCounts(const CellReadCounts &) = default;

    CellReadCounts(CellReadCounts &&) noexcept;

    CellReadCounts(
            const AltNucReadCount &,
            const AltNucReadCount &,
            const AltNucReadCount &,
            u_int32_t
    );

    CellReadCounts(
            AltNucReadCount &&,
            AltNucReadCount &&,
            AltNucReadCount &&,
            u_int32_t
    ) noexcept;

    CellReadCounts &operator=(const CellReadCounts &) = default;

    CellReadCounts &operator=(CellReadCounts &&) noexcept;

    void sortAltNucReadCounts();

    friend std::ostream &operator<<(std::ostream &, const CellReadCounts &);

};


struct ChromosomeLabel {

private:
    static const std::regex pattern;

    std::string fullLabel; // e.g., chr3
    std::string label; // e.g., 3

    static std::string getLabel(const std::string &);

public:
    ChromosomeLabel() = delete;

    ChromosomeLabel(const ChromosomeLabel &) = default;

    ChromosomeLabel(ChromosomeLabel &&) noexcept;

    explicit ChromosomeLabel(const std::string &);

    explicit ChromosomeLabel(std::string &&) noexcept;

    ChromosomeLabel &operator=(const ChromosomeLabel &) = default;

    ChromosomeLabel &operator=(ChromosomeLabel &&) noexcept;

    bool operator==(const ChromosomeLabel &) const;

    bool operator!=(const ChromosomeLabel &) const;

    bool operator>(const ChromosomeLabel &) const;

    bool operator<(const ChromosomeLabel &) const;

    [[nodiscard]] std::string getFullLabel() const;

    friend std::ostream &operator<<(std::ostream &, const ChromosomeLabel &);

};


/**
 * For each cell at each background site, m1, m2, and m3 should be in descending order.
 */
struct ContinuousNoiseCounts {

private:
    std::vector<u_int64_t> m1{};
    std::vector<u_int64_t> m2{};
    std::vector<u_int64_t> m3{};
    std::vector<u_int64_t> ref{};
    std::vector<u_int64_t> cov{};

    static void add(std::vector<u_int64_t> &, size_t);

    static void alignSize(std::vector<u_int64_t> &, std::vector<u_int64_t> &);

    void alignSize(ContinuousNoiseCounts &);

    static void add(std::vector<u_int64_t> &, const std::vector<u_int64_t> &);

public:
    ContinuousNoiseCounts() = default;

    void add(std::array<u_int32_t, 4> &&);

    ContinuousNoiseCounts &operator+=(ContinuousNoiseCounts &);

    [[nodiscard]] const std::vector<u_int64_t> &getM1() const;

    [[nodiscard]] const std::vector<u_int64_t> &getM2() const;

    [[nodiscard]] const std::vector<u_int64_t> &getM3() const;

    [[nodiscard]] const std::vector<u_int64_t> &getRef() const;

    [[nodiscard]] const std::vector<u_int64_t> &getCov() const;
};


struct NoiseCounts {

private:
    std::vector<std::pair<size_t, u_int64_t>> m1{};
    std::vector<std::pair<size_t, u_int64_t>> m2{};
    std::vector<std::pair<size_t, u_int64_t>> m3{};
    std::vector<std::pair<size_t, u_int64_t>> ref{};
    std::vector<std::pair<size_t, u_int64_t>> cov{};
    u_int64_t numPos{};

    static std::ostream &printHelper(
            std::ostream &,
            const std::vector<std::pair<size_t, u_int64_t>> &
    );

public:
    NoiseCounts() = default;

    explicit NoiseCounts(const ContinuousNoiseCounts &);

    NoiseCounts(const NoiseCounts &) = default;

    NoiseCounts(NoiseCounts &&) noexcept;

    NoiseCounts &operator=(const NoiseCounts &) = default;

    NoiseCounts &operator=(NoiseCounts &&) noexcept;

    static void collectCounts(
            std::vector<std::pair<size_t, u_int64_t>> &,
            const std::vector<u_int64_t> &,
            u_int64_t *
    );

    [[nodiscard]] u_int64_t getBackgroundSitesNum() const;

    friend std::ostream &operator<<(std::ostream &, const NoiseCounts &);

};


typedef std::tuple<ChromosomeLabel, u_int32_t, char, std::vector<char>> TPosInfo; // chromosome label, position, reference nucleotide, alternative nucleotides
typedef std::pair<TPosInfo, std::vector<CellReadCounts>> TDataEntry;
typedef std::vector<TDataEntry> TData;


bool compareDataEntry(const TDataEntry &, const TDataEntry &);


/**
 * Collection of output data.
 */
class Data {

private:
    std::vector<std::string> cellNames{};
    TData readCounts{};
    ContinuousNoiseCounts continuousNoiseCounts{};
    NoiseCounts noiseCounts{};

public:
    Data() = default;

    void addCellName(const std::string &);

    void addCellName(std::string &&);

    void addReadCounts(
            const std::vector<std::pair<std::tuple<ChromosomeLabel, unsigned int, char, std::vector<char>>, std::vector<CellReadCounts>>> &);

    void addReadCounts(
            std::vector<std::pair<std::tuple<ChromosomeLabel, unsigned int, char, std::vector<char>>, std::vector<CellReadCounts>>> &&);

    void addContinuousNoiseCounts(ContinuousNoiseCounts &);

    void createNoiseCounts();

    void sortEntries();

    [[nodiscard]] const std::vector<std::string> &getCellNames() const;

    [[nodiscard]] const TData &getCandidateMutatedReadCounts() const;

    [[nodiscard]] u_int32_t getCandidateMutatedSitesNum() const;

    [[nodiscard]] u_int64_t getBackgroundSitesNum() const;

    [[nodiscard]] const NoiseCounts &getBackgroundSitesReadCounts() const;
};

#endif // DATAFILTER_OUTPUT_DATA_H
