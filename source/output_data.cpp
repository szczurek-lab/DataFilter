//
// Created by senbaikang on 15.05.21.
//

#include "output_data.h"

SignificantAltNuc::SignificantAltNuc(
    u_int16_t sigAltNuc,
    double mean,
    double pValue
    ):
    sigAltNuc(sigAltNuc),
    mean(mean),
    pValue(pValue)
{}

bool SignificantAltNuc::operator > (const SignificantAltNuc &v) const {
  if (getMean() != v.getMean())
    return getMean() > v.getMean();

  if (getPValue() != v.getPValue())
    return getPValue() > v.getPValue();

  return getSigAltNuc() < v.getSigAltNuc();
}

u_int16_t SignificantAltNuc::getSigAltNuc() const { return sigAltNuc; }

double SignificantAltNuc::getMean() const { return mean; }

double SignificantAltNuc::getPValue() const { return pValue; }

SignificantAltNucs::SignificantAltNucs(SignificantAltNucs &&v) noexcept:
    sigAltNucs(std::move(v.sigAltNucs)),
    orderedSigAltNucs(std::move(v.orderedSigAltNucs))
{}

SignificantAltNucs & SignificantAltNucs::operator = (SignificantAltNucs &&v) noexcept
{
  if (this != &v)
  {
    sigAltNucs = std::move(v.sigAltNucs);
    orderedSigAltNucs = std::move(v.orderedSigAltNucs);
  }

  return *this;
}

void SignificantAltNucs::addSigAltNuc(const SignificantAltNuc & v) { sigAltNucs.emplace_back(v); }

void SignificantAltNucs::addSigAltNuc(SignificantAltNuc &&v) { sigAltNucs.emplace_back(v); }

void SignificantAltNucs::getOrderedSigAltNucs()
{
  if (this->sigAltNucs.empty()) return;

  std::sort(sigAltNucs.begin(), sigAltNucs.end(), std::greater<>());

  orderedSigAltNucs.resize(sigAltNucs.size());
  for (size_t i = 0; i < sigAltNucs.size(); i++)
    this->orderedSigAltNucs[i] = sigAltNucs[i].getSigAltNuc();
}

int16_t SignificantAltNucs::getAltNucOrder(const u_int16_t & v) const
{
  for (int16_t i = 0; i < orderedSigAltNucs.size(); i++)
    if (v == orderedSigAltNucs[i])
      return i;

  return -1;
}

void SignificantAltNucs::resetSigAltNucs()
{
  sigAltNucs.clear();
  orderedSigAltNucs.clear();
}

AltNucReadCount::AltNucReadCount(AltNucReadCount &&v) noexcept:
    altNuc(v.altNuc),
    altNucChar(v.altNucChar),
    readCount(v.readCount),
    sigAltNucs(std::move(v.sigAltNucs))
{}

AltNucReadCount::AltNucReadCount(
    u_int16_t altNuc,
    char aluNucChar,
    u_int32_t readCount,
    SignificantAltNucs sigAltNucs
    ):
        altNuc(altNuc),
        altNucChar(aluNucChar),
        readCount(readCount),
        sigAltNucs(std::move(sigAltNucs))
{}

AltNucReadCount & AltNucReadCount::operator = (AltNucReadCount &&v) noexcept
{
  if (this != &v)
  {
    altNuc = v.altNuc;
    altNucChar = v.altNucChar;
    readCount = v.readCount;
    sigAltNucs = std::move(v.sigAltNucs);
  }

  return *this;
}

bool AltNucReadCount::operator > (const AltNucReadCount & v) const
{
  if (readCount != v.readCount)
    return readCount > v.readCount;

  const int16_t ord = sigAltNucs.getAltNucOrder(altNuc);
  const int16_t vOrd = v.sigAltNucs.getAltNucOrder(v.altNuc);

  if (ord != vOrd) {
    if (ord == -1 || vOrd == -1)
      return ord > vOrd;
    else
      return ord < vOrd;
  }

  return altNuc < v.altNuc;
}

u_int16_t AltNucReadCount::getAltNuc() const { return altNuc; }

char AltNucReadCount::getAltNucChar() const { return altNucChar; }

u_int32_t AltNucReadCount::getReadCount() const { return readCount; }

const SignificantAltNucs &AltNucReadCount::getSigAltNucs() const { return sigAltNucs; }

CellReadCounts::CellReadCounts(CellReadCounts &&v) noexcept:
    altNucReadCounts(std::move(v.altNucReadCounts)),
    cov(v.cov)
{}

CellReadCounts::CellReadCounts(
    const AltNucReadCount &v1,
    const AltNucReadCount &v2,
    const AltNucReadCount &v3,
    u_int32_t _cov
)
{
  altNucReadCounts[0] = v1;
  altNucReadCounts[1] = v2;
  altNucReadCounts[2] = v3;
  cov = _cov;
}

CellReadCounts::CellReadCounts(
    AltNucReadCount &&v1,
    AltNucReadCount &&v2,
    AltNucReadCount &&v3,
    u_int32_t _cov
) noexcept
{
  altNucReadCounts[0] = std::move(v1);
  altNucReadCounts[1] = std::move(v2);
  altNucReadCounts[2] = std::move(v3);
  cov = _cov;
}

CellReadCounts & CellReadCounts::operator = (CellReadCounts &&v) noexcept
{
  if (this != &v)
  {
    altNucReadCounts = std::move(v.altNucReadCounts);
    cov = v.cov;
  }

  return *this;
}

void CellReadCounts::sortAltNucReadCounts()
{
  std::sort(altNucReadCounts.begin(), altNucReadCounts.end(), std::greater<>());
}

std::ostream & operator << (std::ostream & out, const CellReadCounts & v)
{
  out << v.altNucReadCounts[0].getAltNucChar() << ","
      << v.altNucReadCounts[1].getAltNucChar() << ","
      << v.altNucReadCounts[2].getAltNucChar() << ";";

  out << v.altNucReadCounts[0].getReadCount() << ","
      << v.altNucReadCounts[1].getReadCount() << ","
      << v.altNucReadCounts[2].getReadCount() << ","
      << v.cov;

  return out;
}

const std::regex ChromosomeLabel::pattern(CHROMOSOME_PATTERN);

std::string ChromosomeLabel::getLabel(const std::string &v) {
  if (v.empty())
    throw std::runtime_error("Error! Empty chromosome label.");

  std::smatch match;

  if (std::regex_match(v, match, pattern))
    return match[LABEL_POSITION].str();
  else
    throw std::runtime_error("Error! Illegal format of chromosome label: " + v);
}

ChromosomeLabel::ChromosomeLabel(ChromosomeLabel &&v) noexcept:
    fullLabel(std::move(v.fullLabel)),
    label(std::move(v.label))
{}

ChromosomeLabel::ChromosomeLabel(const std::string &v):
    fullLabel(v),
    label(getLabel(v))
{}

ChromosomeLabel::ChromosomeLabel(std::string &&v) noexcept:
    fullLabel(std::move(v)),
    label(getLabel(fullLabel))
{}

ChromosomeLabel & ChromosomeLabel::operator = (ChromosomeLabel &&v) noexcept
{
  if (this != &v)
  {
    fullLabel = std::move(v.fullLabel);
    label = std::move(v.label);
  }

  return *this;
}

bool ChromosomeLabel::operator == (const ChromosomeLabel &v) const
{
  return label == v.label;
}

bool ChromosomeLabel::operator != (const ChromosomeLabel &v) const
{
  return label != v.label;
}

bool ChromosomeLabel::operator > (const ChromosomeLabel &v) const
{
  u_int16_t i, i1;

  try
  {
    i = std::stoi(label);
    i1 = std::stoi(v.label);
  }
  catch (std::invalid_argument &_v)
  {
    return label > v.label;
  }

  return i > i1;
}

bool ChromosomeLabel::operator < (const ChromosomeLabel &v) const
{
  u_int16_t i, i1;

  try
  {
    i = std::stoi(label);
    i1 = std::stoi(v.label);
  }
  catch (std::invalid_argument &_v)
  {
    return label < v.label;
  }

  return i < i1;
}

std::string ChromosomeLabel::getFullLabel() const { return fullLabel; }

std::ostream & operator << (std::ostream & out, const ChromosomeLabel &v)
{
  out << v.fullLabel;
  return out;
}

void ContinuousNoiseCounts::add(
    std::vector<u_int64_t> &vec,
    size_t val
    )
{
  if (val >= vec.size())
    vec.resize(val + 1, 0);

  vec[val]++;
}

void ContinuousNoiseCounts::alignSize(
    std::vector<u_int64_t> &v1,
    std::vector<u_int64_t> &v2
    )
{
  if (v1.size() < v2.size())
    v1.resize(v2.size(), 0);
  else if (v1.size() > v2.size())
    v2.resize(v1.size(), 0);
}

void ContinuousNoiseCounts::alignSize(ContinuousNoiseCounts &v)
{
  alignSize(m1, v.m1);
  alignSize(m2, v.m2);
  alignSize(m3, v.m3);
  alignSize(ref, v.ref);
  alignSize(cov, v.cov);
}

void ContinuousNoiseCounts::add(
    std::vector<u_int64_t> &v1,
    const std::vector<u_int64_t> &v2
    )
{
  std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::plus<>());
}

/**
 * Add each count to the corresponding vector.
 * @param counts first three members are read counts of alternative nucleotides, the fourth is coverage
 */
void ContinuousNoiseCounts::add(std::array<u_int32_t , 4> &&counts)
{
  std::sort(counts.begin(), counts.end() - 1, std::greater<>());

  add(m1, counts[0]);
  add(m2, counts[1]);
  add(m3, counts[2]);
  add(cov, counts[3]);
  add(ref, counts[3] - counts[0] - counts[1] - counts[2]);
}

ContinuousNoiseCounts &ContinuousNoiseCounts::operator+=(ContinuousNoiseCounts &v) {
  alignSize(v);

  add(m1, v.m1);
  add(m2, v.m2);
  add(m3, v.m3);
  add(ref, v.ref);
  add(cov, v.cov);

  return *this;
}

const std::vector<u_int64_t> &ContinuousNoiseCounts::getM1() const
{
  return m1;
}

const std::vector<u_int64_t> &ContinuousNoiseCounts::getM2() const
{
  return m2;
}

const std::vector<u_int64_t> &ContinuousNoiseCounts::getM3() const
{
  return m3;
}

const std::vector<u_int64_t> &ContinuousNoiseCounts::getRef() const
{
  return ref;
}

const std::vector<u_int64_t> &ContinuousNoiseCounts::getCov() const
{
  return cov;
}

std::ostream & NoiseCounts::printHelper(
    std::ostream &out,
    const std::vector<std::pair<size_t, u_int64_t>> &v
)
{
  out << v[0].first << "," << v[0].second;

  for (size_t i = 1; i < v.size(); i++)
    out << "\t" << v[i].first << "," << v[i].second;

  return out;
}

NoiseCounts::NoiseCounts(const ContinuousNoiseCounts &v)
{
  numPos = 0;

  collectCounts(m1, v.getM1(), nullptr);
  collectCounts(m2, v.getM2(), nullptr);
  collectCounts(m3, v.getM3(), nullptr);
  collectCounts(ref, v.getRef(), nullptr);
  collectCounts(cov, v.getCov(), &numPos);
}

NoiseCounts::NoiseCounts(NoiseCounts &&v) noexcept
{
  m1 = std::move(v.m1);
  m2 = std::move(v.m2);
  m3 = std::move(v.m3);
  ref = std::move(v.ref);
  cov = std::move(v.cov);

  numPos = v.numPos;
}

NoiseCounts & NoiseCounts::operator = (NoiseCounts &&v) noexcept
{
  if (this != &v) {
    m1 = std::move(v.m1);
    m2 = std::move(v.m2);
    m3 = std::move(v.m3);
    ref = std::move(v.ref);
    cov = std::move(v.cov);

    numPos = v.numPos;
  }

  return *this;
}

void NoiseCounts::collectCounts(
    std::vector<std::pair<size_t, u_int64_t>> &vec,
    const std::vector<u_int64_t> &val,
    u_int64_t * const numPos
    )
{
  for (size_t i = 0; i < val.size(); i++) {
    if (val[i] > 0) {
      vec.emplace_back(i, val[i]);

      if (numPos != nullptr)
        (*numPos) += val[i];
    }
  }
}

u_int64_t NoiseCounts::getBackgroundSitesNum() const { return numPos; }

std::ostream &operator<<(
    std::ostream &out,
    const NoiseCounts &v
    )
{
  if (v.numPos > 0) {
    NoiseCounts::printHelper(out, v.m1) << "\n";
    NoiseCounts::printHelper(out, v.m2) << "\n";
    NoiseCounts::printHelper(out, v.m3) << "\n";
    NoiseCounts::printHelper(out, v.ref) << "\n";
    NoiseCounts::printHelper(out, v.cov) << "\n";
  }

  return out;
}

bool compareDataEntry(const TDataEntry &v1, const TDataEntry &v2)
{
  const auto &c1 = std::get<0>(std::get<0>(v1));
  const auto &c2 = std::get<0>(std::get<0>(v2));

  if (c1 == c2)
  {
    const auto p1 = std::get<1>(std::get<0>(v1));
    const auto p2 = std::get<1>(std::get<0>(v2));

    return p1 < p2;
  }

  return c1 < c2;
}

void Data::addCellName(const std::string &v) { cellNames.emplace_back(v); }

void Data::addCellName(std::string &&v) { cellNames.emplace_back(std::move(v)); }

void Data::addReadCounts(
    const std::vector<std::pair<std::tuple<ChromosomeLabel, unsigned int, char, std::vector<char>>,std::vector<CellReadCounts>>> &v
    )
{
  std::copy(v.begin(), v.end(), std::back_inserter(readCounts));
}

void Data::addReadCounts(
    std::vector<std::pair<std::tuple<ChromosomeLabel, unsigned int, char, std::vector<char>>,std::vector<CellReadCounts>>> &&v
    )
{
  std::move(v.begin(), v.end(), std::back_inserter(readCounts));
}

void Data::addContinuousNoiseCounts(ContinuousNoiseCounts &v)
{
  continuousNoiseCounts += v;
}

void Data::createNoiseCounts() {
  noiseCounts = NoiseCounts(continuousNoiseCounts);
}

void Data::sortEntries()
{
  std::sort(
      readCounts.begin(),
      readCounts.end(),
      compareDataEntry
      );
}

const std::vector<std::string> &Data::getCellNames() const { return cellNames; }

const TData &Data::getCandidateMutatedReadCounts() const { return readCounts; }

u_int32_t Data::getCandidateMutatedSitesNum() const { return readCounts.size(); }

u_int64_t Data::getBackgroundSitesNum() const { return noiseCounts.getBackgroundSitesNum(); }

const NoiseCounts &Data::getBackgroundSitesReadCounts() const { return noiseCounts; }
