//
// Created by senbaikang on 19.05.21.
//

#include <algorithm>

#include <datafilter/read_counts_file.h>

ReadCountsFile::ReadCountsFile(std::string v):
    fileName(std::move(v)),
    status(false),
    clearCache(false)
{}

ReadCountsFile::operator bool() const { return status; }

bool ReadCountsFile::operator!() const { return status; }

const std::string &ReadCountsFile::getFileName() const { return fileName; }

void ReadCountsFile::setClearCache(bool val) { clearCache = val; }

MPileupFile::MPileupFile(std::string v):
    ReadCountsFile(std::move(v)),
    iStream(fileName, std::ifstream::in)
{
  clearCache = false;
  if (!iStream)
    throw std::runtime_error("Error: Could not open input file " + fileName);
}

MPileupFile::~MPileupFile() { iStream.close(); }

void MPileupFile::rewind()
{
  iStream.clear();
  iStream.seekg(0);
  status = true;
}

bool MPileupFile::getLine(std::string &line)
{
  getline(iStream, line);

  updateStatus();
  return status;
}

void MPileupFile::updateStatus() { status = static_cast<bool>(iStream); }

bool MPileupFile::isGood() const { return iStream.good(); }

bool MPileupFile::isEOF() const { return iStream.eof(); }

GZFile::GZFile(std::string v, u_int64_t l):
    ReadCountsFile(std::move(v)),
    errorCode(Z_OK),
    eof(false)
{
  clearCache = false;
  readLength = l;
  byteContent = new char[l]();

  file = (gzFile)gzopen(fileName.c_str(), "rb");

  if (file == nullptr)
    throw std::runtime_error("Error: Could not open input file " + fileName);

  gzrewind(file);
}

GZFile::~GZFile()
{
  gzclose(file);
  delete[] byteContent;
}

void GZFile::rewind()
{
  cache.clear();
  cache.seekg(0);
  startPoint = 0;

  errorCode = Z_OK;
  eof = false;
  status = true;
}

bool GZFile::getLine(std::string &line)
{
  line.clear();
  const u_int64_t tmpStartPoint = cache.str().find_first_of('\n', startPoint);
  if (tmpStartPoint != std::string::npos)
  {
    startPoint = tmpStartPoint + 1;
    getline(cache, line);

    return true;
  }
  else
  {
    if (clearCache)
    {
      getline(cache, line);
      cache.clear();
      cache.seekg(0);
      cache = std::stringstream();
      cache << line;
      line.clear();
      startPoint = 0;
    }

    const int len = gzread(file, byteContent, readLength);
    cache << std::string(byteContent, len);

    if (len != readLength && !eof) {
      eof = true;

      if (byteContent[len == 0 ? 0 : len - 1] != '\n')
        cache << '\n';
    }

    updateStatus();
    return status;
  }
}

void GZFile::updateStatus()
{
  const char *errMsg = gzerror(file, &errorCode);

  if (errorCode < 0)
    throw std::runtime_error("Error: Could not read " + fileName + "properly. " + errMsg);

  status = (isGood() && !isEOF()) || !cache.str().empty();
}

bool GZFile::isGood() const { return errorCode == Z_OK; }

bool GZFile::isEOF() const { return eof; }
