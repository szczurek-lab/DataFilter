//
// Created by senbaikang on 19.05.21.
//

#ifndef DATAFILTER_READ_COUNTS_FILE_H
#define DATAFILTER_READ_COUNTS_FILE_H

#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>

class ReadCountsFile {

protected:
    const std::string fileName;

    bool status;

    bool clearCache;

public:
    ReadCountsFile() = delete;

    explicit ReadCountsFile(std::string);

    ReadCountsFile(const ReadCountsFile &) = delete;

    ReadCountsFile(ReadCountsFile &&) = delete;

    virtual ~ReadCountsFile() = default;

    void setClearCache(bool val);

    virtual void rewind() = 0;

    virtual bool getLine(std::string &) = 0;

    virtual void updateStatus() = 0;

    [[nodiscard]] virtual bool isGood() const = 0;

    [[nodiscard]] virtual bool isEOF() const = 0;

    explicit operator bool() const;

    bool operator!() const;

    [[nodiscard]] const std::string &getFileName() const;

};


class MPileupFile : public ReadCountsFile {

protected:
    std::ifstream iStream;

public:
    MPileupFile() = delete;

    explicit MPileupFile(std::string);

    ~MPileupFile() override;

    void rewind() override;

    bool getLine(std::string &) override;

    void updateStatus() override;

    [[nodiscard]] bool isGood() const override;

    [[nodiscard]] bool isEOF() const override;

};


class GZFile : public ReadCountsFile {

protected:
    gzFile file;
    int errorCode; // To receive the error code from ZLIB.
    bool eof; // Whether reached the end of file.

    // Byte content read from the compressed input file.
    char *byteContent;
    u_int64_t readLength;

    // Start point to search for the new line.
    u_int64_t startPoint = 0;

    std::stringstream cache;

public:
    GZFile() = delete;

    explicit GZFile(std::string, u_int64_t);

    ~GZFile() override;

    void rewind() override;

    bool getLine(std::string &) override;

    void updateStatus() override;

    [[nodiscard]] bool isGood() const override;

    [[nodiscard]] bool isEOF() const override;

};

#endif // DATAFILTER_READ_COUNTS_FILE_H
