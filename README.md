# DataFilter

DataFilter identifies candidate somatic variant sites from single-cell DNA sequencing data. The selection process is similar to that in [SCIPhI](https://www.nature.com/articles/s41467-018-07627-7) with the following modifications:

1. The probability of observing data <img src="https://render.githubusercontent.com/render/math?math=D_{ij}"> for the wildtype genotype (0/0) is <img src="https://render.githubusercontent.com/render/math?math=P_{wt}(D_{ij}) = P(s_{ij} | c_{ij}, \frac{1}{3} f_{wt}, w_{wt})">.
2. The probability of observing data <img src="https://render.githubusercontent.com/render/math?math=D_{ij}"> for the single mutant genotype (0/1) is <img src="https://render.githubusercontent.com/render/math?math=P_{a}(D_{ij}) = P(s_{ij} | c_{ij}, \frac{1}{2} - \frac{1}{3} f_{wt}, w_{a})">.
2. It is able to capture sites containing double mutant genotypes (1/1 or 1/1').

## Minimal requirements

1. cmake 3.10
2. C++ 17
3. zlib 1.2.12
4. Boost 1.79
5. dlib 19.24

## Installation

DataFilter uses `conan` to manage dependencies and `cmake` to build the executable file. Please make sure that `conan` and `cmake` are available in your `PATH`.

1. Clone this repository to local:

    ```shell
    $ git clone https://github.com/szczurek-lab/DataFilter.git
    $ cd DataFilter
    $ mkdir build; cd build
    ```

2. Use `conan` to install dependencies:

   ```shell
   $ conan install .. -s build_type=Release
   ```

   This normally goes well under macOS with Apple clang and Linux with gcc. But in case you prefer clang on Linux, you might have to follow the command below instead:

   ```shell
   $ conan install .. -s build_type=Release --profile CONAN_PROFILE_FOR_CLANG --build=missing
   ```

   `CONAN_PROFILE_FOR_CLANG` is the `conan` profile for clang. Please refer to the [official documentation of conan](https://docs.conan.io/en/latest/reference/profiles.html) for more information.

3. Use `cmake` to build the executable:

    ```shell
    $ cmake -D CMAKE_BUILD_TYPE=Release ..; cmake --build .
    ```

    In case of using an undefault compiler, such as clang, on Linux, please use:

    ```shell
    $ cmake -D CMAKE_C_COMPILER=CLANG_BIN -D CMAKE_CXX_COMPILER=CLANG++_BIN -D CMAKE_BUILD_TYPE=Release ..; cmake --build .
    ```

    Where `CLANG_BIN` and `CLANG++_BIN` are the clang executables on your platform.

4. `datafilter` executable is located under `build/bin/`.

### Tested platforms

1. Ubuntu 20.04, gcc 10.3.0, 11.1.0.
2. Ubuntu 20.04, clang 14.0.6.
3. macOS 12.6, Apple clang 14.0.0.

## Usage

```shell
$ ./bin/datafilter -h

Generic options:
  -h [ --help ]          Print help message.

Required options:
  --cellNames arg        File name of names of the BAM files used to create the
                         mpileup.
  -o [ --out ] arg       Output file name.

Required but mutually exclusive options:
  --inFile arg           File name of new line separated paths to input mpileup
                         files.
  --in arg               Multiple input mpileup files.

Optional options:
  -t [ --threads ] arg   Number of threads. [1]
  --cellNameSuf arg      Suffix in regex form of cell names to remove. 
                         [\\..*?$]
  --ex arg               File name of exclusion list (VCF format), containing 
                         loci which should be ignored.
  --me arg               File name of mutations to exclude during the 
                         sequencing error rate estimation (VCF format).
  --inc arg              File name of inclusion list (VCF format) containing 
                         Variants (CHROM, POS) that should be included.
  --mi arg               File name of inclusion list (VCF format) containing 
                         Variants (CHROM, POS, REF, ALT) that should be 
                         included.
  --readLength arg       Read length of gunzipped input mpileup. [524288]
  --mupr arg             Prior mutation rate. [0.0001]
  --germpr arg           Prior germline rate. [0.001]
  --adopr arg            Prior allelic drop-out rate. [0.1]
  --wildMean arg         Mean error rate. If the effective sequencing error 
                         rate should not be learned "--ese 0" one can specify 
                         it. [0.001]
  --wildOverDis arg      Initial overdispersion for wild type. [100.0]
  --muOverDis arg        Initial overdispersion for mutant type. [2.0]
  --ese arg              Estimate the effective sequencing error rate. [on]
  --maxese arg           Max number of sites per input file for estimating 
                         effective sequencing error rate. [100000]
  --cwm arg              Number of tumor cells required to have a mutation in 
                         order to be called. [2]
  --mnp arg              Number of cells which need to pass the filters 
                         described below. [2]
  --mcn arg              Minimum coverage required per normal cell. [5]
  --mc arg               Minimum coverage required per cell. [1]
  --ms arg               Minimum number of reads required to support the 
                         alternative. [3]
  --mf arg               Minimum required frequency of reads supporting the 
                         alternative per cell. [0.0]
  --mff arg              Mean of acceptable variant allele frequency across all
                         cells for a specific locus. Mapping artifacts may 
                         result in low allele frequencies across cells. In 
                         order to filter these events out we apply a 
                         log-likelihood ratio test where the sequencing error 
                         model has a mean of this value. [0.25]
  --bns arg              Loci with up to this number of alternative supporting 
                         reads in the bulk control sample will be skipped as 
                         germline. [2]
  --bnc arg              Minimum required coverage of reads in the bulk control
                         sample. [6]
  --ncf arg              Normal cell filter. Currently there are three options:
                         (0) Do not use the normal cells for filtering; (1) use
                         a simple filtering scheme excluding mutations if the 
                         probability of being mutated is higher than not being 
                         mutated for any cell independently; (2) filter 
                         mutations where the probability that at least one cell
                         is mutated is higher than no cell is mutated. Note 
                         that in contrast to (1) the cells are not independent 
                         and cells with no alternative support need to be 
                         explained via dropout events. [1]
  --mnc arg              Maximum number of control cells allowed to be mutated.
                         [0]
```



### Input data

DataFilter only accepts `mpileup` files. Since they usually occupy a lot of disk space, DataFilter also supports reading gunzipped files (`*.pileup.gz`) directly. In this case, different buffer sizes (specified by `--readLength` flag) can result in distinct running time and memory usage. If it is possible, please use `mpileup` files as input.

To accelerate data processing, the programme includes multithreading, where each thread works on one file at a time. Therefore, the best practice is that after collecting read counts data from all single-cells into a single `mpileup` file, one should extract data by chromosomes into separate files. To run `datafilter`, one can use either `--in` flag, followed by all the independent data files, or `--inFile` flag, followed by a single file containing the paths to data files. 

### Input cell names

The `--cellNames` flag specifying a file containing cell names, which is compatible with [SCIPhI](https://github.com/cbg-ethz/SCIPhI).

### Output file

The output file is specified by `-o` or `--out` flag, which is readily available as the input to [SIEVE](https://github.com/szczurek-lab/SIEVE)'s DataCollector.


## Funding

This project received support from
- the Polish National Science Centre SONATA BIS grant No. 2020/38/E/NZ2/00305,
- European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No. 766030, as well as
- European Research Council (ERC-617457-PHYLOCANCER), the Spanish Ministry of Science and Innovation (PID2019-106247GB-I00), and Xunta de Galicia.

