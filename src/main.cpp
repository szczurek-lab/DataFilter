#include "datafilter.h"

int main(int argc, char* argv[])
{
  Config config;
  parseCommandLineArgs(config, argc, argv);
  process_data(config);

  return 0;
}
