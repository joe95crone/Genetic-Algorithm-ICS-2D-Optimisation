#include "config_parser.hpp"
#include "var_app.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>

std::ofstream diag;

void usage(char* progname)
{
  std::cerr << "Usage: " << progname
	    << "conf_file pisa_prefix poll_interval"  << std::endl
	    << "E.g.: " << progname << "var.conf ../pisa/ 0.2" << std::endl;
}

int main(int argc, char* argv[])
{
  if (argc != 4) {
    usage(argv[0]);
    return 1;
  }

  char *conf_file = argv[1];
  char *pisa_prefix = argv[2];

  double poll_interval;
  std::istringstream strm(argv[3]);
  strm >> poll_interval;
  if (strm.bad()) {
    std::cerr << "Invalid polling interval" << std::endl;
    exit(1);
  }

  diag.open("var_diag.log", std::ios::out);

  var_app app(conf_file, pisa_prefix, poll_interval);
  app.run();

  return 0;
}
