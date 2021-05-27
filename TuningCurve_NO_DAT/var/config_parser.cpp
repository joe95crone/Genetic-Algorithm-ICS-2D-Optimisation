#include "config_parser.hpp"

#include <iostream>
#include <cstdlib> // exit


#define CONFIG_PARSER_READER(type, field)	\
  template <>					\
  void config_parser::opt_reader<type>::read	\
  (std::istream& is, optdata& data)		\
  {						\
    is >> data.field;				\
  }

CONFIG_PARSER_READER(int, i)
CONFIG_PARSER_READER(double, d)


template <>
void config_parser::opt_reader<bool>::read
(std::istream& is, optdata& data)
{
  std::string s;
  is >> s;
  if (s.compare("true") == 0)
    data.b = true;
  else if (s.compare("false") == 0)
    data.b = false;
  else {
    std::cerr << "Error: bad boolean value" << std::endl;
    exit(1);
  }
}


template <>
void config_parser::opt_reader<const char*>::read
(std::istream& is, optdata& data)
{
  // XXX use util string2cstr()
  std::string s;
  getline(is, s);
  int size = s.size();
  char *cs = new char[size+1];	// XXX leak: use auto_ptr instead
  data.s = strncpy(cs, s.data(), size);
  cs[size] = '\0';
}


// ----------------------------------------

using namespace std;
void config_parser::parse(std::istream& is)
{
  while (is.good()) {
    std::string str;
    is >> str;

    if (str.empty()) continue;

    // XXX throw an exception
    if (dict.count(str) == 0) {
      std::cerr << "Error: unknown option \"" << str
		<< "\"" << std::endl;
      exit(1);
    }

    is >> std::ws; // skip whitespace

    // XXX throw an exception
    if (is.eof()) {
      std::cerr << "Error: empty option \"" << str
		<< "\"" << std::endl;
    }

    optentry &opt = dict[str];
    opt.reader->read(is, opt.data);
    opt.set = true;

    // XXX throw an exception
    if (is.fail()) {
      std::cerr << "Error reading option \"" << str
		<< "\"" << std::endl;
      exit(1);
    }
  }
}
