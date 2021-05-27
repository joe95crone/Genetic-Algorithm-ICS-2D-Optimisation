#ifndef VAR_CONFIG_PARSER_HPP
#define VAR_CONFIG_PARSER_HPP

#include <cstdlib>
#include <cstring>
#include <map>
#include <iostream>
#include <fstream>


class config_parser {
public:
  // add option w/o default
  template <typename T>
  config_parser& add(const char* name)
  {
    optentry &entry = dict[name];
    entry.reader = new opt_reader<T>;
    entry.set = false;
    return *this;
  }

  // add option w/ default
  template <typename T>
  config_parser& add(const char* name, T def);

  // parse config
  void parse(std::istream& is);
  void parse(const char* filename)
  {
    std::ifstream ifs(filename);
    parse(ifs);
  }

  // get option
  template <typename T>
  T get(const char* name);

  //private:
  union optdata {
    int i;
    bool b;
    double d;
    const char* s;
  };

  // we store subclasses (opt_reader<...>) of this to parse the options
  class opt_reader_base {
  public:
    virtual void read(std::istream& is, optdata& data) {};
  };

  template <typename T>
  class opt_reader : public opt_reader_base {
  public:
    void read(std::istream& is, optdata& data);
  };

  // one for each option
  struct optentry {
    opt_reader_base *reader;
    bool set;
    optdata data;
  };

  std::map<std::string, optentry> dict;
};

#define CONFIG_PARSER_DEF_ADDER(type, field)	\
  template <> inline				\
  config_parser&				\
  config_parser::add<type>			\
  (const char* name, type def)			\
  {						\
    add<type>(name);				\
    optentry &entry = dict[name];		\
    entry.data.field = def;			\
    entry.set = true;				\
    return *this;				\
  }


#define CONFIG_PARSER_GETTER(type, field)		\
  template <> inline					\
  type config_parser::get<type>				\
  (const char* name)					\
  {							\
    if (dict.count(name) == 0) {			\
      std::cerr << "Error: empty config parameter "	\
		<< name << std::endl;			\
      exit(1);						\
    }							\
    return dict[name].data.field;			\
  }

#define CONFIG_PARSER_DEFINE_FIELD(type, field)	\
  CONFIG_PARSER_DEF_ADDER(type, field)		\
  CONFIG_PARSER_GETTER(type, field)


CONFIG_PARSER_DEFINE_FIELD(int, i)
CONFIG_PARSER_DEFINE_FIELD(bool, b)
CONFIG_PARSER_DEFINE_FIELD(double, d)
CONFIG_PARSER_DEFINE_FIELD(const char*, s)

#endif // VAR_CONFIG_PARSER_HPP
