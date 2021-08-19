// var_app.hpp
#ifndef VARIATOR_VAR_APP_HPP
#define VARIATOR_VAR_APP_HPP

#include "config_parser.hpp"
#include "misc.hpp"
#include "pop.hpp"
#include "pisa_io.hpp"
#include "problem_mgr.hpp"

#include <fstream>

class var_app {
public:
  var_app(const char* conf_file, const char* pisa_prefix,
	  double poll_interval);
  void run();

private:
  void read_var_conf(const char* conf_file);
  void init_output(int rand_seed);
  void config_pop(const char* pisa_prefix);
  bool finished() { return pop.cur_gen == pop.max_gen; }
  void set_state(int st);
  void zzz();
  void state0();
  void state2();
  void dump_gen();

  int state;
  int poll_sec, poll_usec;
  config_parser var_conf;
  population pop;
  problem_mgr problem;
  pisa_io io;
  std::ofstream dump_stream;
};

#endif // VARIATOR_VAR_APP_HPP
