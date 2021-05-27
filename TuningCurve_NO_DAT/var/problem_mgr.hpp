#ifndef VARIATOR_PROBLEM_MGR_HPP
#define VARIATOR_PROBLEM_MGR_HPP

#include "pop.hpp"

#include <string>

class problem_mgr {
public:
  problem_mgr() {};
  problem_mgr(const char* problem_dir, const char* launch_script,
	      const char* extract_script,
	      const char* done_flag_file);
  void init_decvars(const char* variables_spec, population& pop);
  void eval(population& pop);
  bool reap(population& pop);	// returns true when all jobs are done

private:
  std::string get_rundir(int id);
  bool done_file_exists(std::string rundir);

  const char* problem_dir;
  const char* done_flag_file;
  std::string launch_script;
  std::string extract_script;
};

#endif // VARIATOR_PROBLEM_MGR_HPP
