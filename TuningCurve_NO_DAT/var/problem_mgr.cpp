#include "problem_mgr.hpp"
#include "pop.hpp"

#include <unistd.h>		// getcwd
#include <cstdio>		// popen etc.
#include <cstdlib>		// exit
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

static std::string getabs(const char* file)
{
  char cwd[256];
  getcwd(cwd, 256);
  std::string abs = cwd;
  abs.append("/");
  return abs.append(file);
}


problem_mgr::problem_mgr(const char* problem_dir,
			 const char* launch_script,
			 const char* extract_script,
			 const char* done_flag_file)
  : problem_dir(problem_dir), done_flag_file(done_flag_file)
{
  this->launch_script = getabs(launch_script);
  this->extract_script = getabs(extract_script);

  // cleanup old crap
  std::ostringstream cmd;
  cmd << "rm -rf " << problem_dir << "/run.*";
  system(cmd.str().c_str());
}


static void rd_err()
{
  std::cerr << "Invalid variables specification" << std::endl;
  exit(1);
}


// XXX get error handling from old version
void problem_mgr::init_decvars(const char* variables_spec,
			       population& pop)
{
  std::ifstream is(variables_spec);
  if (!is.good()) {
    std::cerr << "Error opening variables specfication" << std::endl;
    exit(1);
  }

  pop.num_dec = 0;

  is >> std::ws;
  while (is.good()) {
    std::string line;
    std::getline(is, line);
    std::istringstream lnstrm(line);

    pop.ini_regions.push_back(population::region());
    population::region& reg = pop.ini_regions.back();

    int here_num_dec = 0;	// how many on this line

    while (lnstrm.good()) {
      population::interval iv;
      lnstrm >> iv.lbound;
      lnstrm >> iv.ubound;
      reg.push_back(iv);
      ++here_num_dec;
    }

    if (pop.num_dec == 0)
      pop.num_dec = here_num_dec;
    else if (here_num_dec != pop.num_dec) {
      std::cerr << "Inconsistent number of variables specified"
		<< std::endl;
      exit(1);
    }

    is >> std::ws;
  }
}



std::string problem_mgr::get_rundir(int id)
{
  std::string ret;
  std::ostringstream strm;

  strm << problem_dir
       << "/run." << std::setw(5)
       << std::setfill('0') << id;

  return strm.str();
}


void problem_mgr::eval(population& pop)
{
  // launch jobs
  for (int id = pop.this_gen_begin(); id < pop.this_gen_end(); ++id) {
    individual& indiv = pop.get(id);

    std::string rundir = get_rundir(id);

    std::string mkdir = "mkdir ";
    mkdir.append(rundir);
    system(mkdir.c_str());

    std::ostringstream launch_strm;
    launch_strm << "cd " << rundir << "; " << launch_script << " ";

    for (int v = 0; v < pop.num_dec; ++v)
      launch_strm << indiv.dec[v] << " ";

    launch_strm << "&";		// so we don't block
    system(launch_strm.str().c_str());
  }

  // harvest jobs
  bool done;
  do {
    sleep(10);
    done = reap(pop);
  } while (!done);
}


static void extract_err()
{
  std::cerr << "Error calling extraction script" << std::endl;
  exit(1);
}


static void read_doubles(FILE* pipe, double* arr, int num)
{
  char buf[256];
  if (fgets(buf, 256, pipe) == NULL)
    extract_err();
  std::istringstream strm(buf);
  for (int i = 0; i < num; ++i) {
    strm >> arr[i];
    if (strm.bad()) extract_err();
  }
}


bool problem_mgr::done_file_exists(std::string rundir)
{
  std::string& donefile = rundir;
  donefile.append("/").append(done_flag_file);
  std::ifstream ifs(donefile.c_str());
  return ifs.good();
}


bool problem_mgr::reap(population& pop)
{
  bool alldone = true;
  
  for (int id = pop.this_gen_begin(); id < pop.this_gen_end(); ++id) {
    individual& indiv = pop.get(id);
    if (indiv.job_finished) continue;

    alldone = false;

    std::string rundir = get_rundir(id);
    if (!done_file_exists(rundir)) continue;
    indiv.job_finished = true;

    std::ostringstream extract_cmd;
    extract_cmd << "cd " << rundir << "; " << extract_script;
    // pass decision vars in case needed for constraint calc
    for (int v = 0; v < pop.num_dec; ++v)
      extract_cmd << " " << indiv.dec[v];

    FILE* pipe = popen(extract_cmd.str().c_str(), "r");
    if (pipe == NULL) extract_err();

    // read objectives
    read_doubles(pipe, indiv.obj, pop.num_obj);

    // read constraints
    if (pop.num_con != 0)
      read_doubles(pipe, indiv.con, pop.num_con);

    pclose(pipe);
  }

  return alldone;
}
