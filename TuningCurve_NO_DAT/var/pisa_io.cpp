#include "pisa_io.hpp"
#include "pop.hpp"
#include "misc.hpp"

#include <cstdlib>		// exit
#include <string>
#include <fstream>
#include <iostream>

#include <sstream>

void pisa_io::init(population& pop)
{
  buf_size = pop.gen_size;
  sel_ids_buf = new int[buf_size];
}


std::string pisa_io::filename(const char* file)
{
  std::string ret = prefix;
  ret.append(file);
  return ret;
}


void pisa_io::write_state(int state)
{
  std::ofstream ofs(filename("sta").c_str(),
		    std::ios::out | std::ios::trunc);
  ofs << state;
}


int pisa_io::read_state()
{
  std::ifstream ifs(filename("sta").c_str());
  int state;
  ifs >> state;
  return state;
}


// check to see if file contains only '0'
// this method is slated for deletion
bool pisa_io::check_file(const char* file)
{
  std::ifstream ifs(filename(file).c_str());
  if (!ifs.good()) return false;
  int val;
  ifs >> val;
  return val == 0;
}


// file is either "ini" or "var"
void pisa_io::write_pop(const char* file, population& pop)
{
  std::ofstream ofs(filename(file).c_str(),
		    std::ios::out | std::ios::trunc);
  ofs << pop.gen_size * (pop.num_obj + pop.num_con + 1) << std::endl;

  for (int id = pop.this_gen_begin(); id < pop.this_gen_end(); ++id) {
    individual& indiv = pop.get(id);
    ofs << id << " ";

    for (int i = 0; i < pop.num_obj; ++i)
      ofs << indiv.obj[i] << " ";

    for (int i = 0; i < pop.num_con; ++i)
      ofs << indiv.con[i] << " ";

    ofs << std::endl;
  }

  ofs << "END";
}


int* pisa_io::read_sel()
{
  std::ifstream ifs(filename("sel").c_str());

  int size;
  ifs >> size;
  if (size != buf_size) {
    std::cerr << "selector isn't giving us enough" << std::endl;
    exit(1);
  }

  for (int i = 0; i < buf_size; ++i)
    ifs >> sel_ids_buf[i];

  std::string tag;
  ifs >> tag;
  if (tag.compare("END")) {
    std::cerr << "couldn't find END tag" << std::endl;
    exit(1);
  }

  ifs.close();

  // save selection history for diagnostics
  std::ostringstream cmd;
  cmd << "cat " << filename("sel") << " >> "
      << filename("sel.his") << "; echo \"\n\" >> "
      << filename("sel.his");	// shell scripting in C++!
  system(cmd.str().c_str());

  // clear sel file for selector
  std::ofstream ofs(filename("sel").c_str(),
		    std::ios::out | std::ios::trunc);
  ofs << 0;

  return sel_ids_buf;
}


void pisa_io::clear_arc()
{
  std::ofstream ofs(filename("arc").c_str(),
		    std::ios::out | std::ios::trunc);
  ofs << 0;
}
