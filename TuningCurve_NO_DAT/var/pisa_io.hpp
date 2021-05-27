#ifndef VARIATOR_PISA_IO_HPP
#define VARIATOR_PISA_IO_HPP

#include "pop.hpp"

#include <fstream>

class pisa_io {
public:
  pisa_io() {}
  pisa_io(const char* prefix)
    : prefix(prefix) {}
  void init(population& pop);

  void write_state(int state);
  int read_state();
  bool check_file(const char* file);
  void write_pop(const char* file, population& pop);
  int* read_sel();
  void clear_arc();
  //void append_his();

private:
  std::string filename(const char* file);

  const char* prefix;
  int *sel_ids_buf;
  int buf_size;
};

#endif // VARIATOR_PISA_IO_HPP
