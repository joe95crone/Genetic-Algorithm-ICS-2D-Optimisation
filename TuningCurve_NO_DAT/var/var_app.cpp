#include "var_app.hpp"
#include "misc.hpp"
#include "config_parser.hpp"
#include "pisa_io.hpp"
#include "pop.hpp"
#include "problem_mgr.hpp"

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include <fstream>

#include <unistd.h>

var_app::var_app(const char* conf_file, const char* pisa_prefix,
		 double poll_interval)
{
  // read global configuration variables
  read_var_conf(conf_file);	// sets var_prop

  // set up job launcher
  problem = problem_mgr
    (var_conf.get<ccstr>("problem_dir"),
     var_conf.get<ccstr>("launch_script"),
     var_conf.get<ccstr>("extract_script"),
     var_conf.get<ccstr>("done_flag_file"));

  config_pop(pisa_prefix);   // set up population (size, etc.)
  problem.init_decvars	     // set bounds etc. for decision variables
    (var_conf.get<ccstr>("variables_spec"),
     pop);
  pop.init();

  // set up PISA io manager
  io = pisa_io(pisa_prefix);
  io.init(pop);

  // PISA polling interval
  poll_sec = int(poll_interval);
  poll_usec = int((poll_interval-poll_sec)*1e6);
}


void var_app::config_pop(const char* pisa_prefix)
{
  pop.max_gen = var_conf.get<int>("generations");
  
  config_parser pisa_conf;
  pisa_conf.add<int>("initial_population_size");
  pisa_conf.add<int>("parent_set_size");
  pisa_conf.add<int>("offspring_set_size");
  pisa_conf.add<int>("objectives");
  pisa_conf.add<int>("constraints");

  std::string cfg_filename = pisa_prefix;
  cfg_filename.append("cfg");
  pisa_conf.parse(cfg_filename.c_str());

  pop.gen_size = pisa_conf.get<int>("initial_population_size");
  if (pop.gen_size != pisa_conf.get<int>("parent_set_size") ||
      pop.gen_size != pisa_conf.get<int>("offspring_set_size")) {
    std::cerr << "initial_population_size must be equal to "
      "parent_set_size and offspring_set_size :(" << std::endl;
    exit(1);
  }

  pop.num_obj = pisa_conf.get<int>("objectives");
  pop.num_con = pisa_conf.get<int>("constraints");
  // num_dec is set by problem_mgr::init_decvars
}


void var_app::init_output(int rand_seed)
{
  dump_stream.open("var_output.txt", std::ios::out | std::ios::app);
  time_t tm = time(0);
  std::string tm_str = asctime(localtime(&tm));
  tm_str.erase(tm_str.size()-1);	// remove newline
  dump_stream << "# Beginning of run: " << tm_str
	      << " (seed: " << rand_seed << ")"
	      << std::endl << "#" << std::endl;
}


void var_app::read_var_conf(const char* conf_file)
{
  // XXX add init data option
  var_conf.add<ccstr>("problem_dir");
  var_conf.add<ccstr>("variables_spec");
  var_conf.add<ccstr>("launch_script");
  var_conf.add<ccstr>("extract_script");
  var_conf.add<ccstr>("done_flag_file");
  var_conf.add<int>("seed");
  var_conf.add<int>("generations");
  var_conf.add<double>("individual_mutation_probability");
  var_conf.add<double>("variable_mutation_probability");
  var_conf.add<double>("variable_swap_probability");
  var_conf.add<double>("individual_recombination_probability");
  var_conf.add<double>("variable_recombination_probability");
  var_conf.add<double>("eta_mutation");
  var_conf.add<double>("eta_recombination");
  var_conf.add<double>("shrink_threshold");
  var_conf.add<double>("shrink_margin");
  var_conf.add<double>("eta_mutation_vary_threshold", 0);
  var_conf.add<double>("eta_mutation_vary_factor", 1);
  var_conf.add<double>("eta_recombination_vary_threshold", 0);
  var_conf.add<double>("eta_recombination_vary_factor", 1);

  var_conf.parse(conf_file);

  int seed = var_conf.get<int>("seed");
  if (seed == 0) seed = time(0);
  srand(seed);
  init_output(seed);

  pop.variation_prop.indiv_mut_prob = var_conf.get<double>
    ("individual_mutation_probability");
  pop.variation_prop.var_mut_prob = var_conf.get<double>
    ("variable_mutation_probability");
  pop.variation_prop.var_swap_prob = var_conf.get<double>
    ("variable_swap_probability");
  pop.variation_prop.indiv_recomb_prob = var_conf.get<double>
    ("individual_recombination_probability");
  pop.variation_prop.var_recomb_prob = var_conf.get<double>
    ("variable_recombination_probability");
  pop.variation_prop.eta_mut = var_conf.get<double>
    ("eta_mutation");
  pop.variation_prop.eta_recomb = var_conf.get<double>
    ("eta_recombination");
  pop.variation_prop.shrink_threshold = var_conf.get<double>
    ("shrink_threshold");
  pop.variation_prop.shrink_margin = var_conf.get<double>
    ("shrink_margin");
  pop.variation_prop.eta_mut_vary_threshold = var_conf.get<double>
    ("eta_mutation_vary_threshold");
  pop.variation_prop.eta_mut_vary_factor = var_conf.get<double>
    ("eta_mutation_vary_factor");
  pop.variation_prop.eta_recomb_vary_threshold = var_conf.get<double>
    ("eta_recombination_vary_threshold");
  pop.variation_prop.eta_recomb_vary_factor = var_conf.get<double>
    ("eta_recombination_vary_factor");
}


void var_app::set_state(int st)
{
  state = st;
  io.write_state(st);
}


void var_app::zzz()
{
  sleep(poll_sec);
  usleep(poll_usec);
}


void var_app::run()
{
  set_state(0);		    // tell selector we're in state 0
  state0();		    // generate and evaluate first generation
  dump_gen();   	    // write this gen to output file
  io.write_pop("ini", pop); // give to selector

  set_state(1);		    // tell selector to select initial parents
  do zzz();		    // wait until it's done
  while ((state = io.read_state()) == 1);

  while (true) {		   // main loop
    if (state == 7 || state == 11) // selector has terminated/reset
      return;

    if (state == 8) {		   // selector wants us to restart
      std::cerr << "Error: restart not implemented" << std::endl;
      exit(1);
    }

    if (state == 2) {	        // selector has chosen parents
      state2();		        // read from selector, variate, evaluate
      dump_gen();

      if (finished()) {		// all generations complete?
	set_state(5);		// tell selector to wind down
	dump_stream << "#" << std::endl;
	return;
      }

      io.write_pop("var", pop);	// give offspring to selector
      io.clear_arc();		// finnicky selector demands this
      set_state(3);		// tell selector to select
    }

    zzz();
    state = io.read_state();
  }
}


// create, evaluate, and write initial generation
void var_app::state0()
{
  diag << "Generation 1" << std::endl;

  // XXX allow pre-initialized
  pop.fill_first_gen();

  problem.eval(pop);
}


// read parent selection, generate offspring, eval, write
void var_app::state2()
{
  int *parents = io.read_sel();

  pop.new_gen(parents);

  diag << "Generation " << pop.cur_gen << std::endl;

  pop.variate();

  problem.eval(pop);
}

// write output
void var_app::dump_gen()
{
  dump_stream << "### GENERATION " << pop.cur_gen << std::endl;

  for (int id = pop.this_gen_begin();
       id < pop.this_gen_end(); ++id) {
    individual& indiv = pop.get(id);

    dump_stream << id << " ";

    for (int obj = 0; obj < pop.num_obj; ++obj)
      dump_stream << indiv.obj[obj] << " ";

    for (int dec = 0; dec < pop.num_dec; ++dec)
      dump_stream << indiv.dec[dec] << " ";

    for (int con = 0; con < pop.num_con; ++con)
      dump_stream << indiv.con[con] << " ";

    dump_stream << std::endl;
  }

  dump_stream << "#" << std::endl;
}
