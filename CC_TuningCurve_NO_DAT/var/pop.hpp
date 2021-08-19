#ifndef VARIATOR_POPULATION_HPP
#define VARIATOR_POPULATION_HPP

#include <vector>
#include <ostream>

struct individual {
  individual()
  : job_finished(false) {};
  double *dec, *con, *obj;
  bool job_finished;
};


class population {
public:
  population() :
    cur_gen(1), decvar_props(0) {}

  ~population()
  {
    if (decvar_props) delete[] decvar_props;
  }

  void init();
  void fill_first_gen();
  void new_gen(int *selected_ids);
  void variate();
  void shrink_ranges();

  int this_gen_begin()
  {
    return (cur_gen-1)*gen_size;
  }

  int this_gen_end()
  {
    return cur_gen*gen_size;
  }

  individual& get(int id)
  {
    return indivs[id];
  }

  individual& get_rel(int i)
  {
    return indivs[this_gen_begin() + i];
  }

  individual& get_relprev(int i)
  {
    return indivs[(cur_gen-2)*gen_size + i];
  }

  
  int gen_size;
  int cur_gen, max_gen;
  int num_dec, num_obj, num_con;


  struct interval {
    interval() {}
    interval(double lbound, double ubound)
      : lbound(lbound), ubound(ubound) {}

    double lbound, ubound;
  };

  typedef std::vector<interval> region;
  std::vector<region> ini_regions, regions;

  struct decvar_prop {
    // overlapping intervals are combined and sorted here
    std::vector<interval> intervals;
    double range;
    double prev_sigma;
    double prev_min_off, prev_max_off;
  };

  decvar_prop* decvar_props;	// array

  double prev_slope;		// for varying eta

  struct {
    double indiv_mut_prob;
    double var_mut_prob;
    double var_swap_prob;
    double indiv_recomb_prob;
    double var_recomb_prob;
    double eta_mut;
    double eta_recomb;
    double shrink_threshold;
    double shrink_margin;
    double eta_mut_vary_threshold;
    double eta_mut_vary_factor;
    double eta_recomb_vary_threshold;
    double eta_recomb_vary_factor;
  } variation_prop;

private:
  void init_intervals();
  void copy_individual(individual& src, individual& dest);
  bool valid_region(individual& indiv);
  void uniform_crossover(individual& ind1, individual& ind2);
  double raw_to_offset(double x, int decvar);
  double offset_to_raw(double x_off, int decvar);
  double get_mean_offset(int decvar);
  double get_offset_sigma(int decvar);
  double get_min_offset(int decvar);
  double get_max_offset(int decvar);
  double get_recomb_spread (double x1, double x2, double range);
  void mix(individual& ind1, individual& ind2);
  void mutate(individual& indiv);
  double get_mut_jiggle(double x_off, double range);
  double get_slope();
  void vary_eta();

  individual *indivs;
};

#endif // VARIATOR_POPULATION_HPP
