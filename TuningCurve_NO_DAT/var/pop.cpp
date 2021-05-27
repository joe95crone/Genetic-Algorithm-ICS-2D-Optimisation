#include "pop.hpp"
#include "misc.hpp"

#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <algorithm>

static double drand()
{
  return double(rand())/double(RAND_MAX);
}

static double drand(double low, double high)
{
  return low + drand()*(high-low);
}

// for sorting
bool operator<(const population::interval& i1, const population::interval& i2)
{
  return i1.lbound < i2.lbound;
}

void population::init()
{
  // allocate individuals and variable arrays
  indivs = new individual[max_gen * gen_size];

  for (int i = 0; i < max_gen * gen_size; ++i) {
    indivs[i].dec = new double[num_dec];
    indivs[i].obj = new double[num_obj];
    indivs[i].con = new double[num_con];
  }


  decvar_props = new decvar_prop[num_dec];

  regions.resize(ini_regions.size());
  std::copy(ini_regions.begin(), ini_regions.end(), regions.begin());
  init_intervals();
}


void population::init_intervals()
{
  for (int v = 0; v < num_dec; ++v) {
    decvar_prop& prop = decvar_props[v];
    prop.intervals.clear();
    prop.range = 0;

    std::vector<interval> intervals;

    for (int r = 0, n_r = regions.size(); r < n_r; ++r)
      intervals.push_back(regions[r][v]);

    std::sort(intervals.begin(), intervals.end()); // sort by lbound

    for (int i = 0, n_i = intervals.size(); i < n_i;) {
      double min_lbound = intervals[i].lbound;
      double max_ubound = intervals[i].ubound;

      // next one overlaps?
      while (i < n_i && intervals[++i].lbound < max_ubound)
	max_ubound = std::max(intervals[i].ubound, max_ubound);

      prop.intervals.push_back( interval(min_lbound, max_ubound) );
      prop.range += max_ubound - min_lbound;
    }

    prop.prev_min_off = 0;
    prop.prev_max_off = prop.range;
  }
}


void population::copy_individual(individual& src, individual& dest)
{
  std::copy(src.dec, src.dec + num_dec, dest.dec);
  std::copy(src.obj, src.obj + num_obj, dest.obj);
  std::copy(src.con, src.con + num_con, dest.con);

  dest.job_finished = false;
}


void population::fill_first_gen()
{
  // store cumulative region volumes (i.e., probabilities) here
  double accum_reg_vol[regions.size()];

  double tot_vol = 0;
  for (int r = 0; r < regions.size(); ++r) {
    double vol = 1;
    for (int v = 0; v < num_dec; ++v)
      vol *= regions[r][v].ubound - regions[r][v].lbound;
    tot_vol += vol;
    accum_reg_vol[r] = tot_vol;
  }
  

  for (int i = 0; i < gen_size; ++i) {
    double u = drand();
    for (int r = 0; r < regions.size(); ++r)
      if (u <= accum_reg_vol[r]/tot_vol) {
	for (int v = 0; v < num_dec; ++v)
	  indivs[i].dec[v] = drand(regions[r][v].lbound,
				   regions[r][v].ubound);
	break;
      }
  }

  for (int v = 0; v < num_dec; ++v)
    decvar_props[v].prev_sigma = get_offset_sigma(v);
}


void population::new_gen(int* selected_ids)
{
  ++cur_gen;

  for (int i = 0; i < gen_size; ++i) {
    individual& dest = get_rel(i);
    individual& src = get(selected_ids[i]);
    copy_individual(src, dest);
  }
}


void population::variate()
{
  shrink_ranges();
   vary_eta();

  // recombination

  int num_to_recomb = gen_size & ~1; // make even

  for (int i = 0; i < num_to_recomb; i += 2) {
    if (drand() <= variation_prop.var_swap_prob)
      uniform_crossover(get_rel(i), get_rel(i+1));

    if (drand() <= variation_prop.indiv_recomb_prob)
      mix(get_rel(i), get_rel(i+1));
  }


  // mutation

  for (int i = 0; i < gen_size; ++i) {
    if (drand() <= variation_prop.indiv_mut_prob)
      mutate(get_rel(i));
  }
}


bool population::valid_region(individual& indiv)
{
  for (int r = 0, n_r = regions.size(); r < n_r; ++r) {
    for (int v = 0; v < num_dec; ++v)
      if (indiv.dec[v] < regions[r][v].lbound
	  || indiv.dec[v] > regions[r][v].ubound)
	goto next_reg;
    // this point only reached if indiv lies in regions[r]
    return true;
next_reg:;
  }

  return false;
}


void population::uniform_crossover
(individual& ind1, individual& ind2)
{
  double saved_dec1[num_dec];
  double saved_dec2[num_dec];
  std::copy(ind1.dec, ind1.dec + num_dec, saved_dec1);
  std::copy(ind2.dec, ind2.dec + num_dec, saved_dec2);

  while (true) {
    for (int v = 0; v < num_dec; ++v)
      if (drand() < 0.5)
	std::swap(ind1.dec[v], ind2.dec[v]);

    if (valid_region(ind1) && valid_region(ind2))
      break;

    // rewind
    std::copy(saved_dec1, saved_dec1 + num_dec, ind1.dec);
    std::copy(saved_dec2, saved_dec2 + num_dec, ind2.dec);
    // and repeat...
  }
}


double population::raw_to_offset(double x, int decvar)
{
  std::vector<interval>& intervals = decvar_props[decvar].intervals;
  double offset = 0;
  for (int i = 0, n_i = intervals.size(); i < n_i; ++i)
    if (x >= intervals[i].lbound && x <= intervals[i].ubound)
      return offset + x - intervals[i].lbound;
    else
      offset += intervals[i].ubound - intervals[i].lbound;
}


double population::offset_to_raw(double x_off, int decvar)
{
  std::vector<interval>& intervals = decvar_props[decvar].intervals;
  double rem = x_off;
  for (int i = 0, n_i = intervals.size(); i < n_i; ++i)
    if (rem <= intervals[i].ubound - intervals[i].lbound)
      return intervals[i].lbound + rem;
    else
      rem -= intervals[i].ubound - intervals[i].lbound;
}


double population::get_mean_offset(int decvar)
{
  double sum = 0;

  for (int i = 0; i < gen_size; ++i) {
    individual& indiv = get_rel(i);
    sum += raw_to_offset(indiv.dec[decvar], decvar);
  }

  return sum / gen_size;
}


double population::get_offset_sigma(int decvar)
{
  double mean = get_mean_offset(decvar);
  double sum = 0;

  for (int i = 0; i < gen_size; ++i) {
    individual& indiv = get_rel(i);
    double diff = raw_to_offset(indiv.dec[decvar], decvar) - mean;
    sum += diff * diff;
  }

  return sqrt(sum / gen_size);
}


double population::get_min_offset(int decvar)
{
  double min = DBL_MAX;

  for (int i = 0; i < gen_size; ++i) {
    individual& indiv = get_rel(i);
    min = std::min(min, raw_to_offset(indiv.dec[decvar], decvar));
  }

  return min;
}


double population::get_max_offset(int decvar)
{
  double max = 0;

  for (int i = 0; i < gen_size; ++i) {
    individual& indiv = get_rel(i);
    max = std::max(max, raw_to_offset(indiv.dec[decvar], decvar));
  }

  return max;
}


void population::shrink_ranges()
{
  regions.resize(ini_regions.size());
  std::copy(ini_regions.begin(), ini_regions.end(), regions.begin());
  init_intervals();		// restore initial intervals to calculate sigma

  for (int v = 0; v < num_dec; ++v) {
    diag << "Resize decvar #" << v << "?" << std::endl;
    decvar_prop& prop = decvar_props[v];
    double range = prop.range;
    double sigma = get_offset_sigma(v);

    diag << "sigma=" << sigma << std::endl;

    double min_off = get_min_offset(v);
    min_off -= min_off * variation_prop.shrink_margin;

    double max_off = get_max_offset(v);
    max_off += (range - max_off) * variation_prop.shrink_margin;

    diag << "min_off=" << min_off << ", "
	 << "max_off=" << max_off  << std::endl;

    if (sigma < variation_prop.shrink_threshold * prop.prev_sigma) {
      // we're resizing (hopefully, shrinking) the window - store this
      // sigma for comparison during next gen, and save the window
      // edges
      diag << "SHRINKING" << std::endl;
      prop.prev_sigma = sigma;
      prop.prev_min_off = min_off;
      prop.prev_max_off = max_off;
    }
    else {
      // we're not supposed to resize, but we'll widen (yes, widen)
      // the window if necessary to accommodate individuals revived
      // from past generations. if widening isn't necessary, use
      // window from previous resizing
      min_off = std::min(min_off, prop.prev_min_off);
      max_off = std::max(max_off, prop.prev_max_off);
      diag << "NOT SHRINKING: "
	   << "min_off=" << min_off << ", "
	   << "max_off=" << max_off << std::endl;
    }

    std::vector<region>::iterator
      it = regions.begin(), end = regions.end();

    // remove regions that lie completely outside the window, and
    // shrink those that overflow its bounds
    int i=0; // dbg
    while (it != end) {
      std::vector<region>::iterator here = it++;
      region& r = *here;
      double lbound_off = raw_to_offset(r[v].lbound, v);
      double ubound_off = raw_to_offset(r[v].ubound, v);

      if (ubound_off < min_off || lbound_off > max_off) {
	regions.erase(here);
	diag << "removing region " << i << std::endl; // dbg
      }
      else {
	if (lbound_off < min_off)
	  r[v].lbound = offset_to_raw(min_off, v);
	if (ubound_off > max_off)
	  r[v].ubound = offset_to_raw(max_off, v);
      }
      ++i;			// dbg
    }
  }

  init_intervals();		// generate new (shrunken) intervals
				// for variate()

  diag << std::endl;
}

// this is called after new_gen(), so we look at generations cur_gen-1
// and cur_gen-2
double population::get_slope()
{
  double this_avg = 0, prev_avg = 0;

  for (int i = 0; i < gen_size; ++i) {
    // (cur_gen-1)-1 b/c array indexing starts at 0
    this_avg += indivs[(cur_gen-2)*gen_size + i].obj[0]; 
    prev_avg += indivs[(cur_gen-3)*gen_size + i].obj[0];
  }

  this_avg /= gen_size;
  prev_avg /= gen_size;

  return log(this_avg) - log(prev_avg);
}


// called after new_gen() but before mutation/recombination
void population::vary_eta()
{
  if (cur_gen < 3) return;	// we need two generations of data

  double slope = get_slope();

  if (cur_gen > 3) {		// we need an old slope to compare to
    double mut_thresh = variation_prop.eta_mut_vary_threshold;
    if (mut_thresh > 1e-20 // safer than comparing to zero
	&& -slope < -prev_slope * mut_thresh) {
      variation_prop.eta_mut *= variation_prop.eta_mut_vary_factor;
      diag << "eta_mutation changed to " << variation_prop.eta_mut
	   << std::endl;
    }

    double recomb_thresh = variation_prop.eta_recomb_vary_threshold;
    if (recomb_thresh > 1e-20
	&& -slope < -prev_slope * recomb_thresh) {
      variation_prop.eta_recomb *= variation_prop.eta_recomb_vary_factor;
      diag << "eta_recombination changed to " << variation_prop.eta_recomb
	   << std::endl;
    }
  }

  prev_slope = slope;
}


void population::mix(individual& ind1, individual& ind2)
{
  double saved_dec1[num_dec];
  double saved_dec2[num_dec];
  std::copy(ind1.dec, ind1.dec + num_dec, saved_dec1);
  std::copy(ind2.dec, ind2.dec + num_dec, saved_dec2);
  
  while (true) {		// attempt until we get valid offspring
    for (int v = 0; v < num_dec; ++v) {
      if (drand() < variation_prop.var_recomb_prob) {
	double x1 = ind1.dec[v];
	double x2 = ind2.dec[v];

	bool swapped = x1 > x2;
	if (swapped) std::swap(x1, x2); // now x1 < x2

	double x1_off = raw_to_offset(x1, v);
	double x2_off = raw_to_offset(x2, v);

	double avg_off = 0.5*(x1_off+x2_off);
	double spread = get_recomb_spread(x1_off, x2_off, decvar_props[v].range);

	// note: +/- were switched before!
	x1 = offset_to_raw(avg_off - spread, v);
	x2 = offset_to_raw(avg_off + spread, v);

	if (swapped) std::swap(x1, x2);
	ind1.dec[v] = x1;
	ind2.dec[v] = x2;
      }
    }

    if (valid_region(ind1) && valid_region(ind2))
      break;			// yay

    // rewind
    std::copy(saved_dec1, saved_dec1 + num_dec, ind1.dec);
    std::copy(saved_dec2, saved_dec2 + num_dec, ind2.dec);
    // and repeat...
  }
}


double population::get_recomb_spread
(double x1, double x2, double range) // expects x1 < x2
{
  double eta = variation_prop.eta_recomb;

  double p_b;

  double dx = x2 - x1;
  if (dx != 0) {
    double bl = 1 + 2*x1/dx;
    double bu = 1 + 2*(range-x2)/dx;

    double bmin = std::min(bl, bu);

    p_b = 1 - 1/(2 * pow(bmin, eta+1));
  } else
    p_b = 1;			// shouldn't happen :)

  double u = drand() * p_b;
  double b;
      
  if (u < 0.5)
    b = pow(2*u, 1/(eta+1));
  else
    b = pow(0.5/(1-u), 1/(eta+1));

  return b*dx/2;
}


void population::mutate(individual& indiv)
{
  double saved_dec[num_dec];
  std::copy(indiv.dec, indiv.dec + num_dec, saved_dec);

  while (true) {
    for (int v = 0; v < num_dec; ++v)
      if (drand() <= variation_prop.var_mut_prob) {
	double x_off = raw_to_offset(indiv.dec[v], v);

	double jiggle = get_mut_jiggle(x_off, decvar_props[v].range);
	indiv.dec[v] = offset_to_raw(x_off + jiggle, v);
      }

    if (valid_region(indiv))
      break;

    // rewind
    std::copy(saved_dec, saved_dec + num_dec, indiv.dec);
    // and repeat...
  }
}


double population::get_mut_jiggle(double x_off, double range)
{
  double rel_max = std::min(x_off, range-x_off) / range; // [0, 0.5]
  double eta = variation_prop.eta_mut;

  double u = drand();
  int sign = 1;

  if (u > 0.5) {
    u -= 0.5;
    sign = -1;
  }

  // u in [0, 0.5]
  double b = 2*u + (1-2*u) * pow(1-rel_max, eta+1);
  double rel_delta = pow(b, 1/(eta+1)) - 1;

  return sign * range * rel_delta;
}
