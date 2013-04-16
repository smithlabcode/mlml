#include <string>
#include <vector>
#include <iostream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include <gsl/gsl_sf_gamma.h>

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::max;
using std::min;



/* NOTATION: 
 * p_m: probability of mC
 * p_h: probability of hmC
 *
 * TAB-seq VARIABLES (for 5hmC)
 * h: counts of C-reads
 * g: counts of T-reads
 *
 * oxBS-seq VARIABLES (for 5mC)
 * m: counts of C-reads
 * l: counts of T-reads
 *
 * BS-seq VARIABLES (for both 5hmC and 5mC)
 * t: counts of C-reads
 * u: counts of T-reads
 *
 * LATENT VARIABLES:
 * t1 (k in expectation): Cs that are from 5mC in BS-seq
 * g1 (j in expectation): Ts that are from 5mC in TAB-seq
 */
 
/*
static double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}
*/

static double
log_L(const size_t h, const size_t g,
      const size_t m, const size_t l,
      const size_t u, const size_t t,
      const double p_h, const double p_m){
  double log_lkhd = gsl_sf_lnchoose(h+g, h) + gsl_sf_lnchoose(m+l, m) +
    gsl_sf_lnchoose(u+t, u);
  if (p_h>0) log_lkhd += h*log(p_h);
  if (p_h<1) log_lkhd += g*log(1-p_h);
  if(p_m>0) log_lkhd += m*log(p_m);
  if(p_m<1) log_lkhd += l*log(1-p_m);
  if(p_h+p_m <1) log_lkhd += u*log(1-p_h-p_m);
  if(p_h+p_m >0)log_lkhd += t*log(p_h+p_m);
  return(log_lkhd);
}

static void
get_start_point(const size_t t, const size_t u,
                const size_t m, const size_t l,
                const size_t h, const size_t g,
                const double tolerance,
                double &p_m, double &p_h) {
  if (t+u==0) {
    p_m = 1.0 * m / l;
    p_h = 1.0 * h / g;
  }
  else if (m+l==0) {
    p_h = 1.0 * h / g;
    p_m = 1.0 - p_h;
  }
  else {
    p_m = 1.0 * m / l;
    p_h = 1.0 - p_m;
  }
  p_m = max(tolerance, min(p_m, 1-tolerance));
  p_h = max(2*tolerance-p_m, min(p_h, 1-p_m-2*tolerance));
}

static void
expectation(const size_t a, const size_t x,
	    const double p, const double q,
	    vector<vector<double> > &coeff) {
  assert(p>0&&q>0);
  assert(p+q<1);
  double ln_p = log(p);
  double ln_q = log(q);
  double p_plus_q = log(p+q);
  double p_u = log(1-p-q);
  double non_q = log(1-q);
  coeff = vector<vector<double> >(x + 1, vector<double>(a + 1));
  for (size_t k = 0; k <= x; ++k)
    for (size_t j = 0; j <= a; ++j)
      coeff[k][j] = exp(gsl_sf_lnchoose(a, j) + 
			ln_q*(a - j)+ ln_p*j - p_plus_q*a + 
			gsl_sf_lnchoose(x, k) +	
			ln_p*k + p_u*(x - k)- non_q*x);
}

static double
maximization(const size_t x, const size_t y,
	     const size_t a, const size_t b,
	     const vector<vector<double> > &coeff){
  double num = y, denom = y+b;
  for (size_t k = 0; k <= x; ++k)
    for (size_t j = 0; j <= a; ++j) {
      num += coeff[k][j]*(a - j);
      denom += coeff[k][j]*(a + x - k - j);
    }
  return num/denom;
}

static double
update_p_m(const size_t x, const size_t y,
	   const size_t z, const size_t w,
	   const size_t a, const size_t b, 
	   const vector<vector<double> > &coeff){
  double num = z;
  for (size_t k = 0; k <= x; ++k)
    for (size_t j = 0; j <= a; ++j)
      num += coeff[k][j]*(k + j);
  return num/(a + b + x + y + z + w);
}

static void
expectation_maximization(const bool VERBOSE,
			 const size_t x, const size_t y,
			 const size_t z, const size_t w,
			 const size_t a, const size_t b,
			 const double tolerance, 
			 double &p, double &q){
  size_t iter = 0;
  double delta = std::numeric_limits<double>::max();
  do {
    vector<vector<double> > coeff;
    expectation( a, x, p, q, coeff);
    const double M = maximization(x, y, a, b, coeff);
    const double p_old = p, q_old = q;
    p = update_p_m(x, y, z, w, a, b, coeff);
    q = M*(1 - p); 
    p = max(tolerance, p);
    p = min(1 - tolerance, p);
    q = max(2 * tolerance - p,  min(q, 1 - p - 2 * tolerance));
    delta = max(fabs(p_old - p), fabs(q_old - q));
    
    if (VERBOSE){
      cerr << iter++ << '\t'
	   << "M=" << M << '\t'
	   << "p_m=" << p << '\t'
	   << "p_h=" << q << '\t'
	   << "delta=" << delta << '\t'
	   <<"log-likelihood=" << log_L(y,x,z,w,b,a,q, p) << endl;
    }
  }
  while (delta > tolerance);
}


// when REV= false : 'a' is the number of C reads, 
//               and 'b' is the number of T reads.
// when REV= true :  'a' is the number of T reads, 
//               and 'b' is the number of C reads.
static void
parse_line(const bool REV, const string &line, 
	   size_t &a, size_t &b, string &chr, size_t &pos) {
  
  std::istringstream is(line);
  is >> chr;
  is >> pos;

  double level = 0.0;
  string dummy, str_count;
  is >> dummy >> str_count >> level;
  size_t count = 
    atoi(str_count.substr(str_count.find_first_of(":") + 1).c_str());
  if (count>50) count = 50;

  if(REV){
    b = static_cast<size_t>(count*level);
    a = count - b;
  }else{
    a = static_cast<size_t>(count*level);
    b = count - a;
  }
}


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string oxbs_seq_file;
    string hydroxy_file;
    string bs_seq_file;
    string outfile; 
    static double tolerance = 1e-10;
    

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "", "");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("bsseq", 'u', "Name of input BS-Seq methcounts file",
		      false , bs_seq_file);
    opt_parse.add_opt("tabseq", 'h', "Name of input TAB-Seq methcounts file",
		      false , hydroxy_file);
    opt_parse.add_opt("oxbsseq", 'm', "Name of input oxBS-Seq methcounts file",
		      false , oxbs_seq_file);
    opt_parse.add_opt("tolerance", 't', "EM convergence threshold. Default 1e-10",
		      false , tolerance);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;

    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() >0) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if ( (oxbs_seq_file.empty() && hydroxy_file.empty() ) ||
	 (oxbs_seq_file.empty() && bs_seq_file.empty() ) ||
	 (bs_seq_file.empty() && hydroxy_file.empty() )) {
      cerr << "Please specify at least 2 bed files as input" << endl;
    }

    /****************** END COMMAND LINE OPTIONS *****************/
    std::ofstream out(outfile.empty() ? "/dev/stdout" : outfile.c_str());
	
    size_t h = 0;
    size_t g = 0;
    size_t m = 0;
    size_t l = 0;
    size_t u = 0;
    size_t t = 0;
    size_t x = 0, y = 0, z = 0, w = 0, a = 0, b= 0;

    if(!hydroxy_file.empty() && !bs_seq_file.empty() && !oxbs_seq_file.empty()){
      std::ifstream h_in(hydroxy_file.c_str());
      std::ifstream b_in(bs_seq_file.c_str());
      std::ifstream o_in(oxbs_seq_file.c_str());
      string hydroxy_line, bs_line, oxbs_line;
      string h_chr, b_chr, o_chr;
      size_t h_pos = 0, b_pos = 0, o_pos = 0;
      while (getline(h_in, hydroxy_line) && getline(b_in, bs_line) &&
	     getline(o_in, oxbs_line)) {  
	parse_line(false, hydroxy_line, h, g, h_chr, h_pos);
	parse_line(false, bs_line, t, u, b_chr, b_pos);
	parse_line(false, oxbs_line, m, l, o_chr, o_pos);
	assert(h_chr == b_chr && h_chr == o_chr && 
	       h_pos == o_pos && h_pos == b_pos);
	double p_m,p_h;
	if((h+g>0 && u+t >0 ) ||
	   (h+g>0 && m+l >0 ) ||
	   (m+l>0 && u+t >0 ) ){
	  x = g; y = h;
	  z = m; w = l;
	  a = t; b = u;
          get_start_point(t,u,m,l,h,g,tolerance,p_m,p_h);
	  expectation_maximization(VERBOSE, x, y, z, w, a, b, tolerance, p_m, p_h); 
	  if (p_h <= 2.0*tolerance) p_h = 0;
	  if (p_m <= 2.0*tolerance) p_m = 0;
	  if (p_m >= 1-2.0*tolerance) p_m = 1;
	  if (p_h >= 1-2.0*tolerance) p_h =1;
	  //cout << h << '\t' << g << '\t' << m << '\t' << l << '\t' << u << '\t' << t << endl;
	  //cout << p_h << '\t' << p_m << endl;
	  out << h_chr << '\t' << h_pos << '\t'
	      << h_pos +1 << '\t' << "mC:"<<p_m << '\t'
	      << p_h << '\t' << '+' << endl;
	}else{
	  out << h_chr << '\t' << h_pos << '\t'
	      << h_pos +1 << '\t' << "mC:nan\tnan\t+" << endl;
	}
      }
    }else{
      std::ifstream f_in;
      std::ifstream s_in;
      string f_line, s_line, f_chr, s_chr;
      size_t f_pos = 0, s_pos = 0;
      bool f_rev = false, s_rev = false;
      size_t x = 0, y=0, z=0, w=0;
      if(oxbs_seq_file.empty()){
	f_in.open(bs_seq_file.c_str());
	s_in.open(hydroxy_file.c_str());
      }else if(hydroxy_file.empty()){
	f_in.open(bs_seq_file.c_str());
	s_in.open(oxbs_seq_file.c_str());
      }else{
	f_rev = true; 
	f_in.open(hydroxy_file.c_str());
	s_in.open(oxbs_seq_file.c_str());
      }

      while (getline(f_in, f_line) && getline(s_in, s_line)) {  
	parse_line(f_rev, f_line, x, y, f_chr, f_pos);
	parse_line(s_rev, s_line, z, w, s_chr, s_pos);
	assert(f_chr == s_chr && f_pos == s_pos);
	double p = 0.5 - tolerance;
	double q = 0.5 - tolerance;
	if(x+y>0 && z+w >0){
	  expectation_maximization(VERBOSE,x, y, z, w, 0, 0, tolerance, p, q); 
	  if (p <= 2.0*tolerance) p = 0;
	  if (q <= 2.0*tolerance) q = 0;
	  if (p >= 1-2.0*tolerance) p = 1;
	  if (q >= 1-2.0*tolerance) q =1;
	  out << f_chr << '\t' << f_pos << '\t'
	      << f_pos +1 << '\t' << "mC:";
	  if(oxbs_seq_file.empty())
	    out << 1.0-p-q << '\t' << p << '\t' << '+' << endl;
	  else if (hydroxy_file.empty())
	    out << p << '\t' << 1.0 - p - q << '\t' << '+' << endl;
	  else 
	    out << p << '\t' << q << '\t' << '+' << endl;
	}else{
	  out << f_chr << '\t' << f_pos << '\t'
	      << f_pos +1 << "\tmC:nan\tnan\t+" << endl;
	}
      } 
    }
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
