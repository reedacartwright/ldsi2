/****************************************************************************
 *  Copyright (C) 2008  Reed A. Cartwright, PhD <reed@scit.us>              *
 *                                                                          *
 *  This program is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 ****************************************************************************/

// src/ldsi -n 50 -m 15 -k 20 -u 1e-5 -p 6 -s 1.5 -c psi -t 1000 -b 10 -x 500 -y 2500  > p06.tsv

#define _USE_MATH_DEFINES

#include "common.h"
#include "pop.h"
#include "join.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_roots.h>
#include <cmath>
#include <cfloat>
#include <utility>
#include <algorithm>
#include <map>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

using namespace std;
using namespace racware;

class gslrand
{
public:
	typedef unsigned int int_type;
	typedef int_type seed_type;
	typedef double real_type;
	gslrand() {
		//gsl_rng_env_setup();
		T = gsl_rng_mt19937;
		r = gsl_rng_alloc(T);
	}
	gslrand(seed_type s) {
		//gsl_rng_env_setup();
		T = gsl_rng_mt19937;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r,s);
	}
	virtual ~gslrand() {
		gsl_rng_free(r);
	}
	
	inline void seed(seed_type s) { gsl_rng_set(r, s); }
	inline int_type get() { return gsl_rng_get(r); }
	inline int_type operator()() { return get(); }
	
	inline int_type uniform() { return get(); }
	inline int_type uniform(int_type max) { return gsl_rng_uniform_int(r, max); } // [0,max)
	inline real_type uniform01() { return gsl_rng_uniform(r); }
	inline bool boolean(double p) { return (uniform01() < p); }
	inline int_type poisson(real_type mu) { return gsl_ran_poisson(r,mu); }
	inline real_type exponential(real_type mu) { return gsl_ran_exponential(r,mu); }
	inline int_type geometric(real_type p) { return gsl_ran_geometric(r,p); }
	inline void multinomial(int_type N, const std::vector<real_type> &p, std::vector<int_type> &n)
	{
		n.resize(p.size(), 0);
		gsl_ran_multinomial(r, p.size(), N, &p[0], &n[0]);
	}

	// assumes n < max/2 for effeciency
	inline void uniform_urn(int_type max, int_type n, std::vector<int_type> &u) {
		u.assign(max,max);
		for(int_type i=0;i<n && i < max;++i)
		{
			int_type r = uniform(max-i);
			int_type t = (u[i] == max) ? i : u[i];
			u[i] = (u[i+r] == max) ? i+r : u[i+r];
			u[i+r] = t;
		}
		u.erase(u.begin()+n, u.end());
	}



private:
	const gsl_rng_type *T;
	gsl_rng *r;
} myrand;

void set_rand_seed(unsigned int u) { myrand.seed(u); }


/****************************************************************************
 * class population                                                         *
 ****************************************************************************/

void population::initialize(size_t w, size_t h, size_t m, size_t k)
{
	if(k < 1) k = 1;
	width = w;
	height = h;
	markers = m+1;
	mallele = k;
	inds.clear();
	inds.reserve(width*height);
	individual I(markers);
	
	// Initialize Dominance
	vector<double> dbsi;
	for(size_t kk=0;kk<k&&compat == BSI;++k)
		dbsi.push_back(myrand.uniform());
	
	for(size_t i=0;i<width*height;++i) {
		I.id = 2*i;
		I.randomize(k);
		if(compat == BSI) {
			I.hdad[0].dom = dbsi[I.hdad[0]];
			I.hmom[0].dom = dbsi[I.hmom[0]];
		}
		inds.push_back(I);
	}
	inds_buf.assign(inds.begin(), inds.end());
}

void population::params(double u, double U, double s, double p, size_t c) {
		mu = u; smu = U; seed = s; pollen = p; compat = c;
		pmut = 1.0-std::pow(1.0-mu, (int)(markers-1))*(1.0-smu);
		psmut = smu/(smu+(markers-1)*mu);
		gmut = myrand.geometric(pmut);
		
}

void population::evolve(size_t g, size_t b) {
	//cerr << "Burnin phase of " << b << " generations..." << endl;
	for(size_t gg=0;gg<b;++gg) {
		for(size_t y=0;y<height;++y) {
			for(size_t x=0;x<width;++x) {
				step(x,y);
			}
		}
		// swap is faster than copy
		std::swap(inds_buf, inds);
	}
	//cerr << "Collection phase of " << g << " generations..." << endl;
	print_stats_header();
	for(size_t gg=0;gg<g;++gg) {
		for(size_t y=0;y<height;++y) {
			for(size_t x=0;x<width;++x) {
				step(x,y);
			}
		}
		// swap is faster than copy
		std::swap(inds_buf, inds);
		if(gg%sample_gen == 0)
			print_stats(b+gg);
	}
}

inline location dispersal(size_t x, size_t y, double mu) {
	double t = 2.0*M_PI*myrand.uniform01();
	double r = myrand.exponential(mu);
	return location(
		(size_t)((int)floor(r*cos(t)+0.5+x)),
		(size_t)((int)floor(r*sin(t)+0.5+y))
	);
}

void population::step(size_t x, size_t y) {
	individual &off = get(x,y, true);
	while(1) {
		// pick a compatible pair
		location momloc = dispersal(x,y,seed);
		// test if momloc is valid
		if(!is_valid(momloc))
			continue;
		location dadloc = dispersal(momloc.first,momloc.second,pollen);
		// test if dadloc is valid by itself and paired with mom
		if(!is_valid(dadloc) || !is_valid(momloc,dadloc))
			continue;
		individual &mom = get(momloc);
		individual &dad = get(dadloc);

		// check for compatibility based on mom and dad genotypes
		if(!is_valid(mom,dad))
			continue;
		// generate pollen and check compatibility with mom
		// use off.hdad as buffer
		dad.gamete(off.hdad);
		mutate(off.hdad);
		if(!is_valid(mom,off.hdad))
			continue;
		mom.gamete(off.hmom);
		mutate(off.hmom);
		
		// set parent locations
		off.pdad = dadloc;
		off.pmom = momloc;
		break;		
	}
}

// assume 0 or 1 mutation per haplotype
void population::mutate(haplotype &h)
{
	if(--gmut > 0)
		return;
	gmut = myrand.geometric(pmut);
	size_t pos = 0;
	if(!myrand.boolean(psmut))
		pos = 1+myrand.uniform(h.size()-1);
	h[pos] = mallele++;
	if(compat == BSI && pos == 0)
		h[pos].dom = myrand.uniform();
}

template<typename T>
inline T sq(const T& t) {
	return t*t;
}

struct theta_kn_params {
	theta_kn_params(double kk, double nn) : k(kk), n(nn) { }
	double k;
	double n;
};

double theta_kn(double th, void *params)
{
	theta_kn_params *p = static_cast<theta_kn_params*>(params);
	return th*(gsl_sf_psi_n(0,p->n+th)-gsl_sf_psi_n(0,th))-p->k;
}

double theta_kn_deriv(double th, void *params)
{
	theta_kn_params *p = static_cast<theta_kn_params*>(params);
	return gsl_sf_psi_n(0,p->n+th)-gsl_sf_psi_n(0,th) +
			th*(gsl_sf_psi_n(1,p->n+th)-gsl_sf_psi_n(1,th));
}

void theta_kn_fdf(double th, void *params, double *y, double *dy)
{
	theta_kn_params *p = static_cast<theta_kn_params*>(params);
	double d = gsl_sf_psi_n(0,p->n+th)-gsl_sf_psi_n(0,th);
	*y = th*d-p->k;
	*dy = d+th*(gsl_sf_psi_n(1,p->n+th)-gsl_sf_psi_n(1,th));
}

class theta_kn_solver {
public:
	theta_kn_solver(double kk, double nn) :
		params(kk,nn)
	{
		T = gsl_root_fdfsolver_newton;
		s = gsl_root_fdfsolver_alloc(T);
		FDF.f = &theta_kn;
		FDF.df = &theta_kn_deriv;
		FDF.fdf = &theta_kn_fdf;
		FDF.params = &params;
	}
	~theta_kn_solver() {
		gsl_root_fdfsolver_free(s);
	}
	
	double operator()(double x, int miter=100) {
		if(-DBL_EPSILON < params.k-1.0 && params.k-1.0 < DBL_EPSILON)
			return 0.0;
		gsl_root_fdfsolver_set(s, &FDF, x);
		int iter = 0;
		do {
			iter++;
			status = gsl_root_fdfsolver_iterate(s);
			double x0 = x;
			x = gsl_root_fdfsolver_root(s);
			status = gsl_root_test_delta(x, x0, 0, 1e-7);
		} while(status == GSL_CONTINUE && iter < miter);
		return x;
	}
	
	theta_kn_params params;
	int status;
	
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	gsl_function_fdf FDF;
};

/*
var of a single theta_ko est is about 

t / (  PolyGamma[0,n+t]-PolyGamma[0,t]
    + t*(PolyGamma[1,n+t]-PolyGamma[1,t])
    )

n is the number of sequences

Calculated from Fisher's Information

*/

struct popstats {
	typedef std::map<size_t,size_t> alleledb;
	alleledb num_allele;
	size_t num_homo;
	size_t num_ibd;
	size_t sum_dist2;

	popstats() : num_allele(), num_homo(0), num_ibd(0), sum_dist2(0)
		{ }
};

void population::print_stats(size_t g)
{
	size_t uM = get_height()*get_width();
	double M = 2.0*sample_size;

	vector<gslrand::int_type> v;
	myrand.uniform_urn(uM, sample_size, v);
	double s2 = 0.0;
	for(size_t m = 0; m < markers; ++m)
	{
		popstats stats;

		foreach(size_t p, v)
		{
			size_t x = p%get_width();
			size_t y = p/get_width();
			const individual &I = get(x,y);
			stats.num_allele[I.hmom[m]] += 1;
			stats.num_allele[I.hdad[m]] += 1;
			if(I.hmom[m] == I.hdad[m])
				stats.num_homo += 1;
			if(I.hmom[m].gpar == I.hdad[m].gpar)
				stats.num_ibd += 1;
			if(m == 0)
			{
				stats.sum_dist2 += sq(I.pdad.first-x)+sq(I.pdad.second-y);
				stats.sum_dist2 += sq(I.pmom.first-x)+sq(I.pmom.second-y);
			}
		}
		size_t dt = 0;
		foreach(popstats::alleledb::value_type &vv, stats.num_allele) {
			dt += sq(vv.second);
		}
		if(m == 0)
			s2 = 0.5*stats.sum_dist2/M;
		double ibd = stats.num_ibd/(1.0*sample_size);
		double hibd = stats.num_homo/(1.0*sample_size);
		double f = dt/sq(M);
		double Ke = 1.0/f;
		double theta_ke = Ke-1.0;
		double N_ke = 0.25*theta_ke/mu;
		double Nb_ke = 4.0*M_PI*s2*N_ke/(uM);
		
		double Ko = (double)stats.num_allele.size();
		theta_kn_solver ts(Ko, M);
		// could fail to converged, ignore--tehe.
		double theta_ko = ts(Ke-0.9);
		double N_ko = 0.25*theta_ko/mu;
		double Nb_ko = 4.0*M_PI*s2*N_ko/(uM);	
				
		cout << join("\t", g, m, ibd, hibd, Ko, Ke, s2) 
		     << "\t" << join("\t", theta_ko, N_ko, Nb_ko)
	         << "\t" << join("\t", theta_ke, N_ke, Nb_ke)
		     << endl;
	}
}

void population::print_stats_header() const
{
	cout << "Gen" "\t"
	        "Mark" "\t"
	        "Ibd" "\t"
	        "Hoz" "\t"
	        "Ko" "\t"
	        "Ke" "\t"
	        "S2" "\t"
	        "ThKo" "\t"
	        "NKo" "\t"
	        "NbKo" "\t"
	        "ThKe" "\t"
	        "NKe" "\t"
	        "NbKe"
	     << endl;
}

/****************************************************************************
 * class individual                                                         *
 ****************************************************************************/
 
individual::individual(size_t m, size_t k) {
	gene g;
	g = myrand.uniform(k);
	hdad.assign(m,g);
	g = myrand.uniform(k);
	hmom.assign(m,g);
	
	pdad = make_pair(-1,-1);
	pmom = pdad;
}

void individual::randomize(size_t k) {
	gene g;
	g = myrand.uniform(k);
	fill(hdad.begin(), hdad.end(), g);
	g = myrand.uniform(k);
	fill(hmom.begin(), hmom.end(), g);
}

gamete_info individual::gamete(haplotype & h) const {
	unsigned int r = myrand.uniform();
	unsigned int w = r & 1; // start with which haplotype
	unsigned int rw = (w+1)&1; //other haplotype
	r >>= 1;
	unsigned int p = 0;
	while(r & 1) {
		r >>= 1;
		++p;
	}
	unsigned int m = hdad.size();
	p = (p >= m) ? 1 : m-p;
	const haplotype *g[2] = {&hdad, &hmom};
	h.assign(g[w]->begin(), g[w]->begin()+p);
	if(p != m)
		h.insert(h.end(), g[rw]->begin()+p, g[rw]->end());
	size_t ii;
	for(ii=0;ii<p;++ii) {
		h[ii].gpar = h[ii].par;
		h[ii].par = id+w;
	}
	for(;ii<m;++ii) {
		h[ii].gpar = h[ii].par;
		h[ii].par = id+rw;
	}
	return make_pair(w,p);
}

