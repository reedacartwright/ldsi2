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

// src/ldsi -n 10 -m 4 -k 6 -u 1e-4 -p 2 -s 1 -c ssi -t 10

#define _USE_MATH_DEFINES

#include "common.h"
#include "pop.h"
#include "join.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
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
	for(size_t i=0;i<width*height;++i) {
		I.id = 2*i;
		I.randomize(k);
		inds.push_back(I);
		I.randomize(k);
		inds_buf.push_back(I);
	}
}

void population::params(double u, double s, double p, size_t c) {
		mu = u; seed = s; pollen = p; compat = c;
		pmut = 1.0-std::pow(1-mu, (int)markers);
		gmut = myrand.geometric(pmut);
}

void population::evolve(size_t g) {
	for(size_t gg=0;gg<g;++gg) {
		for(size_t y=0;y<height;++y) {
			for(size_t x=0;x<width;++x) {
				step(x,y);
			}
		}
		// swap is faster than copy
		std::swap(inds_buf, inds);
		if(gg%sample_gen == 0)
			printstats(gg);
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
	individual &off = get_buf(x,y);
	while(1) {
		// pick a compatible pair
		location momloc = dispersal(x,y,seed);
		if(!is_valid(momloc.first,momloc.second))
			continue;
		location dadloc = dispersal(momloc.first,momloc.second,pollen);
		if((compat != NSI && dadloc == momloc) ||
			!is_valid(dadloc.first,dadloc.second))
			continue;
		individual &mom = get(momloc.first,momloc.second);
		individual &dad = get(dadloc.first,dadloc.second);

		// check for potential compatibility
		if(!is_valid(mom,dad))
			continue;
		// generate pollen and check compatibility
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
	size_t pos = myrand.uniform(h.size());
	h[pos] = mallele++;
}


struct popstats {
	typedef map<size_t,size_t> alleledb;
	alleledb num_allele;
	size_t num_homo;
	size_t num_ibd;
	size_t sum_dist2;

	popstats() : num_allele(), num_homo(0), num_ibd(0), sum_dist2(0)
		{ }
};

template<typename T>
inline T sq(const T& t) {
	return t*t;
}

void population::printstats(size_t g) const
{
	// per locus:
	//   number of alleles
	//   effective number of alleles
	//   true IBD?  over two generations?
	// per ind:
	//   s^2
	size_t uM = get_height()*get_width();
	double M = 2.0*sample_size;

	vector<size_t> v;
	myrand.uniform_urn(uM, sample_size, v);
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
		
		cout << join("\t", g, m, dt/(M*M),
			stats.num_homo/(1.0*sample_size),
			stats.num_ibd/(1.0*sample_size),
			M*M/dt, stats.num_allele.size(),
			0.5*stats.sum_dist2/M
		) << endl;
	}
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

