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

#include "common.h"
#include "pop.h"
#include "join.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <cfloat>
#include <utility>
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
	}
}

inline location dispersal(size_t x, size_t y, double mu) {
	double t = 2.0*M_PI*myrand.uniform01();
	double r = myrand.exponential(mu);
	return location(
		(size_t)((int)floor(t*cos(r)+0.5+x)),
		(size_t)((int)floor(t*sin(r)+0.5+y))
	);
}

void population::step(size_t x, size_t y) {
	//pick a random position to replace
	haplotype hd,hm;
	gamete_info gid, gim;
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
		gid = dad.gamete(hd);
		mutate(hd);
		if(!is_valid(mom,hd))
			continue;
		gim = mom.gamete(hm);
		mutate(hm);
		
		individual &off = get_buf(x,y);
		// update offspring counts in parents
		size_t u;
		if(gid.first == 0) {
			for(u=0;u<gid.second;++u)
				dad.hdad[u].poff++;
			for(;u<dad.hmom.size();++u)
				dad.hmom[u].poff++;
		} else {
			for(u=0;u<gid.second;++u)
				dad.hmom[u].poff++;
			for(;u<dad.hdad.size();++u)
				dad.hdad[u].poff++;		
		}
		if(gim.first == 0) {
			for(u=0;u<gim.second;++u)
				mom.hdad[u].soff++;
			for(;u<mom.hmom.size();++u)
				mom.hmom[u].soff++;
		} else {
			for(u=0;u<gim.second;++u)
				mom.hmom[u].soff++;
			for(;u<mom.hdad.size();++u)
				mom.hdad[u].soff++;		
		}
		// set new genome and stats; should already be zeroed
		std::swap(off.hdad,hd);
		std::swap(off.hmom,hm);
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
	size_t pos = myrand.uniform(h.size());
	h[pos].allele = mallele++;
	gmut = myrand.geometric(pmut);
}

struct popstats {
	popstats() :
		num_pollen2(0.0),num_seed2(0.0), num_off2(0.0), num_allele(),
		num_homo(0.0)
		{ }
	
	double num_pollen2;
	double num_seed2;
	double num_off2;
	typedef map<size_t,double> alleledb;
	alleledb num_allele;
	double num_homo;
};

template<typename T>
inline T sq(const T& t) {
	return t*t;
}

void population::printstats() const
{
	// per locus:
	//   var fec - effective population size (components)
	//   number of alleles
	//   effective number of alleles
	//   conditional inbreeding
	//   r^2 to s
	//   true IBD?  over two generations?
	// per ind:
	//   s^2
	popstats stats;
	double M = 2.0*get_height()*get_width();
	double Mh = get_height()*get_width();	
	size_t U = 1;
	for(size_t y = 0; y < get_height(); ++y)
	{
		for(size_t x = 0; x < get_width(); ++x)
		{
			const individual &I = get_buf(x,y);
			stats.num_seed2   += sq(I.hmom[0].soff) + sq(I.hdad[0].soff);
			stats.num_pollen2 += sq(I.hmom[0].poff) + sq(I.hdad[0].poff);
			stats.num_off2    += sq(I.hmom[0].soff+I.hmom[0].poff)
			                   + sq(I.hdad[0].soff+I.hdad[0].poff);
			stats.num_allele[I.hmom[0].allele] += 1.0;
			stats.num_allele[I.hdad[0].allele] += 1.0;
			if(I.hmom[0].allele == I.hdad[0].allele)
				stats.num_homo += 1.0;
		}
	}
	double dt = 0.0;
	for(popstats::alleledb::iterator it = stats.num_allele.begin();
		it != stats.num_allele.end(); ++it)
	{
		dt += sq(it->second/M);
	}
	
	cout << join("\t",
		stats.num_off2/M - 1.0,
		stats.num_pollen2/M - 0.25,
		stats.num_seed2/M - 0.25,
		1.0/dt, stats.num_allele.size(),
		stats.num_homo/Mh,
		(stats.num_homo/Mh - dt) / (1.0 - dt),
		(stats.num_homo/Mh - dt) / (- dt)
	) << endl;
	
	stats = popstats();
	for(size_t y = 0; y < get_height(); ++y)
	{
		for(size_t x = 0; x < get_width(); ++x)
		{
			const individual &I = get_buf(x,y);
			stats.num_seed2   += sq(I.hmom[U].soff) + sq(I.hdad[U].soff);
			stats.num_pollen2 += sq(I.hmom[U].poff) + sq(I.hdad[U].poff);
			stats.num_off2    += sq(I.hmom[U].soff+I.hmom[U].poff)
			                   + sq(I.hdad[U].soff+I.hdad[U].poff);
			stats.num_allele[I.hmom[U].allele] += 1.0;
			stats.num_allele[I.hdad[U].allele] += 1.0;
			if(I.hmom[0].allele == I.hdad[0].allele)
				stats.num_homo += 1.0;			
		}
	}
	dt = 0.0;
	for(popstats::alleledb::iterator it = stats.num_allele.begin();
		it != stats.num_allele.end(); ++it)
	{
		dt += sq(it->second/M);
	}
	cout << join("\t",
		stats.num_off2/M - 1.0,
		stats.num_pollen2/M - 0.25,
		stats.num_seed2/M - 0.25,
		1.0/dt, stats.num_allele.size(),
		stats.num_homo/Mh,
		(stats.num_homo/Mh - dt) / (1.0 - dt),
		(stats.num_homo/Mh - dt) / (- dt)
	) << endl;
}

/****************************************************************************
 * class individual                                                         *
 ****************************************************************************/
 
individual::individual(size_t m, size_t k) {
	gene g;
	g.allele = myrand.uniform(k);
	hdad.assign(m,g);
	g.allele = myrand.uniform(k);
	hmom.assign(m,g);
	g.allele = -1;
	
	pdad = make_pair(-1,-1);
	pmom = pdad;
}

void individual::randomize(size_t k) {
	gene g;
	g.allele = myrand.uniform(k);
	fill(hdad.begin(), hdad.end(), g);
	g.allele = myrand.uniform(k);
	fill(hmom.begin(), hmom.end(), g);
}

gamete_info individual::gamete(haplotype & h) const {

	unsigned int r = myrand.uniform();
	unsigned int w = r & 1; // start with which haplotype
	r >>= 1;
	unsigned int p = 0;
	while(r & 1) {
		r >>= 1;
		++p;
	}
	unsigned int m = hdad.size();
	p = (p >= m) ? 1 : m-p;
	if(w == 0)
	{
		h.assign(hdad.begin(), hdad.begin()+p);
		if(p != m)
			h.insert(h.end(), hmom.begin()+p, hmom.end());
	}
	else
	{
		h.assign(hmom.begin(), hmom.begin()+p);
		if(p != m)
			h.insert(h.end(), hdad.begin()+p, hdad.end());
	}
	foreach(gene &g, h) {
		g.soff = 0;
		g.poff = 0;
	}
	return make_pair(w,p);
}

