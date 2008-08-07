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

#ifndef POP_H
#define POP_H

#include <vector>
#include <map>
#include <cmath>

class population;
class individual;

struct gene {
	gene() : a((size_t)-1), par((size_t)-1), gpar((size_t)-1) { }
	// true info
	size_t a;
	// meta info;
	size_t par;
	size_t gpar;
	size_t dom;

	bool operator == (const gene & g) const {
		return a == g.a;
	}
	bool operator != (const gene & g) const {
		return a != g.a;
	}
	bool operator < (const gene & g) const {
		return a < g.a;
	}
	
	gene & operator = (size_t u) {
		a = u;
		return *this;
	}

	operator size_t() const {
		return a;
	}
};

typedef std::vector<gene> haplotype;
typedef std::pair<size_t,size_t> location;
typedef location gamete_info;

class individual {
public:
	individual(size_t m, size_t k=1);
	void randomize(size_t k=1);
	
	gamete_info gamete(haplotype & h) const;
	
	size_t id;

	haplotype hdad;
	haplotype hmom;
	
	location pdad;
	location pmom;
};

class population {
public:
	enum {NSI=0,PSI=1,GSI=2,SSI=3,BSI=4};
	typedef std::vector<individual> indy_vec;
	population() { }

	void initialize(size_t w, size_t h, size_t m, size_t k);
	void params(double u, double s, double p, size_t c);
	inline void stat_params(size_t x, size_t y) {
		sample_size = x;
		sample_gen = y;
	}
	void evolve(size_t g=1, size_t b=0);
	void step(size_t x, size_t y);

	inline size_t index(size_t x, size_t y) const {
		return x+y*width;
	}
	inline size_t index(const location &loc) const {
		return index(loc.first, loc.second);
	}
	inline location position(size_t n) const {
		return std::make_pair(n%width,n/width);
	}
	
	inline const individual& get(size_t n, bool b=false) const {
		return b ? inds_buf[n] : inds[n];
	}
	inline individual& get(size_t n, bool b=false) {
		return b ? inds_buf[n] : inds[n];
	}
	inline const individual& get(const location &loc, bool b=false) const {
		return get(index(loc),b);
	}
	inline individual& get(const location &loc, bool b=false) {
		return get(index(loc),b);
	}
	inline const individual& get(size_t x, size_t y, bool b=false) const {
		return get(index(x,y), b);
	}
	inline individual& get(size_t x, size_t y, bool b=false) {
		return get(index(x,y), b);
	}

	
	inline bool is_valid(const location &loc ) const {
		return (loc.first < width && loc.second < height);
	}
	inline bool is_valid(const location &loc1, const location &loc2 ) const {
		return !((compat == PSI || compat == SSI || compat == BSI) && loc1 == loc2);
	}
	inline bool is_valid(const individual &m, const haplotype &h) const {
		return !(compat == GSI && (h[0] == m.hdad[0] || h[0] == m.hmom[0]));
	}
	inline bool is_valid(const individual &m, const individual &d) const {
		if(compat == SSI)
			return (d.hdad[0] == m.hdad[0] ||
		 	        d.hdad[0] == m.hmom[0] ||
				    d.hmom[0] == m.hdad[0] ||
				    d.hmom[0] == m.hmom[0] );
		else if(compat == BSI)
			return (d.hdad[0].dom >= d.hmom[0].dom) ?
				(d.hdad[0] == m.hdad[0] || d.hdad[0] == m.hmom[0]) :
				(d.hmom[0] == m.hdad[0] || d.hmom[0] == m.hmom[0]) ;
		return true;
	}
	
	void mutate(haplotype &h);
	
	inline size_t get_width() const { return width; }
	inline size_t get_height() const { return height; }
	
	void print_stats(size_t g);
	void print_stats_header() const;
	
protected:
	indy_vec inds;
	indy_vec inds_buf;
	size_t width;
	size_t height;
	size_t markers;
	double mu;
	double seed;
	double pollen;
	size_t compat;
	size_t mallele;
	double pmut;
	size_t gmut;
	size_t sample_gen;
	size_t sample_size;
	
};
#endif
