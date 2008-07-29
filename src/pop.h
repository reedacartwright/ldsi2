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
#include <cmath>

class population;
class individual;

//typedef size_t gene;

struct gene {
	gene() : a((size_t)-1), par((size_t)-1), gpar((size_t)-1) { }
	// true info
	size_t a;
	// meta info;
	size_t par;
	size_t gpar;

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
	enum {NSI=0,PSI=1,GSI=2,SSI=3};
	typedef std::vector<individual> indy_vec;
	population() { }

	void initialize(size_t w, size_t h, size_t m, size_t k);
	void params(double u, double s, double p, size_t c);
	inline void stat_params(size_t x, size_t y) {
		sample_size = x;
		sample_gen = y;
	}
	void evolve(size_t g=1);
	void step(size_t x, size_t y);

	inline const individual& get(size_t x, size_t y) const {
		return inds[x+y*width];
	}
	inline individual & get(size_t x, size_t y) {
		return inds[x+y*width];
	}
	inline const individual& get_buf(size_t x, size_t y) const {
		return inds_buf[x+y*width];
	}
	inline individual & get_buf(size_t x, size_t y) {
		return inds_buf[x+y*width];
	}	
	
	inline bool is_valid(size_t x, size_t y) const {
		return (x < width && y < height);
	}
	inline bool is_valid(const individual &m, const haplotype &h) const {
		switch(compat) {
		case GSI:
			return (h[0] != m.hdad[0] && h[0] != m.hmom[0]);
		case SSI:
			return true; //optimization because is_valid(mom,dad) has already been called
		case NSI:
		case PSI:
		default:
			return true;
		};
		return true;
	}
	inline bool is_valid(const individual &m, const individual &d) const {
		switch(compat) {
		case GSI:
			return (d.hdad[0] != m.hdad[0] &&
			        d.hdad[0] != m.hmom[0])||
			       (d.hmom[0] != m.hdad[0] &&
			        d.hmom[0] != m.hmom[0]); 
		case SSI:
			return (d.hdad[0] != m.hdad[0] &&
			        d.hdad[0] != m.hmom[0])&&
			       (d.hmom[0] != m.hdad[0] &&
			        d.hmom[0] != m.hmom[0]); 
		case NSI:
		case PSI:
		default:
			return true;
		};
		return true;
	}
	
	void mutate(haplotype &h);
	
	inline size_t get_width() const { return width; }
	inline size_t get_height() const { return height; }
	
	void printstats(size_t g) const;
	
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
