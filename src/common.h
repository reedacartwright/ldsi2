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

#ifndef RACWARE_COMMON_H
#define RACWARE_COMMON_H

#ifdef NDEBUG
#ifdef _MSC_VER
#pragma warning(disable: 4512 4505 4146 4702)
#endif
#endif

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#ifndef PACKAGE_STRING
#	define PACKAGE_STRING "racware pkg 1-CURRENT"
#endif

#ifndef PACKAGE_BUGREPORT
#	define PACKAGE_BUGREPORT "reed@scit.us"
#endif

#define VERSION_MSG PACKAGE_STRING "\n" \
	"    Copyright (C) 2008  Reed A. Cartwright, PhD <reed@scit.us>\n" \
	"    Report bugs to " PACKAGE_BUGREPORT

// Utility functions
#define CERRORR(err_msg) ((std::cerr << "ERROR: " << err_msg << std::endl), false);
#define CERROR(err_msg) (std::cerr << "ERROR: " << err_msg << std::endl);

#include <ostream>
#include <vector>
#include <iterator>
#include <ctime>

#include "join.h"

namespace std {

template<typename Tp, typename CharT, typename Traits>
basic_ostream<CharT, Traits>&
operator<<(basic_ostream<CharT, Traits>& os, const std::vector<Tp> &v)
{
	return os << racware::join(",", v);
}
} // namespace std

template<typename T>
inline void deleteif(T* &t)
{
	if(t != NULL)
		delete t;
	t = NULL;
}

template<typename T>
inline void deleteifa(T* &t)
{
	if(t != NULL)
		delete[] t;
	t = NULL;
}

template<size_t _N>
size_t key_switch(const std::string &ss, const std::string (&key)[_N])
{
	for(size_t i=0;i<_N;++i)
	{
		if(key[i].find(ss) == 0)
			return i;
	}
	return (size_t)-1;
}

template<size_t _N>
size_t key_switch(const char *cs, const char* (&key)[_N])
{
	for(size_t i=0;i<_N;++i)
	{
		if(strstr(key[i], cs) == key[i])
			return i;
	}
	return (size_t)-1;
}

template<size_t _N>
size_t key_switch2(const char *cs, const char* (&key)[_N])
{
	static const size_t NFND = (size_t)-1;
	size_t ret = NFND;
	for(size_t i=0;i<_N;++i)
	{
		if(strstr(key[i], cs) != key[i])
			continue;
		if(ret != NFND)
			return NFND; // ambigious key
		ret = i;
	}
	return ret;
}

#ifdef _MSC_VER
#	include <process.h>
#endif

inline unsigned int rand_seed_start()
{
	unsigned int u = static_cast<unsigned int>(time(NULL));
	unsigned int v = static_cast<unsigned int>(getpid());
	v += (v << 15) + (v >> 3); // Spread 5-decimal PID over 32-bit number
	return u ^ (v + 0x9e3779b9u + (u<<6) + (u>>2));
}

inline unsigned int rand_seed()
{
	static unsigned int u = rand_seed_start();
	return (u = u*1664525u + 1013904223u);
}

void set_rand_seed(unsigned int u);

#endif
