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

#ifndef RACAPP_H
#define RACAPP_H

#include <boost/program_options.hpp>
#include <vector>
#include <string>
#include <cfloat>
#include <iomanip>

namespace po = boost::program_options;

typedef std::string sstring;

class app  {
public:
	app(int argc, char *argv[]);
	virtual ~app() { }
	
	virtual int run();
	
	po::options_description desc;	
	
	struct arg_t {
	// use X-Macros to specify argument variables
#	define XCMD(lname, sname, desc, type, def) type _V(lname) ;
#	include "app.cmds"
#	undef XCMD
	};
	
	template<typename OS> void print_args(OS &os) const;
protected:
	arg_t arg;
};

template<typename T>
inline bool do_print_args(const T&) { return true;}

template<>
inline bool do_print_args(const std::string& s) {
	return !s.empty();
}

template<typename OS>
void app::print_args(OS &os) const
{
	os << std::setprecision(DBL_DIG);
	#define XCMD(lname, sname, desc, type, def) \
	if(do_print_args(arg._V(lname))) \
		os << _S(lname) << "=" << arg._V(lname) << std::endl;
	#include "app.cmds"
	#undef XCMD
}

#endif
