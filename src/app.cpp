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
 
#include "common.h"

#include <boost/format.hpp>

#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <string>

#include "app.h"
#include "pop.h"
#include "join.h"


using namespace racware;
using namespace std;

void test();

int main(int argc, char *argv[])
{

	int ret = EXIT_FAILURE;
	try {
		app app(argc, argv);
		ret = app.run();
	} catch(exception &e) {
		CERROR(e.what());
	}
	return ret;
}

namespace boost { namespace program_options {
template<>
typed_value<bool>* value() { return bool_switch(); }

template<>
typed_value<bool>* value(bool* v) { return bool_switch(v); }
}}

app::app(int argc, char* argv[]) : desc("Allowed Options")
{
	try {
		desc.add_options()
			#define XCMD(lname, sname, desc, type, def) ( \
				_S(lname) _IFD(sname, "," BOOST_PP_STRINGIZE sname), \
				po::value< type >(&arg._V(lname))->default_value(def), \
				desc )
			#include "app.cmds"
			#undef XCMD
			;
		po::variables_map vm;
		po::positional_options_description pdesc;
		pdesc.add("input", -1);
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
		po::notify(vm);
		
		if(!arg.read_args.empty())
		{
			if(arg.read_args == "-")
			{
				po::store(po::parse_config_file(cin, desc), vm);	
			}
			else
			{
				std::ifstream ifs(arg.read_args.c_str());
				if(!ifs.is_open())
				{
					string sse = "unable to open argument file ";
					sse += arg.read_args;
					throw std::runtime_error(sse);
				}
				po::store(po::parse_config_file(ifs, desc), vm);
			}
			po::notify(vm);
		}
	} catch (exception &e) {
			CERROR(e.what());
			throw std::runtime_error("unable to process command line");
	}
}

// must be the same as the enum in population
const char *sikeys[] = {"nsi","psi","gsi","ssi", "bsi"};

int app::run()
{
	if(arg.help)
	{
		cerr << VERSION_MSG << endl << endl;
		cerr << desc << endl;
		return EXIT_SUCCESS;
	}
	if(arg.version)
	{
		cerr << VERSION_MSG << endl << endl;
		return EXIT_SUCCESS;
	}
		
	while(arg.rand_seed == 0)
		arg.rand_seed = rand_seed();
	set_rand_seed(arg.rand_seed);
	std::cerr << "Simulating with seed " << arg.rand_seed << endl;

	size_t ckey = key_switch(arg.compat.c_str(), sikeys);
	if(ckey == (size_t)-1)
	{
		CERROR("Invalid compatibility switch '" << arg.compat << "'.  Use (psi|nsi|gsi|ssi|bsi).");
		return EXIT_FAILURE;
	}
	
	population pop;
	pop.initialize(arg.size, arg.size, arg.markers, arg.ini_num);
	pop.params(arg.mu, arg.seed, arg.pollen, ckey);
	pop.stat_params(arg.sample, arg.step);
	pop.evolve(static_cast<size_t>(arg.time*2.0*arg.size*arg.size+0.5),
		static_cast<size_t>(arg.burn*2.0*arg.size*arg.size+0.5));
	
	return EXIT_SUCCESS;
}


