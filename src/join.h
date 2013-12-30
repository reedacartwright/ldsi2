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

#ifndef RACWARE_JOIN_H
#define RACWARE_JOIN_H

#include "boost/tuple/tuple.hpp"
#include <iostream>
#include <vector>

namespace racware {
using namespace boost::tuples;

template <typename SEP = null_type,
  typename T0 = null_type, typename T1 = null_type, typename T2 = null_type,
  typename T3 = null_type, typename T4 = null_type, typename T5 = null_type,
  typename T6 = null_type, typename T7 = null_type, typename T8 = null_type,
  typename T9 = null_type>
struct join_t;


template<typename CharT, typename Traits, typename SP, typename T>
inline std::basic_ostream<CharT, Traits>& join_op(std::basic_ostream<CharT, Traits>& os,
	const SP &sp, const std::vector<T> &t);
template<typename CharT, typename Traits, typename SP, typename T, std::size_t N>
inline std::basic_ostream<CharT, Traits>& join_op(std::basic_ostream<CharT, Traits>& os,
	const SP &sp, const T (&t)[N]);
template<typename CharT, typename Traits, typename SP, std::size_t N>
inline std::basic_ostream<CharT, Traits>& join_op(std::basic_ostream<CharT, Traits>& os,
	const SP &sp, const char (&t)[N]);


template<typename CharT, typename Traits, typename SP, typename T>
inline std::basic_ostream<CharT, Traits>&
join_op(std::basic_ostream<CharT, Traits>& os, const SP& /*s*/, const T& t)
{
	return  (os << t);
}

namespace detail {
template<typename SP>
struct base_t
{
	typedef typename boost::tuples::tuple<SP> data_t;
	data_t st;

	base_t(typename access_traits<SP>::parameter_type s) : st(s) { }

	typename access_traits<SP>::const_type sep() const {
		return st.get<0>();
	}

	template<typename CharT, typename Traits, typename T>
	inline std::basic_ostream<CharT, Traits>&
		print(std::basic_ostream<CharT, Traits>& os, 
		const T &t) const {
		return join_op(os, sep(), t);
	}

	template<typename CharT, typename Traits, typename T1>
	inline std::basic_ostream<CharT, Traits>&
		print_tuple(std::basic_ostream<CharT, Traits>& os, 
		const cons<T1,null_type> &t) const {
		return print(os, t.head);
	}
	
	template<typename CharT, typename Traits, typename T1, typename T2>
	inline std::basic_ostream<CharT, Traits>&
		print_tuple(std::basic_ostream<CharT, Traits>& os, 
		const cons<T1,T2> &t) const {
		print(os, t.head) << sep();
		return print_tuple(os, t.tail);
	}	

	template<typename CharT, typename Traits, typename T>
	inline std::basic_ostream<CharT, Traits>&
		print_elements(std::basic_ostream<CharT, Traits>& os, 
		const T &tb, const T &te) const {
		if(tb == te)
			return os;
		T it = tb;
		print(os, *(it++));
		while(it != te) {
			os << sep();
			print(os, *(it++));
		}
		return os;
	}
};

} // detail

template<typename SP, typename T0, typename T1, typename T2, typename T3,
typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
struct join_t : public detail::base_t<SP>
{	
	typedef typename detail::base_t<SP> inherited;
	typedef typename boost::tuples::tuple<T0,T1,T2,T3,T4,T5,T6,T7,T8,T9> data_t;
	
	data_t tt;

	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0)
	  : inherited(sp), tt(t0) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1)
	  : inherited(sp), tt(t0,t1) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1,
		   typename access_traits<T2>::parameter_type t2)
	  : inherited(sp), tt(t0,t1,t2) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1,
		   typename access_traits<T2>::parameter_type t2,
		   typename access_traits<T3>::parameter_type t3)
	  : inherited(sp), tt(t0,t1,t2,t3) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1,
		   typename access_traits<T2>::parameter_type t2,
		   typename access_traits<T3>::parameter_type t3,
		   typename access_traits<T4>::parameter_type t4)
	  : inherited(sp), tt(t0,t1,t2,t3,t4) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1,
		   typename access_traits<T2>::parameter_type t2,
		   typename access_traits<T3>::parameter_type t3,
		   typename access_traits<T4>::parameter_type t4,
		   typename access_traits<T5>::parameter_type t5)
	  : inherited(sp), tt(t0,t1,t2,t3,t4,t5) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1,
		   typename access_traits<T2>::parameter_type t2,
		   typename access_traits<T3>::parameter_type t3,
		   typename access_traits<T4>::parameter_type t4,
		   typename access_traits<T5>::parameter_type t5,
		   typename access_traits<T6>::parameter_type t6)
	  : inherited(sp), tt(t0,t1,t2,t3,t4,t5,t6) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1,
		   typename access_traits<T2>::parameter_type t2,
		   typename access_traits<T3>::parameter_type t3,
		   typename access_traits<T4>::parameter_type t4,
		   typename access_traits<T5>::parameter_type t5,
		   typename access_traits<T6>::parameter_type t6,
		   typename access_traits<T7>::parameter_type t7)
	  : inherited(sp), tt(t0,t1,t2,t3,t4,t5,t6,t7) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1,
		   typename access_traits<T2>::parameter_type t2,
		   typename access_traits<T3>::parameter_type t3,
		   typename access_traits<T4>::parameter_type t4,
		   typename access_traits<T5>::parameter_type t5,
		   typename access_traits<T6>::parameter_type t6,
		   typename access_traits<T7>::parameter_type t7,
		   typename access_traits<T8>::parameter_type t8)
	  : inherited(sp), tt(t0,t1,t2,t3,t4,t5,t6,t7,t8) { }
	join_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<T0>::parameter_type t0,
		   typename access_traits<T1>::parameter_type t1,
		   typename access_traits<T2>::parameter_type t2,
		   typename access_traits<T3>::parameter_type t3,
		   typename access_traits<T4>::parameter_type t4,
		   typename access_traits<T5>::parameter_type t5,
		   typename access_traits<T6>::parameter_type t6,
		   typename access_traits<T7>::parameter_type t7,
		   typename access_traits<T8>::parameter_type t8,
		   typename access_traits<T9>::parameter_type t9)
	  : inherited(sp), tt(t0,t1,t2,t3,t4,t5,t6,t7,t8,t9) { }

	template<typename CharT, typename Traits>
	inline std::basic_ostream<CharT, Traits>&
		operator() (std::basic_ostream<CharT, Traits>& os) const
	{
		inherited::print_tuple(os, tt);
		return os;
	}
};

namespace detail {

template < typename SP = null_type,
  typename T0 = null_type, typename T1 = null_type, typename T2 = null_type,
  typename T3 = null_type, typename T4 = null_type, typename T5 = null_type,
  typename T6 = null_type, typename T7 = null_type, typename T8 = null_type,
  typename T9 = null_type >
struct join_mapper {
	typedef join_t < typename make_tuple_traits<SP>::type,
		typename make_tuple_traits<T0>::type,
		typename make_tuple_traits<T1>::type,
		typename make_tuple_traits<T2>::type,
		typename make_tuple_traits<T3>::type,
		typename make_tuple_traits<T4>::type,
		typename make_tuple_traits<T5>::type,
		typename make_tuple_traits<T6>::type,
		typename make_tuple_traits<T7>::type,
		typename make_tuple_traits<T8>::type,
		typename make_tuple_traits<T9>::type > type;
};

} // namespace detail

template<typename SP, typename T0>
typename detail::join_mapper<SP,T0>::type
join(const SP& sp, const T0& t0) {
	typedef typename detail::join_mapper <SP,T0>::type t;
	return t(sp,t0);
}
template<typename SP, typename T0, typename T1>
typename detail::join_mapper<SP,T0,T1>::type
join(const SP& sp, const T0& t0, const T1& t1) {
	typedef typename detail::join_mapper <SP,T0,T1>::type t;
	return t(sp,t0,t1);
}
template<typename SP, typename T0, typename T1, typename T2>
typename detail::join_mapper<SP,T0,T1,T2>::type
join(const SP& sp, const T0& t0, const T1& t1, const T2& t2) {
	typedef typename detail::join_mapper <SP,T0,T1,T2>::type t;
	return t(sp,t0,t1,t2);
}
template<typename SP, typename T0, typename T1, typename T2, typename T3>
typename detail::join_mapper<SP,T0,T1,T2,T3>::type
join(const SP& sp, const T0& t0, const T1& t1, const T2& t2, const T3& t3) {
	typedef typename detail::join_mapper <SP,T0,T1,T2,T3>::type t;
	return t(sp,t0,t1,t2,t3);
}
template<typename SP, typename T0, typename T1, typename T2, typename T3,
typename T4>
typename detail::join_mapper<SP,T0,T1,T2,T3,T4>::type
join(const SP& sp, const T0& t0, const T1& t1, const T2& t2, const T3& t3, const T4& t4) {
	typedef typename detail::join_mapper <SP,T0,T1,T2,T3,T4>::type t;
	return t(sp,t0,t1,t2,t3,t4);
}
template<typename SP, typename T0, typename T1, typename T2, typename T3,
typename T4, typename T5>
typename detail::join_mapper<SP,T0,T1,T2,T3,T4,T5>::type
join(const SP& sp, const T0& t0, const T1& t1, const T2& t2, const T3& t3, const T4& t4,
	 const T5& t5) {
	typedef typename detail::join_mapper <SP,T0,T1,T2,T3,T4,T5>::type t;
	return t(sp,t0,t1,t2,t3,t4,t5);
}
template<typename SP, typename T0, typename T1, typename T2, typename T3,
typename T4, typename T5, typename T6>
typename detail::join_mapper<SP,T0,T1,T2,T3,T4,T5,T6>::type
join(const SP& sp, const T0& t0, const T1& t1, const T2& t2, const T3& t3, const T4& t4,
	 const T5& t5, const T6& t6) {
	typedef typename detail::join_mapper <SP,T0,T1,T2,T3,T4,T5,T6>::type t;
	return t(sp,t0,t1,t2,t3,t4,t5,t6);
}
template<typename SP, typename T0, typename T1, typename T2, typename T3,
typename T4, typename T5, typename T6, typename T7>
typename detail::join_mapper<SP,T0,T1,T2,T3,T4,T5,T6,T7>::type
join(const SP& sp, const T0& t0, const T1& t1, const T2& t2, const T3& t3, const T4& t4,
	 const T5& t5, const T6& t6, const T7& t7) {
	typedef typename detail::join_mapper <SP,T0,T1,T2,T3,T4,T5,T6,T7>::type t;
	return t(sp,t0,t1,t2,t3,t4,t5,t6,t7);
}
template<typename SP, typename T0, typename T1, typename T2, typename T3,
typename T4, typename T5, typename T6, typename T7, typename T8>
typename detail::join_mapper<SP,T0,T1,T2,T3,T4,T5,T6,T7,T8>::type
join(const SP& sp, const T0& t0, const T1& t1, const T2& t2, const T3& t3, const T4& t4,
	 const T5& t5, const T6& t6, const T7& t7, const T8& t8) {
	typedef typename detail::join_mapper <SP,T0,T1,T2,T3,T4,T5,T6,T7,T8>::type t;
	return t(sp,t0,t1,t2,t3,t4,t5,t6,t7,t8);
}
template<typename SP, typename T0, typename T1, typename T2, typename T3,
typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
typename detail::join_mapper<SP,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9>::type
join(const SP& sp, const T0& t0, const T1& t1, const T2& t2, const T3& t3, const T4& t4,
	 const T5& t5, const T6& t6, const T7& t7, const T8& t8, const T9& t9) {
	typedef typename detail::join_mapper <SP,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9>::type t;
	return t(sp,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9);
}

template<typename SP, typename TI>
class join_elements_t : public detail::base_t<SP>
{
public:
	typedef typename detail::base_t<SP> inherited;
	typedef typename boost::tuples::tuple<TI,TI> data_t;

	data_t tt;

	join_elements_t(typename access_traits<SP>::parameter_type sp,
		   typename access_traits<TI>::parameter_type ta,
		   typename access_traits<TI>::parameter_type tb)
	: inherited(sp), tt(ta,tb) { }


	template<typename CharT, typename Traits>
	inline std::basic_ostream<CharT, Traits>&
		operator() (std::basic_ostream<CharT, Traits>& os) const
	{
		return print_elements(os, tt.get<0>(), tt.get<1>());
	}
};

template<typename SP, typename TI>
join_elements_t<
	typename boost::tuples::make_tuple_traits<SP>::type,
	typename boost::tuples::make_tuple_traits<TI>::type>
join_elements(const SP& sp, const TI& ta, const TI& tb) {
	typedef join_elements_t<
		typename boost::tuples::make_tuple_traits<SP>::type,
		typename boost::tuples::make_tuple_traits<TI>::type> t;
	return t(sp,ta,tb);
}

template<typename CharT, typename Traits, typename SP, typename T>
inline std::basic_ostream<CharT, Traits>& join_op(std::basic_ostream<CharT, Traits>& os,
	const SP &sp, const std::vector<T> &t) {
	return  (os << join_elements(sp, t.begin(), t.end()));
}

template<typename CharT, typename Traits, typename SP, typename T, std::size_t N>
inline std::basic_ostream<CharT, Traits>& join_op(std::basic_ostream<CharT, Traits>& os,
	const SP &sp, const T (&t)[N]) {
	return  (os << join_elements(sp, &t[0], &t[0]+N));
}

template<typename CharT, typename Traits, typename SP, std::size_t N>
inline std::basic_ostream<CharT, Traits>& join_op(std::basic_ostream<CharT, Traits>& os,
	const SP &sp, const char (&t)[N]) {
	return  (os << t);
}

} // namespace racware

namespace std {
template<typename CharT, typename Traits, typename SP, typename T0>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP, T0> &j  )
{
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1> &j) {
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1,
	typename T2>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1,T2> &j) {
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1,
	typename T2, typename T3>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1,T2,T3> &j) {
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1,
	typename T2, typename T3, typename T4>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1,T2,T3,T4> &j) {
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1,
	typename T2, typename T3, typename T4, typename T5>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1,T2,T3,T4,T5> &j) {
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1,
	typename T2, typename T3, typename T4, typename T5, typename T6>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1,T2,T3,T4,T5,T6> &j) {
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1,
	typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1,T2,T3,T4,T5,T6,T7> &j) {
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1,
	typename T2, typename T3, typename T4, typename T5, typename T6, typename T7,
	typename T8>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1,T2,T3,T4,T5,T6,T7,T8> &j) {
	return j(os);
}
template<typename CharT, typename Traits, typename SP, typename T0, typename T1,
	typename T2, typename T3, typename T4, typename T5, typename T6, typename T7,
	typename T8, typename T9>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_t<SP,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9> &j) {
	return j(os);
}

template<typename CharT, typename Traits, typename SP, typename TI>
inline basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& os,
	const racware::join_elements_t<SP, TI> &j ) {
	return j(os);
}
} // namespace std

#endif

