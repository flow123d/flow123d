/*!@file doxy_dummy_defs.hh Dummy definitions for doxygen (DO NOT INCLUDE!!)  */

#error Do not include.

namespace std
{
  template<class T> class vector { public: T data; };
  template<class T> class deque { public: T data; };
  template<class T> class list { public: T data; };
  template<class T> class slist { public: T data; };
  template<class Key, class Data> class map { public: Key key; Data data; };
};


namespace arma
{
  class vec3 {};
  class mat33 {};
};
