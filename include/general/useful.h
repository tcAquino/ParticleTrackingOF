//
// useful.h
//
// Created on: Mar 15, 2011
// Author: tomas
//

// Miscelaneous collection of useful objects and algorithms

#ifndef USEFUL_H_
#define USEFUL_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <ios>
#include <iostream>
#include <list>
#include <regex>
#include <stdexcept>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace useful
{
  // Convert string to bool
  bool stob(std::string const& string)
  {
    if (string == "true" || string == "1")
      return true;
    if (string == "false" || string == "0")
      return false;
    throw std::runtime_error{
      "Expected true or false, got"
      + string };
  }
  
  // Convert bool to string
  std::string btos(bool val, bool boolalpha = 1)
  {
    std::ostringstream stream;
    if (boolalpha)
      stream << std::boolalpha;
    else
      stream << std::noboolalpha;
    stream << val;
    
    return stream.str();
  }
  
  // Check if string is empty
  // Strings are considered empty if they hold:
  // - Nothing
  // - ''
  // - ""
  bool empty(std::string const& str)
  {
    return str == ""
      || str == "''"
      || str == R"("")";
  }
  
  // Change extension after last dot
  // Note: new_extension should include dot if wanted
  std::string change_extension
  (std::string const& filename,
   std::string const& new_extension)
  {
    return filename.substr(0, filename.find_last_of('.'))
      + new_extension;
  }
  
  // Expand environment variables
  // Based on Toby Speight's answer here:
  // https://codereview.stackexchange.com/questions/172644/c-environment-variable-expansion
  inline std::string expand_env(std::string text)
  {
    static const std::regex env_re{R"--(\$\{([^}]+)\})--"};
    std::smatch match;
    while (std::regex_search(text, match, env_re))
    {
      auto const from = match[0];
      text.replace(from.first, from.second,
                   std::getenv(match[1].str().c_str()));
    }
    return text;
  }
  
  // Get widths of bins having given values as midpoints
  // and given minimum (left edge)
  template <typename Value_type = double>
  std::vector<Value_type> get_bin_widths
  (std::vector<Value_type> const& values, Value_type minimum)
  {
    std::vector<Value_type> widths;
    widths.reserve(values.size());
    double left_edge = minimum;
    for (auto const& val : values)
    {
      widths.push_back(2.*(val-left_edge));
      if (widths.back() < 0.)
        throw std::runtime_error{ "Bad minimum value, could not determine PDF bin widths" };
      left_edge += widths.back();
    }
    
    return widths;
  }
  
  // Get widths of bins having given values as midpoints
  // and given minimum (left edge)
  template <typename Value_type = double>
  std::vector<Value_type> get_bin_widths
  (std::vector<Value_type> const& values)
  {
    if (values.size() < 2)
      throw std::runtime_error{ "Less than two values, could not determine PDF bin widths" };
    double midpoint = (values[1]+values[0])/2.;
    return get_bin_widths(values, values[0]-(midpoint-values[0]));
  }
  
  // Check if container contains val
  // Warning: container must be sorted
  template <typename T, typename U>
  bool contains(const std::vector<T>& container, U const& val)
  {
    auto it = std::lower_bound(
          container.begin(),
          container.end(),
          val,
          [](T const& elem, U const& val){ return elem < val; });
    
    return it != container.end() && *it == val;
  }
  
  // Split string
  // Adapted from Beder Acosta Borges's answer here:
  // https://stackoverflow.com/questions/14265581/
  // parse-split-a-string-in-c-using-string-delimiter-standard-c
  bool endsWith(const std::string& string, const std::string& suffix)
  {
    return string.size() >= suffix.size() &&
      string.substr(string.size() - suffix.size()) == suffix;
  }
  
  // Split string into vector of strings
  // Split at instances of delimeter
  // Empty entries can be included or discarded
  std::vector<std::string> split
  (std::string const& string, std::string const& delimiter = " ",
   bool empty_entries = false)
  {
    std::vector<std::string> tokens;

    for (std::size_t start = 0, end; start < string.length();
         start = end+delimiter.length())
    {
      std::size_t position = string.find(delimiter, start);
      end = position != std::string::npos? position : string.length();

      std::string token = string.substr(start, end-start);
      if (empty_entries || !token.empty())
        tokens.push_back(token);
    }

    if (empty_entries &&
        (string.empty() || endsWith(string, delimiter)))
      tokens.push_back("");

    return tokens;
  }
  
  // Throw exception for parsing a file line
  auto parse_error
  (std::string const& filename, std::string const& line)
  {
    return std::runtime_error{
      "Could not parse line\n" + line +
      "\nin file " + filename };
  }
  
  // Throw exception for parsing a file
  auto parse_error_file(std::string const& filename)
  {
    return std::runtime_error{
      "Could not parse file " + filename };
  }
  
  // Throw exception for parsing a line
  auto parse_error_line(std::string const& line)
  {
    return std::runtime_error{
      "Could not parse line\n" + line };
  }
  
  // Throw exception for opening a file for reading
  auto open_read_error(std::string const& filename)
  {
    return std::runtime_error{
      "Could not open file " + filename + " for reading" };
  }
  
  // Throw exception for opening a file for writing
  auto open_write_error(std::string const& filename)
  {
    return std::runtime_error{
      "Could not open file " + filename + " for writing" };
  }
  
  // Throw exception for unexpected contents in file
  auto bad_file_contents(std::string const& filename)
  {
    return std::runtime_error{
      "Innapropriate contents in file " + filename };
  }
  
  // Throw exception for finding end of a file
  // before expected string
  auto bad_eof
  (std::string const& filename, std::string const& string)
  {
    return std::runtime_error{
    "Reached end of " +
      filename + " before " + string + " was found" };
  }
  
  // Throw exception for bad parameters
  auto bad_parameters()
  {
    return std::invalid_argument{ "Inappropriate parameters" };
  }
  
  // Throw exception for bad parameters and suggest help
  auto bad_parameters_help()
  {
    return std::invalid_argument{ "Inappropriate parameters"
      " (-h or --help for help)" };
  }
  
  // Check command line options for help flag
  bool check_options_help(int argc, const char* const argv[])
  {
    return argc == 2 &&
      (std::string(argv[1]) == "--help" ||
       std::string(argv[1]) == "-h");
  }
  
  // Open file for reading
  std::ifstream open_read(std::string const& filename)
  {
    std::ifstream file(filename);
    if (!file.is_open())
      throw useful::open_read_error(filename);
    
    return file;
  }
  
  // Open file for writing
  std::ofstream open_write(std::string const& filename)
  {
    std::ofstream file(filename);
    if (!file.is_open())
      throw useful::open_write_error(filename);
    
    return file;
  }
  
  //Load 1-column file into vector of doubles
  auto load_1
  (std::string const& filename, std::size_t nr_estimate = 0,
   std::size_t header_lines = 0,
   std::string const& delim = " ")
  {
    using Value = double;
    using Container = std::vector<Value>;
    Container values;
    values.reserve(nr_estimate);
    
    auto file = open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line = split(line, delim);
      if (split_line.size() != 1)
        throw parse_error(filename, line);
      values.push_back(std::stod(split_line[0]));
    }
    file.close();
    
    return values;
  }
  
  //Load 2-column file into pair of vectors of doubles
  auto load_2
  (std::string const& filename, std::size_t nr_estimate = 0,
   std::size_t header_lines = 0, std::string const& delim = " ")
  {
    using Value = double;
    using Container = std::vector<Value>;
    std::pair<Container, Container> values;
    values.first.reserve(nr_estimate);
    values.second.reserve(nr_estimate);
    
    auto file = open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line = split(line, delim);
      if (split_line.size() != 2)
        throw parse_error(filename, line);
      
      values.first.push_back(std::stod(split_line[0]));
      values.second.push_back(std::stod(split_line[1]));
    }
    file.close();
    
    return values;
  }
  
  //Load file into vector of vectors of doubles
  auto load(std::string const& filename, std::size_t nr_columns,
            std::size_t nr_estimate = 0,
              std::size_t header_lines = 0,
              std::string const& delim = " ")
  {
    using Value = double;
    using Container = std::vector<Value>;
    std::vector<Container> values(nr_columns);
    for (auto& val : values)
      val.reserve(nr_estimate);
    
    auto file = open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line = split(line, delim);
      if (split_line.size() != nr_columns)
        throw parse_error(filename, line);
      
      for (std::size_t cc = 0; cc < nr_columns; ++cc)
        values[cc].push_back(std::stod(split_line[cc]));
    }
    file.close();
    
    return values;
  }
  
  // Read next value from file
  template <typename Type>
  void read(std::ifstream& input, Type& val)
  {
    input >> val;
    if (input.fail())
      throw std::runtime_error{ "Could not read value" };
  }
  
  // From Passer By's answer here:
  // https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain
  // Check whether class can call sqrt
  template <typename = void, typename... Args>
  struct can_call_sqrt : std::false_type {};
  template <typename... Args>
  struct can_call_sqrt<
  std::void_t<decltype(std::sqrt(std::declval<Args>()...))>, Args...>
  : std::true_type {};
  template <typename... Args>
  inline constexpr bool can_call_sqrt_v = can_call_sqrt<void, Args...>::value;
  // Check whether class can call abs
  template <typename = void, typename... Args>
  struct can_call_abs : std::false_type {};
  template <typename... Args>
  struct can_call_abs<
  std::void_t<decltype(std::abs(std::declval<Args>()...))>, Args...>
  : std::true_type {};
  template <typename... Args>
  inline constexpr bool can_call_abs_v = can_call_abs<void, Args...>::value;
  
  
  // From Richard Hodges's answer here:
  // https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain
  // Check whether classes can be used with binary operator
  template<class X, class Y, class Op>
  struct op_valid_impl
  {
      template<class U, class L, class R>
      static auto test(int) -> decltype(std::declval<U>()(std::declval<L>(), std::declval<R>()),
                                        void(), std::true_type());

      template<class U, class L, class R>
      static auto test(...) -> std::false_type;

      using type = decltype(test<Op, X, Y>(0));

  };
  template<class X, class Y, class Op> using op_valid = typename op_valid_impl<X, Y, Op>::type;
  // Operators not in std::
  namespace notstd {

      struct left_shift {

          template <class L, class R>
          constexpr auto operator()(L&& l, R&& r) const
          noexcept(noexcept(std::forward<L>(l) << std::forward<R>(r)))
          -> decltype(std::forward<L>(l) << std::forward<R>(r))
          {
              return std::forward<L>(l) << std::forward<R>(r);
          }
      };

      struct right_shift {

          template <class L, class R>
          constexpr auto operator()(L&& l, R&& r) const
          noexcept(noexcept(std::forward<L>(l) >> std::forward<R>(r)))
          -> decltype(std::forward<L>(l) >> std::forward<R>(r))
          {
              return std::forward<L>(l) >> std::forward<R>(r);
          }
      };

  }
  template<class X, class Y> using has_equality = op_valid<X, Y, std::equal_to<>>;
  template<class X, class Y> using has_inequality = op_valid<X, Y, std::not_equal_to<>>;
  template<class X, class Y> using has_less_than = op_valid<X, Y, std::less<>>;
  template<class X, class Y> using has_less_equal = op_valid<X, Y, std::less_equal<>>;
  template<class X, class Y> using has_greater_than = op_valid<X, Y, std::greater<>>;
  template<class X, class Y> using has_greater_equal = op_valid<X, Y, std::greater_equal<>>;
  template<class X, class Y> using has_bit_xor = op_valid<X, Y, std::bit_xor<>>;
  template<class X, class Y> using has_bit_or = op_valid<X, Y, std::bit_or<>>;
  template<class X, class Y> using has_left_shift = op_valid<X, Y, notstd::left_shift>;
  template<class X, class Y> using has_right_shift = op_valid<X, Y, notstd::right_shift>;
  template<class X, class Y> using has_plus = op_valid<X, Y, std::plus<>>;
  template<class X, class Y> using has_minus = op_valid<X, Y, std::minus<>>;
  template<class X, class Y> using has_multiplies = op_valid<X, Y, std::multiplies<>>;
  template<class X, class Y> using has_divides = op_valid<X, Y, std::divides<>>;
  
  // Check whether container contains value
  // Container must be sorted according to comp_less,
  // equality is checked using comp_eq
  template <typename T, typename U, typename Comp_less, typename Comp_eq>
  bool contains
  (std::vector<T> const& container,
   U const& val,
   Comp_less comp_less = std::less<void>{},
   Comp_eq comp_eq = std::equal_to<void>{})
  {
    auto it = std::lower_bound(
          container.begin(),
          container.end(),
          val,
          comp_less);
    
    return it != container.end() && comp_eq(*it, val);
  }

  // Sign of val
  template <typename T> int sgn(T val)
  { return (T(0) < val) - (val < T(0)); }

  // Class holding const object
  template <typename Object_Type, typename Return_Type = Object_Type>
  struct StoreConst
  {
    using value_type = Return_Type;
    
    Return_Type operator()() const
    { return obj; }
    const Object_Type obj;
  };

  // Class holding object
  template <typename Object_Type, typename Return_Type = Object_Type>
  struct Store
  {
    using value_type = Return_Type;
    
    Return_Type operator()() const
    { return obj; }
    Object_Type obj;
  };

  // Class holding nothing
  struct Empty
  {
    template <typename ...Args>
    Empty(Args...){}
    Empty(){}
  };
  
  // Print container
  template <typename Stream, typename Container>
  void print
  (Stream& stream, Container const& container,
   bool delimit_first = 0, std::string delimiter = "\t")
  {
    // TODO: Choose this specialization when stream << container exists
    if constexpr (std::is_pod<Container>::value)
    {
      if (delimit_first)
        stream << delimiter;
      stream << container;
    }
    else
    {
      std::string delim = delimit_first ? delimiter : "";
      for (auto const& val : container)
      {
        stream << delim << val;
        delim = delimiter;
      }
    }
  }
  
  // Read contents of file as a sequence of doubles
  std::vector<double> read
  (std::string const& filename)
  {
    std::ifstream file{ filename };
    if (!file.is_open())
      throw useful::open_read_error(filename);
    std::vector<double> vals;
    double val;
    while(file >> val)
      vals.push_back(val);
    file.close();

    return vals;
  }
  
  // Functor that does nothing
  struct DoNothing
  {
    template <typename ...Args>
    void operator()(Args...) const
    {}
  };

  // Functor that always returns false
  struct DoFalse
  {
    template <typename ...Args>
    bool operator()(Args...) const
    { return false; }
  };
  
  // Functor that always returns true
  struct DoTrue
  {
    template <typename ...Args>
    bool operator()(Args...) const
    { return true; }
  };

  // Make an object passing a parameters object
  template <typename Object, typename Params>
  Object Create(Params const& parameters = useful::Empty())
  {
    Object tracker;
    tracker.Initialize(parameters);

    return tracker;
  }
  
  // Functor to make objects of given type,
  // passing arbitrary parameters
  template <typename Object>
  class Maker
  {
    template <typename ...Args>
    Object operator()(Args... args)
    { return Object{ args... }; }
  };

  // Functor to make an object on the heap,
  // passing arbitrary parameters
  template <typename T>
  struct Creator
  {
    using value_type = T;
    using pointer = T*;

    pointer operator() (void) const
    { return (new T ()); }

    pointer operator() (double param) const
    { return (new T(param)); }
  };

  // Functor to return back an argument without change
  template <typename Object_Type>
  struct Forward
  {
    Object_Type operator()(Object_Type const& object)
    { return object; }
  };

  // Functor to return const reference to an argument
  template <typename Object_Type>
  struct Forward_ref
  {
    Object_Type const& operator()(Object_Type const& object)
    { return object; }
  };

  // Types for selecting function implementations
  // at compile time
  template <typename TT> struct Selector_t{};
  template <typename TT, TT(val)> struct Selector{};

  // Ensure reference type
  template <typename T>
  T& ensure_ref(T const* obj)
  { return *obj; }

  // Ensure reference type
  template <typename T>
  T& ensure_ref(T const& obj)
  { return obj; }

  // Simple hash for a container by combining element hashes
  template <typename Container>
  struct hash_container
  {
    std::size_t operator()(Container const& container) const
    {
      using value_type = typename Container::value_type;
      std::size_t hash_val = std::hash<value_type>()(container[ 0 ]);
      for (size_t ii = 1; ii < container.size() - 1; ++ii)
      {
        hash_val = (hash_val ^ (std::hash<value_type>()(container[ii]) << 1)) >> 1;
      }
      hash_val = hash_val ^ (std::hash<value_type>()(container[container.size()-1]) << 1);

      return hash_val;
    }
  };

  // Combining seed with object hash
  template <typename T>
  void hash_combine(std::size_t& seed, T const& v)
  {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  // Hash for std::pair
  template <typename S, typename T>
  struct hash_pair
  {
    std::size_t operator()(std::pair<S,T > const& v) const
    {
      std::size_t seed = 0;
      hash_combine(seed, v.first);
      hash_combine(seed, v.second);
      return seed;
    }
  };

  // Replace NaNs
  template < typename T >
  void deNaN(T& to_replace, T replace_with = T())
  {
    if (to_replace != to_replace) to_replace = replace_with;
  }

  // Replace values below tolerance
  template <typename T>
  void chop(T& to_chop, T tolerance, T replace_with = T())
  {
    if (std::abs(to_chop - replace_with) < tolerance) to_chop = replace_with;
  }

  // Copy end element of container to position, then erase end element
  template <typename Container>
  void swap_erase(Container& container, std::size_t position)
  {
    if (position != container.size() - 1)
      container[position] = container.back();
    container.pop_back();
  }
  
  // Swap-erase all elements at positions in list
  template <typename Container, typename List>
  void swap_erase(Container& container, List& positions)
  {
    positions.sort(std::greater<std::size_t>{});
    for (auto const position : positions)
      swap_erase(container, position);
  }
  
  // Swap-erase all elements satisfying criterium
  template <typename Container, typename Criterium>
  void swap_erase_if(Container& container, Criterium criterium)
  {
    std::list<std::size_t> to_delete;
    for (std::size_t ii = 0; ii < container.size(); ++ii)
      if (criterium(container[ii]))
        to_delete.push_back(ii);
    
    swap_erase(container, to_delete);
  }

  // Copy position to end element of vector to position, then delete end element
  template <typename Container>
  void swap_delete(Container& container, std::size_t position)
  {
    if (position != container.size() - 1)
      container[position] = container.back();
    delete container.back();
    container.pop_back();
  }
  
  // Swap-delete all elements at positions in list
  template <typename Container, typename List>
  void swap_delete(Container& container, List& positions)
  {
    positions.sort(std::greater<std::size_t>{});
    for (auto const position : positions)
      swap_delete(container, position);
  }
  
  // Count lines in file
  std::size_t countlines(FILE *fin)
  {
    std::size_t lines = 0;
    int maxlength = 255;
    char buffer[ maxlength + 1 ];

    while (fgets(buffer , maxlength + 1 , fin) != NULL)
    {
      if (buffer[std::strlen(buffer) - 1] != '\n')
        throw "Line too long.";
      ++lines;
    }
    if (!feof(fin))
      throw "Could not reach end of file.";

    rewind(fin);

    return lines;
  }

  //	To check if a class has a method
  //	From here: https://stackoverflow.com/questions/29772601/why-is-sfinae-causing-failure-when-there-are-two-functions-with-different-signat
  //	bundle of types
  template <class...> struct types{using type=types;};
  template <class...> struct voider{using type=void;};
  //	Some dists still do not to have void_t
  template <class...Ts> using void_t=typename voider<Ts...>::type;
  //	hide the SFINAE stuff in a details namespace:
  namespace details
  {
    template <template <class...> class Z, class types, class=void>
    struct has_method : std::false_type {};
    template <template <class...> class Z, class...Ts>
    struct has_method<Z,types<Ts...>,void_t<Z<Ts...>>>:std::true_type{};
  }
  // has_method<template, types...> is true iff template <types...> is valid
  template <template <class...> class Z, class...Ts>
  using has_method=details::has_method<Z,types<Ts...>>;
  // Implemented here because some versions of gcc have trouble with the standard one
  template <typename T>
  bool isnan(T const& val)
  { return val != val; }
  
  // Apply function to each element of tuple
  template <typename Tuple, typename F, std::size_t ...Indices>
  void for_each_impl(Tuple const& tuple, F f, std::index_sequence<Indices...>)
  {
    using swallow = int[];
    (void)swallow{ 1, (f(std::get<Indices>(tuple)), void(), int{})... };
  }
  template <typename Tuple, typename F>
  void for_each(Tuple const& tuple, F f)
  {
    constexpr std::size_t N
      = std::tuple_size<std::remove_reference_t<Tuple>>::value;
    for_each_impl(tuple, f, std::make_index_sequence<N>{});
  }
  template <typename Tuple, typename F, std::size_t ...Indices>
  void for_each_impl
  (Tuple&& tuple, F&& f, std::index_sequence<Indices...>)
  {
    using swallow = int[];
    (void)swallow{ 1, (f(std::get<Indices>(std::forward<Tuple>(tuple))),
                       void(),
                       int{})... };
  }
  template <typename Tuple, typename F>
  void for_each(Tuple&& tuple, F&& f)
  {
    constexpr std::size_t N = std::tuple_size<std::remove_reference_t<Tuple>>::value;
    for_each_impl(std::forward<Tuple>(tuple),
                  std::forward<F>(f),
                  std::make_index_sequence<N>{});
  }

  // Indices for template metamagic
  template <std::size_t... Indices>
  struct indices
  { using next = indices<Indices..., sizeof...(Indices)>; };
  template <std::size_t size>
  struct build_indices
  { using type = typename build_indices<size-1>::type::next; };
  template <>
  struct build_indices< 0 >
  { using type = indices<>; };
}

#endif /* USEFUL_H_ */
