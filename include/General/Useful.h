/**
 \file General/Useful.h
 \author Tomás Aquino
 \date 03/15/2011
*/

// Miscelaneous collection of useful objects and algorithms

#ifndef GENERAL_USEFUL_H
#define GENERAL_USEFUL_H

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
  /** \brief Display execution time in human-readable format.
  \details Adapted from Howard Hinnant's answer here:
   https://stackoverflow.com/questions/22590821/convert-stdduration-to-human-readable-time */
  template <typename Clock>
  std::ostream& display_duration
  (std::ostream& stream,
   std::chrono::time_point<Clock> start_time,
   std::chrono::time_point<Clock> end_time)
  {
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
    using days = std::chrono::duration<int, std::ratio<86400>>;
    char fill = stream.fill();
    stream.fill('0');
    auto d = std::chrono::duration_cast<days>(ns);
    ns -= d;
    auto h = std::chrono::duration_cast<std::chrono::hours>(ns);
    ns -= h;
    auto m = std::chrono::duration_cast<std::chrono::minutes>(ns);
    ns -= m;
    auto s = std::chrono::duration_cast<std::chrono::seconds>(ns);
    ns -= s;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(ns);
    stream << std::setw(2) << d.count() << "d:"
           << std::setw(2) << h.count() << "h:"
           << std::setw(2) << m.count() << "m:"
           << std::setw(2) << s.count() << "s:"
           << std::setw(2) << ms.count() << "ms";
    stream.fill(fill);
    
    return stream;
  };
  
  /** \return Convert \p string to boolean. */
  bool stob(std::string const& string)
  {
    if (string == "true" || string == "1"
        || string == "True" || string == "TRUE")
      return true;
    if (string == "false" || string == "0"
        || string == "False" || string == "FALSE")
      return false;
    throw std::runtime_error{
      "Expected true or false, got"
      + string };
  }
  
  /**
   \param val Boolean value.
   \param boolalpha Use true/false if true or 1/0 if false to represent booleans.
   \return Convert boolean \p val to string. */
  std::string btos(bool val, bool boolalpha = true)
  {
    std::ostringstream stream;
    if (boolalpha)
      stream << std::boolalpha;
    else
      stream << std::noboolalpha;
    stream << val;
    
    return stream.str();
  }
  
  /** \brief Check if string is empty.
   \details
   Strings are considered empty if they hold:
   - Nothing.
   - ''
   - "" */
  bool empty(std::string const& str)
  {
    return str == ""
      || str == "''"
      || str == R"("")";
  }
  
  /** \brief Remove extension after last dot, including dot. */
  std::string remove_extension
  (std::string const& filename)
  {
    return filename.substr(0, filename.find_last_of('.'));
  }

  /** \brief Change extension after last dot.
   \note: \p new_extension should include dot if wanted. */
  std::string change_extension
  (std::string const& filename,
   std::string const& new_extension)
  {
    return remove_extension(filename) + new_extension;
  }

  /** \brief Expand environment variables.
   \details
   Based on Toby Speight's answer here:
   https://codereview.stackexchange.com/questions/172644/c-environment-variable-expansion . */
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
  
  /** \brief Compute widths of bins having given \p values as midpoints and given \p minimum (left edge).
   \details Widths are computed sequentially as twice the distance between the next midpoint and the current left edge. */
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
  
  /** \brief Compute widths of bins having given \p values as midpoints.
   \details Leftmost bin edge is the minimum value minus half the width of the first bin.
   
   Widths are computed sequentially as twice the distance between the next midpoint and the current left edge. */
  template <typename Value_type = double>
  std::vector<Value_type> get_bin_widths
  (std::vector<Value_type> const& values)
  {
    if (values.size() < 2)
      throw std::runtime_error{ "Less than two values, could not determine PDF bin widths" };
    double midpoint = (values[1]+values[0])/2.;
    return get_bin_widths(values, values[0]-(midpoint-values[0]));
  }
  
  /** \brief Check if sorted \p container contains \p val.
  \note: Container must be sorted. */
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
  
  /** \brief Check whether \p container contains \p value.
  \details
   Container must be sorted according to \p comp_less, equality is checked using \p comp_eq. */
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
  

  /** \brief Check if \p string ends with \p suffix. */
  bool endsWith(const std::string& string, const std::string& suffix)
  {
    return string.size() >= suffix.size() &&
      string.substr(string.size() - suffix.size()) == suffix;
  }
  
  /** \brief Split \p string into vector of strings at instances of \p delimeter.
  \details Empty entries can be included or discarded.
   
  Adapted from Beder Acosta Borges's answer here:
  https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c . */
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
  
  /** \return Exception for failing to parse a file line. */
  auto parse_error
  (std::string const& filename, std::string const& line)
  {
    return std::runtime_error{
      "Could not parse line\n" + line +
      "\nin file " + filename };
  }
  
  /** \return Exception for failing to parse a file. */
  auto parse_error_file(std::string const& filename)
  {
    return std::runtime_error{
      "Could not parse file " + filename };
  }
  
  /** \return Exception for failing to parse a line. */
  auto parse_error_line(std::string const& line)
  {
    return std::runtime_error{
      "Could not parse line\n" + line };
  }
  
  /** \return Exception for failing to open a file for reading. */
  auto open_read_error(std::string const& filename)
  {
    return std::runtime_error{
      "Could not open file " + filename + " for reading" };
  }
  
  /** \return Exception for failing to open a file for writing. */
  auto open_write_error(std::string const& filename)
  {
    return std::runtime_error{
      "Could not open file " + filename + " for writing" };
  }
  
  /** \return Exception for unexpected contents in file. */
  auto bad_file_contents(std::string const& filename)
  {
    return std::runtime_error{
      "Innapropriate contents in file " + filename };
  }
  
  /** \return TException for finding the end of a file before expected string. */
  auto bad_eof
  (std::string const& filename, std::string const& string)
  {
    return std::runtime_error{
    "Reached end of " +
      filename + " before " + string + " was found" };
  }
  
  /** \return Exception for bad parameters. */
  auto bad_parameters()
  {
    return std::invalid_argument{ "Inappropriate parameters" };
  }
  
  /** \return Exception for bad parameters and suggest help. */
  auto bad_parameters_help()
  {
    return std::invalid_argument{ "Inappropriate parameters"
      " (-h or --help for help)" };
  }
  
  /** \brief Check command line options for help flags. */
  bool check_options_help(int argc, const char* const argv[])
  {
    return argc == 2 &&
      (std::string(argv[1]) == "--help" ||
       std::string(argv[1]) == "-h");
  }
  
  /** \brief Open file for reading. */
  std::ifstream open_read(std::string const& filename)
  {
    std::ifstream file(filename);
    if (!file.is_open())
      throw useful::open_read_error(filename);
    
    return file;
  }
  
  /** \brief Open file for writing. */
  std::ofstream open_write(std::string const& filename)
  {
    std::ofstream file(filename);
    if (!file.is_open())
      throw useful::open_write_error(filename);
    
    return file;
  }
  
  /** \brief Load 1-column file into vector of doubles. */
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
  
  /** \brief Load 2-column file into pair of vectors of doubles. */
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
  
  /** \brief Load 3-column file,  first two columns into vector of pairs, last column into vector of doubles. */
  auto load_pair_1
  (std::string const& filename, std::size_t nr_estimate = 0,
   std::size_t header_lines = 0, std::string const& delim = " ")
  {
    using Value = double;
    using Container_pair = std::vector<std::pair<Value, Value>>;
    using Container_scalar = std::vector<Value>;
    using Output = std::tuple<Container_pair, Container_scalar>;
    
    Output output;
    std::get<0>(output).reserve(nr_estimate);
    std::get<1>(output).reserve(nr_estimate);
    
    auto file = open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line = split(line, delim);
      if (split_line.size() != 3)
        throw parse_error(filename, line);
      
      std::get<0>(output).push_back({ std::stod(split_line[0]), std::stod(split_line[1]) });
      std::get<1>(output).push_back(std::stod(split_line[2]));
    }
    file.close();
    
    return output;
  }
  
  /** \brief Load file into vector of vectors of doubles. */
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
  
  /** \brief Get number of numbers in first line of file. */
  std::size_t nr_numbers_in_first_line(std::string const& filename)
  {
    auto input = open_read(filename);
    std::string line;
    std::getline(input, line);
    std::stringstream stream(line);
    double number;
    std::size_t nr_numbers = 0;
    while (stream >> number)
      ++nr_numbers;
    
    return nr_numbers;
  }
  
  /** \brief Read next value from file. */
  template <typename Type>
  void read(std::ifstream& input, Type& val)
  {
    input >> val;
    if (input.fail())
      throw std::runtime_error{ "Could not read value" };
  }
  
  /** \class can_call_sqrt General/Useful.h "General/Useful.h"
   \brief To check whether class can call sqrt
   \details From Passer By's answer here:
   https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain . */
  template <typename = void, typename... Args>
  struct can_call_sqrt : std::false_type {};
  /** \class can_call_sqrt General/Useful.h "General/Useful.h"
  \brief To check whether class can call sqrt.
  \details From Passer By's answer here:
  https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain . */
  template <typename... Args>
  struct can_call_sqrt<
  std::void_t<decltype(std::sqrt(std::declval<Args>()...))>, Args...>
  : std::true_type {};
  /**
  \brief To check whether class can call sqrt.
  \details From Passer By's answer here:
  https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain . */
  template <typename... Args>
  inline constexpr bool can_call_sqrt_v = can_call_sqrt<void, Args...>::value;
  
  template <typename = void, typename... Args>
  struct can_call_abs : std::false_type {};
  template <typename... Args>
  struct can_call_abs<
  std::void_t<decltype(std::abs(std::declval<Args>()...))>, Args...>
  : std::true_type {};
  /**
  \brief To check whether class can call abs.
  \details From Passer By's answer here:
  https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain . */
  template <typename... Args>
  inline constexpr bool can_call_abs_v = can_call_abs<void, Args...>::value;
  
  
  /** Check whether classes can be used with binary operator
  \details
    From Richard Hodges's answer here:
    https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain . */
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
  
  /** \namespace notstd Operators not in std. */
  namespace notstd {

      /** \struct left_shift General/Useful.h "General/Useful.h"
       \brief Left shift operator. */
      struct left_shift {

          template <class L, class R>
          constexpr auto operator()(L&& l, R&& r) const
          noexcept(noexcept(std::forward<L>(l) << std::forward<R>(r)))
          -> decltype(std::forward<L>(l) << std::forward<R>(r))
          {
              return std::forward<L>(l) << std::forward<R>(r);
          }
      };

      /** \struct right_shift General/Useful.h "General/Useful.h"
      \brief Right shift operator. */
      struct right_shift {

          template <class L, class R>
          constexpr auto operator()(L&& l, R&& r) const
          noexcept(noexcept(std::forward<L>(l) >> std::forward<R>(r)))
          -> decltype(std::forward<L>(l) >> std::forward<R>(r))
          {
              return std::forward<L>(l) >> std::forward<R>(r);
          }
      };
    
      /** \struct operator_and General/Useful.h "General/Useful.h"
      \brief And operator. */
      struct operator_and {

          template <class L, class R>
          constexpr auto operator()(L&& l, R&& r) const
          noexcept(noexcept(std::forward<L>(l) & std::forward<R>(r)))
          -> decltype(std::forward<L>(l) & std::forward<R>(r))
          {
              return std::forward<L>(l) & std::forward<R>(r);
          }
      };

  }
  /** \brief Check whether classes can be used with equality operator. */
  template<class X, class Y> using has_equality = op_valid<X, Y, std::equal_to<>>;
  /** \brief Check whether classes can be used with inequality operator. */
  template<class X, class Y> using has_inequality = op_valid<X, Y, std::not_equal_to<>>;
  /** \brief Check whether classes can be used with less than operator. */
  template<class X, class Y> using has_less_than = op_valid<X, Y, std::less<>>;
  /** \brief Check whether classes can be used with less than or equal operator. */
  template<class X, class Y> using has_less_equal = op_valid<X, Y, std::less_equal<>>;\
  /** \brief Check whether classes can be used with greater than operator. */
  template<class X, class Y> using has_greater_than = op_valid<X, Y, std::greater<>>;
  /** \brief Check whether classes can be used with greater than or equal operator. */
  template<class X, class Y> using has_greater_equal = op_valid<X, Y, std::greater_equal<>>;
  /** \brief Check whether classes can be used with bit xor operator. */
  template<class X, class Y> using has_bit_xor = op_valid<X, Y, std::bit_xor<>>;
  /** \brief Check whether classes can be used with bit or operator. */
  template<class X, class Y> using has_bit_or = op_valid<X, Y, std::bit_or<>>;
  /** \brief Check whether classes can be used with left shift operator. */
  template<class X, class Y> using has_left_shift = op_valid<X, Y, notstd::left_shift>;
  /** \brief Check whether classes can be used with right shift operator. */
  template<class X, class Y> using has_right_shift = op_valid<X, Y, notstd::right_shift>;
  /** \brief Check whether classes can be used with plus operator. */
  template<class X, class Y> using has_plus = op_valid<X, Y, std::plus<>>;
  /** \brief Check whether classes can be used with minus operator. */
  template<class X, class Y> using has_minus = op_valid<X, Y, std::minus<>>;
  /** \brief Check whether classes can be used with multiplication operator. */
  template<class X, class Y> using has_multiplies = op_valid<X, Y, std::multiplies<>>;
  /** \brief Check whether classes can be used with division operator. */
  template<class X, class Y> using has_divides = op_valid<X, Y, std::divides<>>;
  /** \brief Check whether classes can be used with and operator. */
  template<class X, class Y> using has_and = op_valid<X, Y, notstd::operator_and>;
  
  /** \class has_value_type General/Useful.h "General/Useful.h"
   \brief Check if class defines the type value_type
   \details
   From here: https://gist.github.com/ilya-biryukov/887b7e543b72b49376ed . */
  template <class T>
  class has_value_type
  {
      struct One { char a[1]; };
      struct Two { char a[2]; };

      template <class U>
      static One foo(typename U::type*);

      template <class U>
      static Two foo(...);

  public:
      static const bool value = sizeof(foo<T>(nullptr)) == sizeof(One);
  };

  /** \return Sign of \p val. */
  template <typename T> int sgn(T val)
  { return (T(0) < val) - (val < T(0)); }

  /** \class StoreConst General/Useful.h "General/Useful.h"
   \brief Class holding const object. */
  template <typename Object_Type, typename Return_Type = Object_Type>
  struct StoreConst
  {
    using value_type = Return_Type;
    
    Return_Type operator()() const
    { return obj; }
    const Object_Type obj;
  };

  /** \class Store General/Useful.h "General/Useful.h"
  \brief Class holding object. */
  template <typename Object_Type, typename Return_Type = Object_Type>
  struct Store
  {
    using value_type = Return_Type;
    
    Return_Type operator()() const
    { return obj; }
    Object_Type obj;
  };

  /** \class Empty General/Useful.h "General/Useful.h"
  \brief Class holding nothing and doing nothing. */
  struct Empty
  {
    template <typename ...Args>
    Empty(Args...){}
    Empty(){}
  };
  
  /** \brief Print values in container. */
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
  
  /** \brief Read contents of file as a sequence of doubles. */
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
  
  /** \struct DoNothing General/Useful.h "General/Useful.h"
   \brief Functor that does nothing. */
  struct DoNothing
  {
    template <typename ...Args>
    void operator()(Args...) const
    {}
  };

  /** \struct DoFalse General/Useful.h "General/Useful.h"
  \brief Functor that always returns false. */
  struct DoFalse
  {
    template <typename ...Args>
    bool operator()(Args...) const
    { return false; }
  };
  
  /** \struct DoTrue General/Useful.h "General/Useful.h"
  \brief Functor that always returns true. */
  struct DoTrue
  {
    template <typename ...Args>
    bool operator()(Args...) const
    { return true; }
  };

  /** \brief Make an object and initialize it given \p parameters.
   \note Object must implement:
   - initialize(Params const&). */
  template <typename Object, typename Params>
  Object create(Params const& parameters = useful::Empty())
  {
    Object object;
    object.initialize(parameters);

    return object;
  }
  
  /** \class Maker General/Useful.h "General/Useful.h"
   \brief Functor to make objects of given type, passing arbitrary parameters. */
  template <typename Object>
  struct Maker
  {
    template <typename ...Args>
    Object operator()(Args... args) const
    { return Object{ args... }; }
  };

  /** \class Creator General/Useful.h "General/Useful.h"
   \brief Functor to make an object on the heap,  passing arbitrary parameters. */
  template <typename T>
  struct Creator
  {
    using value_type = T;
    using pointer = T*;

    template <typename ...Args>
    pointer operator()(Args... args) const
    { return new T{ args... }; }
  };

  /** \class Forward General/Useful.h "General/Useful.h"
   \brief Functor to forward an argument. */
  template <typename Object_Type>
  struct Forward
  {
    Object_Type operator()(Object_Type const& object)
    { return object; }
  };

  /** \class Forward_ref General/Useful.h "General/Useful.h"
   \brief Functor to forward a const reference. */
  template <typename Object_Type>
  struct Forward_ref
  {
    Object_Type const& operator()(Object_Type const& object)
    { return object; }
  };

  /** \class Selector_t General/Useful.h "General/Useful.h"
   \brief Type for selecting function implementations at compile time. */
  template <typename TT> struct Selector_t{};
  
  /** \class Selector General/Useful.h "General/Useful.h"
   \brief Type for selecting function implementations at compile time. */
  template <typename TT, TT(val)> struct Selector{};

  /** \brief Used along with ensure_ref(T const& obj) to ensure \p T is a reference type. */
  template <typename T>
  T& ensure_ref(T const* obj)
  { return *obj; }

  /** \brief Used along with ensure_ref(T const* obj) to ensure \p T is a reference type. */
  template <typename T>
  T& ensure_ref(T const& obj)
  { return obj; }

  /** \class hash_container General/Useful.h "General/Useful.h"
   \brief Simple hash for a container by combining element hashes. */
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

  /** \brief Combine \p seed with the hash of an object \p v. */
  template <typename T>
  void hash_combine(std::size_t& seed, T const& v)
  {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  /** \class hash_pair General/Useful.h "General/Useful.h"
   \brief Hash for std::pair. */
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

  /** \brief If \p to_replace is NaN, replace it with \p replace_with. */
  template <typename T>
  void deNaN(T& to_replace, T replace_with = T())
  {
    if (to_replace != to_replace) to_replace = replace_with;
  }

  /** \brief If \p to_chop is smaller than \p tolerance, replace it with \p replace_with. */
  template <typename T>
  void chop(T& to_chop, T tolerance, T replace_with = T())
  {
    if (std::abs(to_chop - replace_with) < tolerance) to_chop = replace_with;
  }

  /** \brief Copy end element of \p container to position \p position, then erase end element. */
  template <typename Container>
  void swap_erase(Container& container, std::size_t position)
  {
    if (position != container.size() - 1)
      container[position] = container.back();
    container.pop_back();
  }
  
  /** \brief Swap-erase elements in \p container at positions \p positions.
   \note \p positions will be sorted, and must implement \c sort(std::greater<std::size_t>).*/
  template <typename Container, typename List>
  void swap_erase(Container& container, List&& positions)
  {
    positions.sort(std::greater<std::size_t>{});
    for (auto const position : positions)
      swap_erase(container, position);
  }
  
  /** \brief Swap-erase elements of \p container satisfying \p criterion. */
  template <typename Container, typename Criterion>
  void swap_erase_if(Container& container, Criterion&& criterion)
  {
    for (std::size_t ii = container.size(); ii --> 0;)
      if (criterion(container[ii]))
        swap_erase(container, ii);
  }

  /** \brief Like swap_erase(Container&, std::size_t), but call \c delete on the removed element. */
  template <typename Container>
  void swap_delete(Container& container, std::size_t position)
  {
    if (position != container.size() - 1)
      container[position] = container.back();
    delete container.back();
    container.pop_back();
  }
  
  /** \brief Swap-delete elements in \p container at positions \p positions. */
  template <typename Container, typename List>
  void swap_delete(Container& container, List& positions)
  {
    positions.sort(std::greater<std::size_t>{});
    for (auto const position : positions)
      swap_delete(container, position);
  }
              
  /** \brief Swap-delete elements of \p container satisfying \p criterion. */
  template <typename Container, typename Criterion>
  void swap_delete_if(Container& container, Criterion&& criterion)
  {
    for (std::size_t ii = container.size(); ii --> 0;)
      if (criterion(container[ii]))
        swap_delete(container, ii);
  }

  /** \class types General/Useful.h "General/Useful.h"
   \brief Bundle of types. */
  template <class...> struct types{using type=types;};
  /** \class voider General/Useful.h "General/Useful.h"
   \brief Metafunction to turn a bundle of type arguments into \c void. */
  template <class...> struct voider{using type=void;};
  /** \brief void_t General/Useful.h "General/Useful.h"
   \brief For dists that do not to have \c void_t. */
  template <class...Ts> using void_t = typename voider<Ts...>::type;
  /** \namespace useful::details Implementation details. */
  namespace details
  {
    template <template <class...> class Z, class types, class=void>
    struct has_method : std::false_type {};
    template <template <class...> class Z, class...Ts>
    struct has_method<Z,types<Ts...>,void_t<Z<Ts...>>>:std::true_type{};
  }
  /** \brief types General/Useful.h "General/Useful.h"
  \brief To check if a class has a method.
  \details has_method<template, types...> is true iff template <types...> is valid.
   
   From here: https://stackoverflow.com/questions/29772601/why-is-sfinae-causing-failure-when-there-are-two-functions-with-different-signat . */
  template <template <class...> class Z, class...Ts>
  using has_method=details::has_method<Z, types<Ts...>>;
  
  /** \brief begin_result General/Useful.h "General/Useful.h"
  \brief Helper to check if a class has begin() method. */
  template<class X>
  using begin_result = decltype(std::declval<X>().begin());
  /** \brief has_begin General/Useful.h "General/Useful.h"
  \brief Check if a class has begin() method. */
  template<class X>
  using has_begin = has_method<begin_result, X>;
  
  /** \brief Implemented here because some versions of gcc have trouble with std::isnan. */
  template <typename T>
  bool isnan(T const& val)
  { return val != val; }
  
  /** \brief Helper to apply function to each element of tuple. */
  template <typename Tuple, typename F, std::size_t ...Indices>
  void for_each_impl(Tuple const& tuple, F f, std::index_sequence<Indices...>)
  {
    using swallow = int[];
    (void)swallow{ 1, (f(std::get<Indices>(tuple)), void(), int{})... };
  }
  
  /** \brief Apply function to each element of tuple. */
  template <typename Tuple, typename F>
  void for_each(Tuple const& tuple, F f)
  {
    constexpr std::size_t N
      = std::tuple_size<std::remove_reference_t<Tuple>>::value;
    for_each_impl(tuple, f, std::make_index_sequence<N>{});
  }
  
  /** \brief Helper to apply function to each element of tuple. */
  template <typename Tuple, typename F, std::size_t ...Indices>
  void for_each_impl
  (Tuple&& tuple, F&& f, std::index_sequence<Indices...>)
  {
    using swallow = int[];
    (void)swallow{ 1, (f(std::get<Indices>(std::forward<Tuple>(tuple))),
                       void(),
                       int{})... };
  }
  
  /** \brief Apply function to each element of tuple. */
  template <typename Tuple, typename F>
  void for_each(Tuple&& tuple, F&& f)
  {
    constexpr std::size_t N = std::tuple_size<std::remove_reference_t<Tuple>>::value;
    for_each_impl(std::forward<Tuple>(tuple),
                  std::forward<F>(f),
                  std::make_index_sequence<N>{});
  }

  /** \class indices General/Useful.h "General/Useful.h"
   \brief Indices for template metamagic. */
  template <std::size_t... Indices>
  struct indices
  { using next = indices<Indices..., sizeof...(Indices)>; };

  /** \class build_indices General/Useful.h "General/Useful.h"
  \brief Indices for template metamagic. */
  template <std::size_t size>
  struct build_indices
  { using type = typename build_indices<size-1>::type::next; };
  template <>
  struct build_indices<0>
  { using type = indices<>; };

  /** \class Parameters_or_empty General/Useful.h "General/Useful.h"
   \brief  Get Parameters type for true or Empty type for false. */
  template <bool, typename>
  struct Parameters_or_empty{ using type = useful::Empty; };
  template <typename T>
  struct Parameters_or_empty<true, T>
  { using type = typename T::Parameters; };
  template <bool condition, typename T>
  using Parameters_or_empty_t = typename Parameters_or_empty<condition, T>::type;
  
  /** \class has_time_step_setter General/Useful.h "General/Useful.h"
   \brief Check if type T has member function void time_step(double).
   \details Adapted from kispaljr's answer here: https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function .
  */
  template <typename T> struct has_time_step_setter
  {
      typedef char (&Yes)[1];
      typedef char (&No)[2];

      template<class U>
      static Yes test(U* data,
                      typename std::enable_if<std::is_void<
                        decltype(data->time_step(0.))>::value>::type* = 0);
      static No test(...);
      static const bool value =
        sizeof(Yes) ==
          sizeof(has_time_step_setter::test((typename std::remove_reference<T>::type*)0));
  };
  
  /** \class has_time_step_getter General/Useful.h "General/Useful.h"
   \brief Check if type T has member function double time_step().
   \details Adapted from kispaljr's answer here: https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function .
  */
  template <typename T> struct has_time_step_getter
  {
      typedef char (&Yes)[1];
      typedef char (&No)[2];

      template<class U>
      static Yes test(U* data,
                      typename std::enable_if<std::is_same<
                        double,
                        decltype(data->time_step())>::value>::type* = 0);
      static No test(...);
      static const bool value = sizeof(Yes) == sizeof(has_time_step_getter::test((typename std::remove_reference<T>::type*)0));
  };
  
  /** \class has_time_step General/Useful.h "General/Useful.h"
   \brief Check if type T has member variable double time_step.
   \details Adapted from kispaljr's answer here: https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function .
  */
  template <typename T> struct has_time_step
  {
      typedef char (&Yes)[1];
      typedef char (&No)[2];

      template<class U>
      static Yes test(U* data,
                      typename std::enable_if<std::is_same<
                        double,
                        decltype(data->time_step)>::value>::type* = 0);
      static No test(...);
      static const bool value = sizeof(Yes) == sizeof(has_time_step::test((typename std::remove_reference<T>::type*)0));
  };
}

#endif /* GENERAL_USEFUL_H */
