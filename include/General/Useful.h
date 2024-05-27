/**
 \file General/Useful.h
 \author Tomás Aquino
 \date 03/15/2011
 \brief Miscelaneous collection of useful objects and algorithms
*/

#ifndef GENERAL_USEFUL_H
#define GENERAL_USEFUL_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <list>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/** \namespace useful Miscellaneous utilities. */
namespace useful {
// Adapted from Howard Hinnant's answer here:
// https://stackoverflow.com/questions/22590821/convert-stdduration-to-human-readable-time
/** \brief Display execution time in human-readable format.  */
inline std::ostream &display_duration(std::ostream &stream,
                                      std::chrono::nanoseconds ns) {
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
  stream << std::setw(2) << d.count() << "d:" << std::setw(2) << h.count()
         << "h:" << std::setw(2) << m.count() << "m:" << std::setw(2)
         << s.count() << "s:" << std::setw(3) << ms.count() << "ms";
  stream.fill(fill);

  return stream;
};

/** \brief Display execution time in human-readable format.  */
template <typename Clock>
std::ostream &display_duration(std::ostream &stream,
                               std::chrono::time_point<Clock> start_time,
                               std::chrono::time_point<Clock> end_time) {
  return display_duration(stream,
                          std::chrono::duration_cast<std::chrono::nanoseconds>(
                              end_time - start_time));
}

/** \return Convert \p string to boolean. */
inline bool stob(std::string const &string) {
  if (string == "true" || string == "1" || string == "True" || string == "TRUE")
    return true;
  if (string == "false" || string == "0" || string == "False" ||
      string == "FALSE")
    return false;
  throw std::runtime_error{std::string("Expected true or false, got") + string};
}

/**
 \param val Boolean value.
 \param boolalpha Use true/false if true or 1/0 if false to represent booleans.
 \return Convert boolean \p val to string. */
inline std::string btos(bool val, bool boolalpha = true) {
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
inline bool is_empty(std::string const &str) {
  return str == "" || str == "''" || str == R"("")";
}

/** \brief Remove extension after last dot, including dot. */
inline std::string remove_extension(std::string const &filename) {
  return filename.substr(0, filename.find_last_of('.'));
}

/** \brief Change extension after last dot.
 \note: \p new_extension should include dot if wanted. */
inline std::string change_extension(std::string const &filename,
                                    std::string const &new_extension) {
  return remove_extension(filename) + new_extension;
}

/** \brief Remove \verbatim \r \endverbatim character */
inline std::string remove_carriage_return(std::string const &str) {
  std::string str_copy;
  if (!str_copy.empty() && str_copy.back() == '\r')
    str_copy.pop_back();
  return str_copy;
}

/** \brief Remove \verbatim \r \endverbatim character */
inline std::string &remove_carriage_return_in_place(std::string &str) {
  if (!str.empty() && str.back() == '\r')
    str.pop_back();
  return str;
}

// Based on Toby Speight's answer here:
// https://codereview.stackexchange.com/questions/172644/c-environment-variable-expansion
/** \brief Expand environment variables. */
inline std::string expand_env(std::string text) {
  static const std::regex env_re{R"--(\$\{([^}]+)\})--"};
  std::smatch match;
  while (std::regex_search(text, match, env_re)) {
    auto const from = match[0];
    text.replace(from.first, from.second, std::getenv(match[1].str().c_str()));
  }
  return text;
}

/** \brief Check if sorted \p container contains \p val.
\note: Container must be sorted. */
template <typename T, typename U>
bool contains(const std::vector<T> &container, U const &val) {
  auto it =
      std::lower_bound(container.begin(), container.end(), val,
                       [](T const &elem, U const &val) { return elem < val; });

  return it != container.end() && *it == val;
}

/** \brief Check whether \p container contains \p value.
\details
 Container must be sorted according to \p comp_less, equality is checked using
\p comp_eq. */
template <typename T, typename U, typename Comp_less, typename Comp_eq>
bool contains(std::vector<T> const &container, U const &val,
              Comp_less comp_less = std::less<void>{},
              Comp_eq comp_eq = std::equal_to<void>{}) {
  auto it =
      std::lower_bound(container.begin(), container.end(), val, comp_less);

  return it != container.end() && comp_eq(*it, val);
}

/** \brief Check if \p string ends with \p suffix. */
inline bool ends_with(const std::string &string, const std::string &suffix) {
  return string.size() >= suffix.size() &&
         string.substr(string.size() - suffix.size()) == suffix;
}

/** \brief Split \p string into vector of strings at instances of \p
 delimeter. \tparam empty_entries Whether to keep empty entries between
 delimiting characters */
template <bool empty_entries = false>
std::vector<std::string> split(std::string const &str,
                               std::string const &delims = "\t,|\r ") {
  std::vector<std::string> tokens;
  std::size_t pos_last_non_delim = str.find_last_not_of(delims);
  std::size_t upper_limit = pos_last_non_delim != std::string::npos
                                ? pos_last_non_delim + 1
                                : str.length();
  for (std::size_t start = str.find_first_not_of(delims), end;
       start < upper_limit; start = end + 1) {
    std::size_t position = str.find_first_of(delims, start);
    end = position != std::string::npos ? position : str.length();

    std::string token = str.substr(start, end - start);
    if (!token.empty())
      tokens.push_back(token);
  }

  return tokens;
}

template <>
inline std::vector<std::string> split<true>(std::string const &str,
                                            std::string const &delims) {
  std::vector<std::string> tokens;
  for (std::size_t start = 0, end; start < str.length(); start = end + 1) {
    std::size_t position = str.find_first_of(delims, start);
    end = position != std::string::npos ? position : str.length();
    std::string token = str.substr(start, end - start);
    tokens.push_back(token);
  }
  if (str.empty() ||
      std::any_of(delims.begin(), delims.end(),
                  [&str](auto const &delim) { return delim == str.back(); }))
    tokens.push_back("");

  return tokens;
}

/** \return Exception for failing to parse a file line. */
inline auto parse_error(std::string const &filename, std::string const &line) {
  return std::runtime_error{std::string("Could not parse line\n") + line +
                            "\nin file " + filename};
}

/** \return Exception for failing to parse a file. */
inline auto parse_error_file(std::string const &filename) {
  return std::runtime_error{std::string("Could not parse file ") + filename};
}

/** \return Exception for failing to parse a line. */
inline auto parse_error_line(std::string const &line) {
  return std::runtime_error{std::string("Could not parse line\n") + line};
}

/** \return Exception for failing to open a file for reading. */
inline auto open_read_error(std::string const &filename) {
  return std::runtime_error{std::string("Could not open file ") + filename +
                            " for reading"};
}

/** \return Exception for failing to open a file for writing. */
inline auto open_write_error(std::string const &filename) {
  return std::runtime_error{std::string("Could not open file ") + filename +
                            " for writing"};
}

/** \return Exception for unexpected contents in file. */
inline auto bad_file_contents(std::string const &filename) {
  return std::runtime_error{std::string("Innapropriate contents in file ") +
                            filename};
}

/** \return TException for finding the end of a file before expected string. */
inline auto bad_eof(std::string const &filename, std::string const &string) {
  return std::runtime_error{std::string("Reached end of ") + filename +
                            " before " + string + " was found"};
}

/** \return Exception for bad parameters. */
inline auto bad_parameters() {
  return std::invalid_argument{"Inappropriate parameters"};
}

/** \return Exception for bad parameters and suggest help. */
inline auto bad_parameters_help() {
  return std::invalid_argument{"Inappropriate parameters"
                               " (-h or --help for help)"};
}

/** \brief Check command line options for help flags. */
inline bool check_options_help(int argc, const char *const argv[]) {
  return argc == 2 &&
         (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h");
}

/** \brief Open file for reading. */
inline std::ifstream open_read(std::string const &filename) {
  std::ifstream file(filename);
  if (!file.is_open())
    throw useful::open_read_error(filename);

  return file;
}

/** \brief Open file for writing. */
inline std::ofstream open_write(std::string const &filename) {
  std::ofstream file(filename);
  if (!file.is_open())
    throw useful::open_write_error(filename);

  return file;
}

/** \brief Load 1-column file into vector of doubles. */
inline auto load_1(std::string const &filename, std::size_t nr_estimate = 0,
            std::size_t header_lines = 0,
            std::string const &delims = "\t,|\r ") {
  using Value = double;
  using Container = std::vector<Value>;
  Container values;
  values.reserve(nr_estimate);

  auto file = open_read(filename);
  std::string line;
  for (std::size_t ll = 0; ll < header_lines; ++ll)
    std::getline(file, line);

  while (std::getline(file, line)) {
    remove_carriage_return_in_place(line);
    std::vector<std::string> split_line = split(line, delims);
    if (split_line.size() != 1)
      throw parse_error(filename, line);
    values.push_back(std::stod(split_line[0]));
  }
  file.close();

  return values;
}

/** \brief Load 2-column file into pair of vectors of doubles. */
inline auto load_2(std::string const &filename, std::size_t nr_estimate = 0,
            std::size_t header_lines = 0,
            std::string const &delims = "\t,|\r ") {
  using Value = double;
  using Container = std::vector<Value>;
  std::pair<Container, Container> values;
  values.first.reserve(nr_estimate);
  values.second.reserve(nr_estimate);

  auto file = open_read(filename);
  std::string line;
  for (std::size_t ll = 0; ll < header_lines; ++ll)
    std::getline(file, line);

  while (std::getline(file, line)) {
    remove_carriage_return_in_place(line);
    std::vector<std::string> split_line = split(line, delims);
    if (split_line.size() != 2)
      throw parse_error(filename, line);

    values.first.push_back(std::stod(split_line[0]));
    values.second.push_back(std::stod(split_line[1]));
  }
  file.close();

  return values;
}

/** \brief Load 3-column file,  first two columns into vector of pairs, last
 * column into vector of doubles. */
inline auto load_pair_1(std::string const &filename, std::size_t nr_estimate = 0,
                 std::size_t header_lines = 0,
                 std::string const &delims = "\t,|\r ") {
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
    std::getline(file, line);

  while (std::getline(file, line)) {
    remove_carriage_return_in_place(line);
    std::vector<std::string> split_line = split(line, delims);
    if (split_line.size() != 3)
      throw parse_error(filename, line);

    std::get<0>(output).push_back(
        {std::stod(split_line[0]), std::stod(split_line[1])});
    std::get<1>(output).push_back(std::stod(split_line[2]));
  }
  file.close();

  return output;
}

/** \brief Load file into vector of vectors of doubles. */
inline auto load(std::string const &filename, std::size_t nr_columns,
          std::size_t nr_estimate = 0, std::size_t header_lines = 0,
          std::string const &delims = "\t,|\r ") {
  using Value = double;
  using Container = std::vector<Value>;
  std::vector<Container> values(nr_columns);
  for (auto &val : values)
    val.reserve(nr_estimate);

  auto file = open_read(filename);
  std::string line;
  for (std::size_t ll = 0; ll < header_lines; ++ll)
    std::getline(file, line);

  while (std::getline(file, line)) {
    remove_carriage_return_in_place(line);
    std::vector<std::string> split_line = split(line, delims);
    if (split_line.size() != nr_columns)
      throw parse_error(filename, line);

    for (std::size_t cc = 0; cc < nr_columns; ++cc)
      values[cc].push_back(std::stod(split_line[cc]));
  }
  file.close();

  return values;
}

/** \brief Get number of numbers in first line of file. */
inline std::size_t nr_numbers_in_first_line(std::string const &filename) {
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

/** \brief Convert string to type.
 \note \c Type must be default-constructible. */
template <typename Type> Type convert_to(std::string const &str) {
  std::istringstream ss(str);
  Type val;
  ss >> val;
  return val;
}

/** \brief Remove comments after escape sequence */
inline std::string clear_escape(std::string const &str,
                         std::string const &escape_sequence) {
  return str.substr(0, str.find(escape_sequence));
}

/** \brief Remove comments after escape sequence */
inline std::string &clear_escape_in_place(std::string &str,
                                   std::string const &escape_sequence) {
  std::size_t pos = str.find(escape_sequence);
  if (pos != std::string::npos)
    str.erase(pos);
  return str;
}

/** \brief Read next value from file. */
template <typename Type> void read(std::ifstream &input, Type &val) {
  input >> val;
  if (input.fail())
    throw std::runtime_error{"Could not read value"};
}

/** \brief Extract first value from line, discarding escaped lines and the rest
 * of the line */
template <typename Type>
void read_first_from_line(std::ifstream &input, Type &val,
                          std::string const &escape_sequence,
                          std::string const &delims = "\t,|\r ") {
  std::string line;
  while (std::getline(input, line))
    if (line.find(escape_sequence) != 0)
      break;
  remove_carriage_return_in_place(line);
  clear_escape_in_place(line, escape_sequence);
  auto split_line = split(line, delims);
  if (split_line.empty())
    throw std::runtime_error{"Could not read value"};

  val = convert_to<Type>(split_line[0]);
}

/** \brief Extract first value from line, discarding escaped lines and the rest
 of the line \note Type must be default-constructible. */
template <typename Type>
Type read_first_from_line(std::ifstream &input,
                          std::string const &escape_sequence,
                          std::string const &delims = "\t,|\r ") {
  Type val;
  read_first_from_line(input, val, escape_sequence, delims);
  return val;
}

/** \return Sign of \p val. */
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

/** \class StoreConst General/Useful.h "General/Useful.h"
 \brief Class holding const object. */
template <typename Object_Type, typename Return_Type = Object_Type>
struct StoreConst {
  using value_type = Return_Type;

  Return_Type operator()() const { return obj; }
  const Object_Type obj;
};

/** \class Store General/Useful.h "General/Useful.h"
\brief Class holding object. */
template <typename Object_Type, typename Return_Type = Object_Type>
struct Store {
  using value_type = Return_Type;

  Return_Type operator()() const { return obj; }
  Object_Type obj;
};

/** \class Empty General/Useful.h "General/Useful.h"
\brief Class holding nothing and doing nothing. */
struct Empty {
  template <typename... Args> Empty(Args...) {}
  Empty() {}
};

/** \brief Print values in container. */
template <typename Stream, typename Container>
void print(Stream &stream, Container const &container, bool delimit_first = 0,
           std::string delimiter = "\t") {
  // TODO: Choose this specialization when stream << container exists
  if constexpr (std::is_pod<Container>::value) {
    if (delimit_first)
      stream << delimiter;
    stream << container;
  } else {
    std::string delim = delimit_first ? delimiter : "";
    for (auto const &val : container) {
      stream << delim << val;
      delim = delimiter;
    }
  }
}

/** \brief Print values in container. */
template <typename Stream, typename Container>
void print(Stream &stream, Container const &container, int width,
           std::ios_base &(*alignment)(std::ios_base &) = std::right) {
  if constexpr (std::is_pod<Container>::value)
    stream << alignment << std::setw(width) << container;
  else
    for (auto const &val : container)
      stream << alignment << std::setw(width) << val;
}

/** \brief Read contents of file as a sequence of doubles. */
inline std::vector<double> read(std::string const &filename) {
  std::ifstream file{filename};
  if (!file.is_open())
    throw useful::open_read_error(filename);
  std::vector<double> vals;
  double val;
  while (file >> val)
    vals.push_back(val);
  file.close();

  return vals;
}

/** \struct DoNothing General/Useful.h "General/Useful.h"
 \brief Functor that does nothing. */
struct DoNothing {
  template <typename... Args> void operator()(Args...) const {}
};

/** \struct DoFalse General/Useful.h "General/Useful.h"
\brief Functor that always returns false. */
struct DoFalse {
  template <typename... Args> bool operator()(Args...) const { return false; }
};

/** \struct DoTrue General/Useful.h "General/Useful.h"
\brief Functor that always returns true. */
struct DoTrue {
  template <typename... Args> bool operator()(Args...) const { return true; }
};

/** \brief Make an object and initialize it given \p parameters.
 \note Object must implement:
 - initialize(Params const&). */
template <typename Object, typename Params>
Object create(Params const &parameters = useful::Empty()) {
  Object object;
  object.initialize(parameters);

  return object;
}

/** \class Maker General/Useful.h "General/Useful.h"
 \brief Functor to make objects of given type, passing arbitrary parameters. */
template <typename Object> struct Maker {
  template <typename... Args> Object operator()(Args... args) const {
    return Object{args...};
  }
};

/** \class Creator General/Useful.h "General/Useful.h"
 \brief Functor to make an object on the heap,  passing arbitrary parameters. */
template <typename T> struct Creator {
  using value_type = T;
  using pointer = T *;

  template <typename... Args> pointer operator()(Args... args) const {
    return new T{args...};
  }
};

/** \class Forward General/Useful.h "General/Useful.h"
 \brief Functor to forward an argument. */
template <typename Object_Type> struct Forward {
  Object_Type operator()(Object_Type const &object) { return object; }
};

/** \class Forward_ref General/Useful.h "General/Useful.h"
 \brief Functor to forward a const reference. */
template <typename Object_Type> struct Forward_ref {
  Object_Type const &operator()(Object_Type const &object) { return object; }
};

/** \class hash_container General/Useful.h "General/Useful.h"
 \brief Simple hash for a container by combining element hashes. */
template <typename Container> struct hash_container {
  std::size_t operator()(Container const &container) const {
    using value_type = typename Container::value_type;
    std::size_t hash_val = std::hash<value_type>()(container[0]);
    for (size_t ii = 1; ii < container.size() - 1; ++ii) {
      hash_val =
          (hash_val ^ (std::hash<value_type>()(container[ii]) << 1)) >> 1;
    }
    hash_val = hash_val ^
               (std::hash<value_type>()(container[container.size() - 1]) << 1);

    return hash_val;
  }
};

/** \brief Combine \p seed with the hash of an object \p v. */
template <typename T> void hash_combine(std::size_t &seed, T const &v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/** \class hash_pair General/Useful.h "General/Useful.h"
 \brief Hash for std::pair. */
template <typename S, typename T> struct hash_pair {
  std::size_t operator()(std::pair<S, T> const &v) const {
    std::size_t seed = 0;
    hash_combine(seed, v.first);
    hash_combine(seed, v.second);
    return seed;
  }
};

/** \brief If \p to_replace is NaN, replace it with \p replace_with. */
template <typename T> void deNaN(T &to_replace, T replace_with = T()) {
  if (to_replace != to_replace)
    to_replace = replace_with;
}

/** \brief If \p to_chop is smaller than \p tolerance, replace it with \p
 * replace_with. */
template <typename T> void chop(T &to_chop, T tolerance, T replace_with = T()) {
  if (std::abs(to_chop - replace_with) < tolerance)
    to_chop = replace_with;
}

/** \brief Copy end element of \p container to position \p position, then erase
 * end element. */
template <typename Container>
void swap_erase(Container &container, std::size_t position) {
  if (position != container.size() - 1)
    container[position] = container.back();
  container.pop_back();
}

/** \brief Swap-erase elements in \p container at positions \p positions.
 \note \p positions will be sorted, and must implement \c
 sort(std::greater<std::size_t>).*/
template <typename Container, typename List>
void swap_erase(Container &container, List &&positions) {
  positions.sort(std::greater<std::size_t>{});
  for (auto const position : positions)
    swap_erase(container, position);
}

/** \brief Swap-erase elements of \p container satisfying \p criterion. */
template <typename Container, typename Criterion>
void swap_erase_if(Container &container, Criterion &&criterion) {
  for (std::size_t ii = container.size(); ii-- > 0;)
    if (criterion(container[ii]))
      swap_erase(container, ii);
}

/** \brief Like swap_erase(Container&, std::size_t), but call \c delete on the
 * removed element. */
template <typename Container>
void swap_delete(Container &container, std::size_t position) {
  if (position != container.size() - 1)
    container[position] = container.back();
  delete container.back();
  container.pop_back();
}

/** \brief Swap-delete elements in \p container at positions \p positions. */
template <typename Container, typename List>
void swap_delete(Container &container, List &positions) {
  positions.sort(std::greater<std::size_t>{});
  for (auto const position : positions)
    swap_delete(container, position);
}

/** \brief Swap-delete elements of \p container satisfying \p criterion. */
template <typename Container, typename Criterion>
void swap_delete_if(Container &container, Criterion &&criterion) {
  for (std::size_t ii = container.size(); ii-- > 0;)
    if (criterion(container[ii]))
      swap_delete(container, ii);
}

/** \brief Get value of shell environment variable \c name (do not pass \$ in \c
 * name). */
inline std::string getenv(std::string const &name) {
  char const *val = std::getenv(name.c_str());
  return val == nullptr ? std::string{} : std::string{val};
}

/** \brief Get $HOME environment variable. */
inline std::string getenv_home() { return getenv("HOME"); }

/** \brief Get $HOME environment variable, throw if empty */
inline std::string getenv_home_check() {
  std::string home = getenv_home();
  if (home.empty())
    throw std::runtime_error{"$HOME environmental variable is undefined"};
  return home;
}

/** \brief Expand ~ if present at the beginning of directory \c dir. */
inline std::string expand_home_dir(std::string dir) {
  if (dir[0] == '~')
    dir.replace(0, 1, getenv_home_check());
  return dir;
}

/** \brief Expand ~ if present at the beginning of directory \c dir. */
inline std::string &expand_home_dir_inplace(std::string &dir) {
  if (dir[0] == '~')
    dir.replace(0, 1, getenv_home_check());
  return dir;
}
} // namespace useful

#endif /* GENERAL_USEFUL_H */
