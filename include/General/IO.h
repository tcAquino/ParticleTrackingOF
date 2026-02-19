/**
 \file General/IO.h
 \author Tomas Aquino
 \date 29/03/2024
 \brief Utilities for input and output.
*/

#ifndef GENERAL_IO_H
#define GENERAL_IO_H

#include "General/Meta.h"
#include "General/Useful.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/** \namespace Utilities for input and output. */
namespace io {
/**
   \class StreamScopeFormat General/IO.h "General/IO.h"
   \brief Construct from stream object to restrict formatting to local scope.
*/
// From here:
// https://github.com/gelldur/gcpp/blob/master/src/gcpp/stream/StreamScopeFormat.hpp
template <typename T> struct StreamScopeFormat {
  explicit StreamScopeFormat(std::basic_ios<T> &stream)
      : init{nullptr}, stream{stream} {
    init.copyfmt(stream);
  }

  ~StreamScopeFormat() {
    stream.copyfmt(init); // restore format
  }

  std::ios init;
  std::basic_ios<T> &stream;
};

/** \brief Remove extension after last dot, including dot. */
inline std::string remove_extension(std::string const &filename) {
  return filename.substr(0, filename.find_last_of('.'));
}

/**
   \brief Change extension after last dot.
   \note: \p new_extension should include dot if wanted.
*/
inline std::string change_extension(std::string const &filename,
                                    std::string const &new_extension) {
  return remove_extension(filename) + new_extension;
}

/** \brief Remove \verbatim \r \endverbatim character */
inline std::string &remove_carriage_return_in_place(std::string &str) {
  if (!str.empty() && str.back() == '\r') {
    str.pop_back();
  }
  return str;
}

/** \brief Remove \verbatim \r \endverbatim character */
inline std::string remove_carriage_return(std::string str) {
  remove_carriage_return_in_place(str);
  return str;
}

/**
   \brief Get value of shell environment variable \c name.
   \note Do not pass \$ in \c name.
*/
template <bool check = true> std::string getenv(std::string const &name) {
  char const *env = std::getenv(name.c_str());
  return env == nullptr ? std::string{} : std::string{env};
}

/**
   \brief Specialization of \c getenv() that throws if empty.
*/
template <> std::string getenv<true>(std::string const &name) {
  std::string env = getenv<false>(name);
  if (env.empty()) {
    throw std::runtime_error{"$HOME environmental variable is undefined"};
  }
  return env;
}

/** \brief Get $HOME environment variable. */
template <bool check = true> std::string getenv_home() {
  return getenv<check>("HOME");
}

/** \brief Expand ~ if present at the beginning of directory \c dir. */
template <bool check = true>
std::string &expand_home_dir_in_place(std::string &dir) {
  if (dir[0] == '~') {
    dir.replace(0, 1, getenv_home<check>());
  }
  return dir;
}

/** \brief Expand ~ if present at the beginning of directory \c dir. */
template <bool check = true> std::string expand_home_dir(std::string dir) {
  expand_home_dir_in_place<check>(dir);
  return dir;
}

// Based on Toby Speight's answer here:
// https://codereview.stackexchange.com/questions/172644/c-environment-variable-expansion
/** \brief Expand environment variables. */
template <bool check = true>
std::string &expand_env_in_place(std::string &str) {
  static const std::regex env_re{R"--(\$\{([^}]+)\})--"};
  std::smatch match;
  while (std::regex_search(str, match, env_re)) {
    auto const from = match[0];
    str.replace(from.first, from.second, getenv<check>(match[1].str()));
  }
  return str;
}

/** \brief Expand environment variables. */
template <bool check = true> std::string expand_env(std::string str) {
  expand_env_in_place<true>(str);
  return str;
}

/** \brief Remove comments after escape sequence. */
inline std::string clear_escape(std::string const &str,
                                std::string const &escape_sequence = "#") {
  return str.substr(0, str.find(escape_sequence));
}

/** \brief Remove comments after escape sequence. */
inline std::string &
clear_escape_in_place(std::string &str,
                      std::string const &escape_sequence = "#") {
  std::size_t pos = str.find(escape_sequence);
  if (pos != std::string::npos) {
    str.erase(pos);
  }
  return str;
}

/**
   \brief Check if string is empty.
   \details Strings are considered empty if they hold:
   - Nothing.
   - ''
   - ""
*/
inline bool is_empty(std::string const &str) {
  return str == "" || str == "''" || str == R"("")";
}

/** \return Ordinal string (1st, 2nd, 3rd, 4th, ...) for unsigned integer. */
inline std::string ordinal(std::size_t number) {
  switch (number) {
  case 1: {
    return "1st";
  }
  case 2: {
    return "2nd";
  }
  case 3: {
    return "3rd";
  }
  default: {
    return std::to_string(number) + "th";
  }
  }
}

/** \return Convert \p string to boolean. */
inline bool stob(std::string const &string) {
  if (string == "true" || string == "1" || string == "True" ||
      string == "TRUE") {
    return true;
  }
  if (string == "false" || string == "0" || string == "False" ||
      string == "FALSE") {
    return false;
  }
  throw std::runtime_error{std::string("Expected true or false, got") + string};
}

/**
   \param val Boolean value.
   \param boolalpha Use true/false if true or 1/0 if false to represent
   booleans.
   \return Convert boolean \p val to string.
*/
inline std::string btos(bool val, bool boolalpha = true) {
  std::ostringstream stream;
  if (boolalpha) {
    stream << std::boolalpha;
  } else {
    stream << std::noboolalpha;
  }
  stream << val;
  return stream.str();
}

/**
   \brief Convert string to type.
*/
template <typename Type> Type &stotype(std::string const &str, Type &val) {
  std::istringstream ss(str);
  ss >> val;
  return val;
}

/**
   \brief Convert string to type.
*/
template <typename Type> Type stotype(std::string const &str) {
  if constexpr (std::is_same_v<Type, double>) {
    return std::stod(str);
  } else if constexpr (std::is_same_v<Type, int>) {
    return std::stoi(str);
  } else if constexpr (std::is_same_v<Type, bool>) {
    return stob(str);
  } else if constexpr (std::is_same_v<Type, long>) {
    return std::stol(str);
  } else if constexpr (std::is_same_v<Type, long long>) {
    return std::stoll(str);
  } else if constexpr (std::is_same_v<Type, float>) {
    return std::stof(str);
  } else if constexpr (std::is_same_v<Type, long double>) {
    return std::stold(str);
  } else if constexpr (std::is_same_v<Type, long int>) {
    return std::stol(str);
  } else if constexpr (std::is_same_v<Type, unsigned long>) {
    return std::stoul(str);
  } else if constexpr (std::is_same_v<Type, unsigned long long>) {
    return std::stoull(str);
  } else if constexpr (std::is_same_v<Type, unsigned>) {
    return static_cast<unsigned>(std::stoul(str));
  } else if constexpr (std::is_same_v<Type, std::string>) {
    return str;
  } else if constexpr (std::is_same_v<Type, const char *>) {
    return str.c_str();
  } else {
    Type val;
    stotype(str, val);
    return val;
  }
}

/**
   \brief Split \p string into vector of strings at instances of \p delimeter
   and append to \c tokens.
   \tparam empty_entries Whether to keep empty entries between delimiting
   characters.
*/
template <bool empty_entries = false>
std::vector<std::string> &split(std::string const &str,
                                std::vector<std::string> &tokens,
                                std::string const &delims = "\t,|\r ") {
  std::size_t pos_last_non_delim = str.find_last_not_of(delims);
  std::size_t upper_limit = pos_last_non_delim != std::string::npos
                                ? pos_last_non_delim + 1
                                : str.length();
  for (std::size_t start = str.find_first_not_of(delims), end;
       start < upper_limit; start = end + 1) {
    std::size_t position = str.find_first_of(delims, start);
    end = position != std::string::npos ? position : str.length();

    std::string token = str.substr(start, end - start);
    if (!token.empty()) {
      tokens.push_back(token);
    }
  }
  return tokens;
}

template <>
inline std::vector<std::string> &split<true>(std::string const &str,
                                             std::vector<std::string> &tokens,
                                             std::string const &delims) {
  for (std::size_t start = 0, end; start < str.length(); start = end + 1) {
    std::size_t position = str.find_first_of(delims, start);
    end = position != std::string::npos ? position : str.length();
    std::string token = str.substr(start, end - start);
    tokens.push_back(token);
  }
  if (str.empty() ||
      std::any_of(delims.begin(), delims.end(),
                  [&str](auto const &delim) { return delim == str.back(); })) {
    tokens.push_back("");
  }
  return tokens;
}

/**
   \brief Split \p string into vector of strings at instances of \p delimeter.
   \tparam empty_entries Whether to keep empty entries between delimiting
   characters.
*/
template <bool empty_entries = false>
std::vector<std::string> split(std::string const &str,
                               std::string const &delims = "\t,|\r ") {
  std::vector<std::string> tokens;
  split<empty_entries>(str, tokens, delims);
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

/** \return Exception for finding the end of a file before expected string. */
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

/** \brief Open file for reading. */
inline std::ifstream open_read(std::string const &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw open_read_error(filename);
  }
  return file;
}

/** \brief Open file for writing. */
inline std::ofstream
open_write(std::string const &filename,
           std::ios_base::openmode mode = std::ios_base::out) {
  std::ofstream file(filename, mode);
  if (!file.is_open()) {
    throw open_write_error(filename);
  }
  return file;
}

/** \brief Check command line options for help flags. */
inline bool check_options_help(int argc, const char *const argv[]) {
  return argc >= 2 &&
         (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h");
}

/**
   \brief Extract first value from line, discarding escaped lines and the rest
   of the line.
*/
template <typename Type>
Type &read_first_from_line(Type &val, std::ifstream &input,
                           std::string const &escape_sequence = "#",
                           std::string const &delims = "\t,|\r ") {
  std::string line;
  while (std::getline(input, line)) {
    if (line.find(escape_sequence) != 0) {
      break;
    }
  }
  clear_escape_in_place(remove_carriage_return_in_place(line), escape_sequence);

  auto split_line = split(line, delims);
  if (split_line.empty()) {
    throw std::runtime_error{"Could not read value"};
  }
  val = stotype<Type>(split_line[0]);
  return val;
}

/**
   \brief Extract first value from line, discarding escaped lines and the rest
   of the line.
   \note Type must be default-constructible.
*/
template <typename Type>
Type read_first_from_line(std::ifstream &input,
                          std::string const &escape_sequence = "#",
                          std::string const &delims = "\t,|\r ") {
  Type val;
  read_first_from_line(val, input, escape_sequence, delims);
  return val;
}

/** \brief Read next value from file. */
template <typename Type> Type &read(std::ifstream &input, Type &val) {
  input >> val;
  if (input.fail()) {
    throw std::runtime_error{"Could not read value"};
  }
  return val;
}

/** \brief Read contents of file as a sequence of values of given type. */
template <typename Type = double>
std::vector<Type> read(std::string const &filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    throw open_read_error(filename);
  }

  std::vector<Type> vals;
  Type val;
  while (file >> val) {
    vals.push_back(val);
  }
  return vals;
}

/**
   \class tuple_types General/IO.h "General/IO.h"
   \brief Helper object to read tuples.
*/
template <typename... Args> struct tuple_types {
  tuple_types(std::tuple<Args...>) {}

  static auto read_tuple(std::vector<std::string> const &strings,
                         std::size_t &index) {
    return std::tuple{[&]() { return stotype<Args>(strings[index++]); }()...};
  }
};

/**
   \brief Extract value from string.
   \param string String to extract from.
   \param error Base error message in case of failure.
   \return Converted value.
*/
template <typename Type>
Type read(std::string const &string, std::string const &error) {
  try {
    return stotype<Type>(string);
  } catch (const std::exception &except) {
    throw std::runtime_error{error.empty() ? except.what()
                                           : error + " : " + except.what()};
  }
}

/**
   \brief Extract value from string or use default.
   \param string String to extract from.
   \param error Base error message in case of failure.
   \param result Variable to place result.
   \return Reference to result.
*/
template <typename Type>
Type &read(std::string const &string, std::string const &error, Type &result) {
  try {
    stotype<Type>(string, result);
  } catch (const std::exception &except) {
    throw std::runtime_error{error.empty() ? except.what()
                                           : error + " : " + except.what()};
  }
}

/**
   \brief Extract value from string or use default.
   \param string String to extract from.
   \param default_value Default value, used if \c string is empty according to
   io::is_empty().
   \param error Base error message in case of failure.
   \return Converted value or default.
*/
template <typename Type>
Type read_or_default(std::string const &string, Type default_value,
                     std::string const &error) {
  if (io::is_empty(string)) {
    return default_value;
  } else {
    return read<Type>(string, error);
  }
}

/**
   \brief Extract value from string or use default.
   \param string String to extract from.
   \param default_value Default value, used if \c string is empty according to
   io::is_empty().
   \param error Base error message in case of failure.
   \param result Variable to place result.
   \return Reference to result.
*/
template <typename Type>
Type &read_or_default(std::string const &string, Type default_value,
                      std::string const &error, Type &result) {
  if (io::is_empty(string)) {
    result = default_value;
  } else {
    read<Type>(string, result, error);
  }
  return result;
}

/**
   \brief Extract value from container of strings and increment \p index.
   \param strings Strings to extract from.
   \param index Index of first string to read; is incremented by one per value
   read.
   \param error Base error message in case of failure.
   \return Extracted value.
*/
template <typename Type>
Type read(std::vector<std::string> const &strings, std::size_t &index,
          std::string const &error) {
  try {
    if constexpr (meta::is_instance_of_v<Type, std::pair>) {
      useful::check_size(strings, index, 2);
      std::size_t start = index;
      index += 2;
      return {stotype<typename Type::first_type>(strings[start]),
              stotype<typename Type::second_type>(strings[start + 1])};
    } else if constexpr (meta::is_instance_of_v<Type, std::tuple>) {
      useful::check_size(strings, index, std::tuple_size_v<Type>);
      return decltype(tuple_types(std::declval<Type>()))::read_tuple(strings,
                                                                     index);
    } else {
      useful::check_size(strings, index, 1);
      return stotype<Type>(strings[index++]);
    }
  } catch (const std::exception &except) {
    throw std::runtime_error{error.empty() ? except.what()
                                           : error + " : " + except.what()};
  }
}

/**
   \brief Extract value from container of strings.
   \param strings Strings to extract from.
   \param index Index of first string to read; is incremented by one per value
   read.
   \param nr_to_read Number of elements to read.
   \param error Base error message in case of failure.
   \return Vector of extracted values.
*/
template <typename Type>
std::vector<Type> read(std::vector<std::string> const &strings,
                       std::size_t &index, std::size_t nr_to_read,
                       std::string const &error) {
  std::vector<Type> result;
  result.reserve(nr_to_read);
  for (std::size_t ii = 0; ii < nr_to_read; ++ii) {
    result.push_back(read<Type>(strings, index, error));
  }
  return result;
}

/**
 \brief Extract value from container of strings and increment \p index.
 \param strings Strings to extract from.
 \param index Index of first string to read; is incremented by one per value
 read.
 \param error Base error message in case of failure.
 \param result Variable to extract into, converting based on type.
 \return Reference to extracted value.
*/
template <typename Type>
Type &read(std::vector<std::string> const &strings, std::size_t &index,
           std::string const &error, Type &result) {
  try {
    if constexpr (meta::is_instance_of_v<Type, std::pair>) {
      useful::check_size(strings, index, 2);
      std::size_t start = index;
      index += 2;
      result = {stotype<typename Type::first_type>(strings[start]),
                stotype<typename Type::second_type>(strings[start + 1])};
    } else if constexpr (meta::is_instance_of_v<Type, std::tuple>) {
      check_size(strings, index, std::tuple_size_v<Type>);
      result = decltype(tuple_types(std::declval<Type>()))::read_tuple(strings,
                                                                       index);
    } else {
      useful::check_size(strings, index, 1);
      stotype<Type>(strings[index++], result);
    }
  } catch (const std::exception &except) {
    throw std::runtime_error{error.empty() ? except.what()
                                           : error + " : " + except.what()};
  }

  return result;
}

/**
   \brief Extract value from container of strings or use default.
   \param strings Strings to extract from.
   \param index Index of first string to read; is incremented by one per value
   read unless there are no more strings.
   \param error Base error message in case of failure.
   \param default_value Default value, used if not enough strings or read value
   is empty according to io::is_empty().
   \return Extracted value.
*/
template <typename Type>
Type read_or_default(std::vector<std::string> const &strings,
                     std::size_t &index, Type default_value,
                     std::string const &error) {
  std::string value = index < strings.size()
                          ? io::read<std::string>(strings, index, error)
                          : "";
  if (io::is_empty(value)) {
    return default_value;
  } else {
    return read<Type>(strings, --index, error);
  }
}

/**
 \brief Extract value from container of strings and increment \p index.
 \param strings Strings to extract from.
 \param index Index of first string to read; is incremented by one per value
 read unless there are no more strings.
 \param error Base error message in case of failure.
 \param result Variable to extract into, converting based on type.
 \param default_value Default value, used if not enough strings or read value
 is empty according to io::is_empty().
 \return Reference to extracted value.
*/
template <typename Type>
Type &read_or_default(std::vector<std::string> const &strings,
                      std::size_t &index, Type default_value,
                      std::string const &error, Type &result) {
  std::string value = index < strings.size()
                          ? io::read<std::string>(strings, index, error)
                          : "";
  if (io::is_empty(value)) {
    result = default_value;
  } else {
    read(strings, --index, error, result);
  }
  return result;
}

/**
   \brief Extract values from container of strings and increment \p index.
   \param strings Strings to extract from.
   \param index Index of first string to read; is incremented by one per value
   read.
   \param error Base error message in case of failure.
   \param result Vector of variables to extract into, converting based on type;
   should be passed with the intended size.
   \return Reference to vector of extracted values.
*/
template <typename Type>
std::vector<Type> &read(std::vector<std::string> const &strings,
                        std::size_t &index, std::string const &error,
                        std::vector<Type> &result) {
  for (auto &val : result) {
    read(strings, index, error, val);
  }
  return result;
}

/**
   \brief Extract values from container of strings and increment \p index.
   \param strings Strings to extract from.
   \param index Index of first string to read; is incremented by one per value
   read.
   \param nr_to_read Number of elements to read.
   \param error Base error message in case of failure.
   \param result Vector of variables to append to, converting based on type.
   \return Reference to vector of extracted values.
*/
template <typename Type>
std::vector<Type> &read(std::vector<std::string> const &strings,
                        std::size_t &index, std::size_t nr_to_read,
                        std::string const &error, std::vector<Type> &result) {
  result.reserve(result.size() + nr_to_read);
  for (std::size_t ii = 0; ii < nr_to_read; ++ii) {
    result.push_back(read<Type>(strings, index, error));
  }
  return result;
}

/**
   \brief Extract values from container of strings and increment \p index.
   \param strings Strings to extract from.
   \param index Index of first string to read; is incremented by one per value
   read.
   \param error Base error message in case of failure.
   \param results variables to extract into, converting based on type.
*/
template <typename... Types>
void read(std::vector<std::string> const &strings, std::size_t &index,
          std::string const &error, Types &...results) {
  ([&]() { read(strings, index, error, results); }(), ...);
}

/**
   \brief Read and split next line, discarding escaped lines.
   \tparam empty_entries Whether to keep empty entries between delimiting
   characters.
*/
template <bool empty_entries = false>
auto split_line(std::ifstream &input, std::string const &escape_sequence = "#",
                std::string const &delims = "\t,|\r ") {
  std::string line;
  while (std::getline(input, line)) {
    if (line.find(escape_sequence) != 0) {
      break;
    }
  }
  clear_escape_in_place(remove_carriage_return_in_place(line), escape_sequence);
  return split<empty_entries>(line, delims);
}

/**
   \brief Read and split next line, discarding escaped lines, appending result
   to to \c tokens. \tparam empty_entries Whether to keep empty entries between
   delimiting characters. \return \c true if an unescaped line was read, \c
   false otherwise. */
template <bool empty_entries = false>
bool split_line(std::ifstream &input, std::vector<std::string> &tokens,
                std::string const &escape_sequence = "#",
                std::string const &delims = "\t,|\r ") {
  std::string line;
  bool read = false;
  while (std::getline(input, line)) {
    if (line.find(escape_sequence) != 0) {
      read = true;
      break;
    }
  }
  if (!read) {
    return false;
  }
  clear_escape_in_place(remove_carriage_return_in_place(line), escape_sequence);
  split<empty_entries>(line, tokens, delims);
  return true;
}

/** \brief Get number of numbers in first line of file. */
inline std::size_t nr_numbers_in_first_line(std::string const &filename) {
  auto input = open_read(filename);
  std::string line;
  std::getline(input, line);
  std::stringstream stream(line);
  double number;
  std::size_t nr_numbers = 0;
  while (stream >> number) {
    ++nr_numbers;
  }
  return nr_numbers;
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
  for (std::size_t ll = 0; ll < header_lines; ++ll) {
    std::getline(file, line);
  }

  while (std::getline(file, line)) {
    remove_carriage_return_in_place(line);
    std::vector<std::string> split_line = split(line, delims);
    if (split_line.size() != 1) {
      throw parse_error(filename, line);
    }
    values.push_back(std::stod(split_line[0]));
  }

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
  for (std::size_t ll = 0; ll < header_lines; ++ll) {
    std::getline(file, line);
  }

  while (std::getline(file, line)) {
    remove_carriage_return_in_place(line);
    std::vector<std::string> split_line = split(line, delims);
    if (split_line.size() != 2) {
      throw parse_error(filename, line);
    }
    values.first.push_back(std::stod(split_line[0]));
    values.second.push_back(std::stod(split_line[1]));
  }

  return values;
}

/**
   \brief Load 3-column file, first two columns into vector of pairs, last
   column into vector of doubles.
*/
inline auto load_pair_1(std::string const &filename,
                        std::size_t nr_estimate = 0,
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
    if (split_line.size() != 3) {
      throw parse_error(filename, line);
    }

    std::get<0>(output).push_back(
        {std::stod(split_line[0]), std::stod(split_line[1])});
    std::get<1>(output).push_back(std::stod(split_line[2]));
  }

  return output;
}

/** \brief Load file into vector of vectors of doubles. */
inline auto load(std::string const &filename, std::size_t nr_columns,
                 std::size_t nr_estimate = 0, std::size_t header_lines = 0,
                 std::string const &delims = "\t,|\r ") {
  using Value = double;
  using Container = std::vector<Value>;
  std::vector<Container> values(nr_columns);
  for (auto &val : values) {
    val.reserve(nr_estimate);
  }

  auto file = open_read(filename);
  std::string line;
  for (std::size_t ll = 0; ll < header_lines; ++ll) {
    std::getline(file, line);
  }

  while (std::getline(file, line)) {
    remove_carriage_return_in_place(line);
    std::vector<std::string> split_line = split(line, delims);
    if (split_line.size() != nr_columns) {
      throw parse_error(filename, line);
    }

    for (std::size_t cc = 0; cc < nr_columns; ++cc) {
      values[cc].push_back(std::stod(split_line[cc]));
    }
  }

  return values;
}

/**
   \brief Print values in container.
   \note Pass true or false, not 0 or 1, in \c delimit_first to avoid choosing a
   different overload.
*/
template <typename Container>
std::ostream &print(std::ostream &stream, Container const &container,
                    bool delimit_first = false, std::string delimiter = "\t") {
  if constexpr (std::is_scalar_v<Container>) {
    if (delimit_first) {
      stream << delimiter;
    }
    stream << container;
  } else {
    std::string delim = delimit_first ? delimiter : "";
    for (auto const &val : container) {
      stream << delim << val;
      delim = delimiter;
    }
  }
  return stream;
}

/** \brief Print values in container. */
template <typename Container>
std::ostream &print(std::ostream &stream, Container const &container, int width,
                    std::ios_base &(*alignment)(std::ios_base &) = std::right) {
  io::StreamScopeFormat guard{stream};
  stream << alignment;
  if constexpr (std::is_scalar_v<Container>) {
    stream << std::setw(width) << container;
  } else {
    for (auto const &val : container) {
      stream << std::setw(width) << val;
    }
  }
  return stream;
}

/**
   \struct Logger General/IO.h "General/IO.h"
   \brief Base object to handle log messages.
*/
struct Logger {
  virtual void nonewline(std::string const &message) = 0;
  virtual void operator()(std::string const &) = 0;
};

/**
   \struct StreamLogger General/IO.h "General/IO.h"
   \brief Logger object that uses provided stream.
*/
template <typename Stream = std::ostream &>
struct StreamLogger : public Logger {
  StreamLogger(Stream &&stream) : _stream{std::forward<Stream>(stream)} {}

  void nonewline(std::string const &message) override { _stream << message; }

  void operator()(std::string const &message) override {
    nonewline(message);
    _stream << std::endl;
  }

protected:
  Stream _stream;
};
template <typename Stream> StreamLogger(Stream &&) -> StreamLogger<Stream>;

/**
   \struct NullLogger General/IO.h "General/IO.h"
   \brief Logger object that ignores messages.
*/
struct NullLogger : public Logger {
  void nonewline(std::string const &) override{};
  void operator()(std::string const &) override{};
};

/**
   \struct FileLogger General/IO.h "General/IO.h"
   \brief Logger object that logs to a file.
*/
struct FileLogger : public StreamLogger<std::ofstream> {
  FileLogger(std::string const &filename,
             std::ios_base::openmode mode = std::ios_base::out)
      : StreamLogger{open_write(filename, mode)} {}

protected:
  using StreamLogger::_stream;
};

/**
   \struct CoutLogger General/IO.h "General/IO.h"
   \brief Logger object that outputs to standard output.
*/
struct CoutLogger : public StreamLogger<> {
  CoutLogger() : StreamLogger{std::cout} {}
};

/** \brief Line of hyphens. */
inline std::string line(std::size_t nr_characters = 80, bool endl = true) {
  return std::string(nr_characters, '-') + (endl ? "\n" : "");
}
} // namespace io

#endif /* GENERAL_IO_H */
