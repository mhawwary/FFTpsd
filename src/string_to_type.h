#ifndef __string_to_type_h__
#define __string_to_type_h__

// C++ includes
#include <string>

  /**
   * Takes the string \p s and returns the matching
   * enumeration of type \p T.
   */
  template <typename T>
  T string_to_enum (const std::string& s);

  /**
   * Takes the enumeration \p e of type \p T
   * and returns the matching string.
   */
  template <typename T>
  std::string enum_to_string (const T e);
  
#endif
