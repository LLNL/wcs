WCS coding style mostly follows the standard C++ libraries and Boost projects.

* overriding style is lowercase separated with underbar
* member fields: m_*
* function names: lowercase with underbar
* class names: start with uppercase (camel case after that e.g. DataType)
* templated types: uppercase
  * derived typedef types:  using value_type = SomeType::value_type_t
* header preprocessor guard: __<NAMESPACE_PATH_NAME>_HPP__
* 2 space, no tabs
* comments:
  * doxygen:  
    * /// single line comment 
    * /** multi-line comment */
    * @TODO  - TODO note
  * inside of a function use //
  * outside use a doxygen comment
    * minimize blocks of // or /* */ comments that are not picked up by doxygen
* for variable names, spell out the full dictionary words as much as possible
  and refrain from using undocumented or unfamiliar acronyms
* use const keyword when a variable is not supposed to be modified and a member
  function does not update any member state.
* avoid relying on implicit type promotion or conversion,
  and use explicit casting such as static_cast
* do not convert an unsigned type to a signed type unless it is algorithmically 
  necessary
* check and fix compiler warnings
* avoid using "using namespace xxx", and explicitly specify the
  namespace of each symbol of a different namespace referenced
* do not put implementation details inside of a class definition body
