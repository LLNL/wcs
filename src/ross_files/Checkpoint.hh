#pragma once

#include <list>
#include <string>
#include <sstream>

// This could be done any number of ways faster/better than a list,
// but keep it for now.  We can optimize this later, so long as we
// have the following interface.
class CheckpointQueue {
 public:
   void push_back(std::string newString) { queue.push_back(newString); }
   void pop_back() { queue.pop_back(); }
   std::string back() { return queue.back(); }
   std::string front() { return queue.front(); }
   void pop_front() { queue.pop_front(); }

   std::list<std::string> queue;
};

class InputCheckpoint {
 public:
   std::string str() { return sss.str(); }
   
   template <typename TTT>
   InputCheckpoint& operator&(const TTT& ttt) {
      sss.write(reinterpret_cast<const char*>(&ttt), sizeof(TTT));
      return *this;
   }

   std::stringstream sss;
};

class OutputCheckpoint {
 public:
   OutputCheckpoint(std::string newString) {
      sss.str(newString);
   }
   template <typename TTT>
   OutputCheckpoint& operator&(TTT& ttt) {
      sss.read(reinterpret_cast<char *>(&ttt), sizeof(TTT));
      return *this;
   }

   std::stringstream sss;
};
