// hash.h

#ifndef HASH_HEADER
#define HASH_HEADER

#include "fileread.hh"

// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;

class hash_table
// Simple hash table for pairs of ints
//  Nfind   "index" number stored with Nstore 
//  Nstore  number to be used to retrieve Nfind (must be > 0)
//
//  Typically Nstore is a randomly large-range positive integer
//  (eg batch number) and Nfind is a running index (eg 0->n-1)
//  hash_table.add(Nstore, Nfind)    adds a pair to the list
//  hash_table.lookup(Nstore)        retrieves the index number
//                             belonging to Nstore
//
//
// Code translated from Fortran functions ccp4_hash_setup/lookup/zeroit
// This is a simple-minded hash function suitable for batch-number hashing
// but not necessarily for other purposes
//

//  hash_table.add(Nstore, Nfind) sets up a value so that when
//  hash_table.lookup(Nstore) is later evaluated it will return Nfind. This
//  function will allow the efficient retrieval of an identifier for a
//  large range variable (such as a batch number). The values of the
//  function hash_table.lookup(Nstore) are stored in the arrays Nstore_list,
//  Nfind_list of length table_size, which is the prime number used to
//  generate the function.
//  
//  NOTES: A hash table is a way of storing information so that it easily
//  be retrieved without the need for indexing or long searches. Nstore is
//  referred to as the "key", which is "hashed" (computer- science speak
//  for "messed up") by the hashing function (in this case (Nstore %
//  table_size) to determine where the value pair will be stored. The
//  function lookup can then search on the same basis when supplied with
//  the key, to retreive the pair in (at most) 3 calculations. Note that
//  the table size MUST BE A PRIME in order for this method to work.
//  


{
public:

// Constructor:-
//   hash_table()                       dummy, must be followed by a
//                                      call to hash_table.set_size
//   hash_table(const int& size)        the table will be the next biggest prime
//                                      larger than size
  hash_table() {table_size = 0;}  
  hash_table(const int& size);

  void set_size(const int& size);
  // get size
  int Size() const {return table_size;}

  // Add pair to list
  void add(const int& Nstore, const int& Nfind);

  // Return index Nstore for stored number
  //   -1 if not found
  int lookup(const int& Nstore) const;
  // Returns actual stored number for index
  // = -1 if out of range
  int number(const int& index) const;

  // Returns smallest prime <= argument
  static int prime(const int& min_value);

  void Clear() {setup();}

  //! format table for save/restore
  std::string FormatSave() const;
  //! restore
  void Restore(Fileread& FR);

private:
  static const int MinPrime;  // this hash method doesn't work well
                        //  on small numbers so give a minimum value

  int table_size;
  std::vector<int> Nstore_list;
  std::vector<int> Nfind_list;

  void setup();
};

#endif
