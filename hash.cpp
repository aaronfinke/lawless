// hash.cpp

// Has table & prime number things


#include "hash.hh"
#include "util.hh"

const int hash_table::MinPrime = 1009;

  //--------------------------------------------------------------
int hash_table::prime(const int& min_value)
  //  Returns smallest prime >= argument 
{ 
  // Force odd and >=  MinPrime
  int prime_number = Max((min_value/2)*2+1, MinPrime);

  for( ; ; )
    { 
      bool prime_found = true;
      int max_times = ( prime_number / 2 );
      for ( int i = Min(3, max_times-1); i < max_times; ++i )
	{
	  if ( ( ( prime_number / i ) * i ) == prime_number )
	    {
	      prime_found = false;
	      break;
	    }
	}
      if ( prime_found ) 
	break;
      prime_number += 2;
    }
  
  return (prime_number);
}
//--------------------------------------------------------------
void hash_table::setup()
{
  Nstore_list.clear();
  Nstore_list.reserve(table_size);
  Nstore_list.assign(table_size, -1);
  Nfind_list.clear();
  Nfind_list.reserve(table_size);
  Nfind_list.assign(table_size, -1);
}
//--------------------------------------------------------------
hash_table::hash_table (const int& size) 
  // Set up hash table with size = the smallest prime >= size
{
  table_size = prime(size);
  setup();
} 
//--------------------------------------------------------------
void hash_table::set_size(const int& size) 
  // Set size for hash table with size = the smallest prime >= size
{
  table_size = prime(size);
  setup();
} 
//--------------------------------------------------------------
// Add pair to list
void hash_table::add(const int& Nstore, const int& Nfind)
{
  if (table_size <= 0)
    {Message::message(Message_fatal( "hash_table::add - table not initialised"));}

  if (Nstore <= 0)
    {Message::message(Message_fatal( "hash_table::add - Nstore non-positive"));}

  int n = Nstore;
  int index;
  while (true)
    {
      index = n % table_size;
      if ((n-Nstore) >= 3*table_size)
	{Message::message(Message_fatal( "hash_table::setup - overflowed hash table"));}
      if (Nstore_list[index] < 0) break;
      n += 3;
    }
  Nstore_list[index] = Nstore;
  Nfind_list[index] = Nfind;
  return;
}
//--------------------------------------------------------------
int hash_table::lookup(const int& Nstore) const
  // Return index Nstore for stored number
  //   -1 if not found
{
  int n = Nstore;
  int index;
  int count = 0;

  while (count++ < table_size)
    {
      index = n % table_size;
      if (Nstore == Nstore_list[index])
	return Nfind_list[index];
      n += 3;
    }
  return -1;
}
//--------------------------------------------------------------
// Returns actual stored number for index
// = -1 if out of range
int hash_table::number(const int& index) const 
{
  if (index < 0 || index > table_size) {return -1;}
  for (int i=1;i<table_size;i++) {
    if (index == Nfind_list[i])
	return Nstore_list[i];
  }
  return -1;
}
//--------------------------------------------------------------
std::string hash_table::FormatSave() const
// Save data for restore
{
  std::string dump = "HashTable V1 {\n";
  dump += "TableSize "+ clipper::String(table_size)+"\n";
  // Don't dump empty slots
  for (int i=0;i<table_size;++i) {
    if (Nstore_list[i] >= 0) {
      dump += "HashIndexStoreFind "+clipper::String(i)+" "+
	clipper::String(Nstore_list[i])+
	" "+clipper::String(Nfind_list[i]);
      dump += "\n";
    }
  }
  return dump+"}\n";
}
//--------------------------------------------------------------
void hash_table::Restore(Fileread& FR) 
// Restore from file
{
  FR.ReadTag("HashTable"); // fails if tag does not match
  if (FR.GetTag() != "V1") {  // version check
    clipper::Message::message(Message_fatal
      ("hash_table::Restore incompatible version in "+FR.Filename()));
  }
  FR.Skip();
  FR.ReadTag("TableSize"); table_size = FR.Int();
  setup();  // clear all entries
  while (true) { // loop reading entry lines
    std::string tag = FR.GetTag();
    if (tag == "}") break;  // end of list
    if (tag != "HashIndexStoreFind") {
      clipper::Message::message(Message_fatal
	("hash_table::Restore unrecognised tag "+tag));
    }
    int i = FR.Int();
    Nstore_list.at(i) = FR.Int();
    Nfind_list.at(i)  = FR.Int();
  }
}
