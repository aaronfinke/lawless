// CCP4base.cpp
//
// Written by Airlie McCoy for Phaser
//(c) 2000-2005 Cambridge University Technical Services Ltd
//All rights reserved
#include "jiffy.hh"
#include "CCP4base.hh"

namespace phaser_io {

void CCP4base::Add_Key(std::string card)
{ possible_keys.push_back(card); }

Token_value CCP4base::skip_line(std::istringstream& input_stream) 
{
  //This routine returns either END or ENDLINE after skipping to the end 
  //of a line, warning of any unexpected tokens in the process.
  size_t keywords_size = keywords.size();
  if (curr_tok == ENDLINE || curr_tok == END)
    return curr_tok;
  get_token(input_stream);
  if (curr_tok == ENDLINE || curr_tok == END)
    return curr_tok;
  else
  {
    do { get_token(input_stream); }
    while (curr_tok != ENDLINE && curr_tok != END);
    size_t keywords_size2 = keywords.size();
    std::string remainder = keywords.substr(keywords_size,keywords_size2-keywords_size);
    if (annotate_keywords)
      keywords += "\tWarning - Additional tokens ignored: " + remainder;
  }
  return curr_tok;
}

Token_value CCP4base::get_token(std::istringstream& input_stream)
{
  std::string tmp;
  bool test_number = true;
  bool found_digit = false;

  char ch;
  do { if (!input_stream.get(ch)) 
       {
         keywords += '\n';
         return curr_tok = END;
       }
       else
         keywords += ch;
     } while ((std::isspace)(ch) && ch != '\n'); //skip space, tabs, feedforms


  //if line continuation, skip all spaces including \n
  if ((ch == '&') || (ch == '\\'))
    do { if (!input_stream.get(ch)) 
         {
           keywords += '\n';
           return curr_tok = END;
         }
         else
           keywords += ch;
       } while ((std::isspace)(ch)); 

  switch (ch)
  {
    case '#': case '!':
      std::getline(input_stream,string_value,'\n');
      keywords += string_value;
      keywords += '\n';
      return curr_tok=ENDLINE;
      break;

    case '\n':
      return curr_tok=ENDLINE;
      break;
    
    case '=': 
      return curr_tok=ASSIGN;
      break; 
   
    case '+': case '-': case '.':
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':    
      found_digit = std::isdigit(ch);
      string_value = ch;
      // Obsolete: allow digits and dots only for subsequent chars with 
      // while(input_stream.get(ch) && (isdigit(ch) || ch=='.')) here
      while(input_stream.get(ch) && (std::isprint)(ch) && !(std::isspace)(ch) && ch != '=' && ch !='#') 
      {
//rest of string is all digits or . or e or E or + or -(exponent)
        if (!((std::isdigit)(ch) || ch=='.' || ch=='e' || ch=='E' || ch=='+'|| ch=='-')) 
          test_number = false;
        string_value += ch; 
	if (std::isdigit(ch)) found_digit = true;
      }
      input_stream.putback(ch); // oops - read one too far
      tmp = string_value;
      keywords += tmp.erase(0,1);// take off already added ch from front of keywords string
      if (!found_digit) test_number = false;  // a number must contain at least one digit
      if (test_number)
      {
        number_value = atof(string_value.c_str());
        return curr_tok = NUMBER; //at end, string_value is the number too
      }
      else //this allows std::strings starting with a any character except = or #
        return curr_tok = NAME;
      break; 

    default:
      if ((std::isprint)(ch))
      {
        string_value = ch;
        //allow printable characters except spaces and ASSIGN
        while(input_stream.get(ch) && (std::isprint)(ch) && !(std::isspace)(ch) && ch != '=') 
          string_value += ch;
        input_stream.putback(ch); // oops - read one too far
        tmp = string_value;
        keywords += tmp.erase(0,1);
        return curr_tok = NAME;
      }
    }

    throw SyntaxError(keywords,"Unrecognised token (" + std::string(1,ch) + ")");
    return curr_tok = ENDLINE;
}

Token_value CCP4base::get_key(std::istringstream& input_stream)
{
  get_token(input_stream);
  
  switch (curr_tok)
  {
    case NAME:
    {
      Preprocessor p; //temporary, just to call end_keys
      std::vector<std::string> end_keys = p.getEndKeys();
      for (size_t i = 0; i < end_keys.size(); i++)
        if (!stoup(string_value).find(end_keys[i])) 
          return curr_tok = END;
      for (size_t i = 0; i < possible_keys.size(); i++) 
        if (keyIs(possible_keys[i]))
          return curr_tok = possible_fns[i]->parse(input_stream);
    }
    break;
     
    case ENDLINE:
      return ENDLINE;
      break;
   
    case END:
      return END;
      break;
    
    default:
      ; 
  }
  //Keyword not recognised
  std::string key = stoup(string_value);
  //This routine returns either END or ENDLINE after skipping to the end
  //of a line
  do { get_token(input_stream); }
  while (curr_tok != ENDLINE && curr_tok != END);
 // if (annotate_keywords && key.find('@')) //i.e. ignore preprocessor
  if (annotate_keywords)
  {
    size_t erase_end = keywords.rfind("\n");
    size_t erase_start =  keywords.rfind("\n",erase_end-1);
    if (erase_end != std::string::npos && erase_start != std::string::npos)
    {
      size_t erase_length = erase_end-erase_start;
      keywords.erase(erase_start,erase_length);
   // keywords += "\tWarning - Keyword " + key + " not relevant\n";
    }
  }
  return curr_tok;
}

bool CCP4base::getBoolean(std::istringstream& input_stream)
{
  compulsoryKey(input_stream,2,"ON","OFF");
  if (keyIs("ON")) return true;
  return false;
}

std::string CCP4base::getFileName(std::istringstream& input_stream)
{
 //allows =,#,! in filenames
  //skip leading spaces (there must be leading spaces)
  for (;;) 
  {
    char ch;
    if (input_stream.get(ch)) 
    {
      if (ch == '\n')
      {
        input_stream.putback(ch); //so that ENDLINE token is found
        return "";
      }
      else if ((std::isspace)(ch))
      {
        keywords += ch; //accept and ignore space
      }
      else //it is an interesting char and needs to be read again
      {
        input_stream.putback(ch);
        break;
      }
    }
    else
    {
      input_stream.putback(ch); //so that END can be found
      return "";
    }
  } 
  std::string filename;
  for (;;) 
  {
    char ch;
    if (input_stream.get(ch)) 
    {
      if ((std::isspace)(ch))
      {
        input_stream.putback(ch); //so that ENDLINE token is found
        return filename;
      }
      else 
      {
        keywords += ch; 
        filename += ch;
      }
    }
    else
    {
      input_stream.putback(ch); //so that END can be found
      return filename;
    }
  } 
  //  return filename; //should never be reached
}

std::string CCP4base::getLine(std::istringstream& input_stream)
{
  std::string line;
  std::getline(input_stream,line,'\n');
  keywords += line + '\n'; //have to add this explicitly as it is not going through the tokenizer
  return line;
}

void CCP4base::getAssign(std::istringstream& input_stream)
{
  if (get_token(input_stream) != ASSIGN)
  throw SyntaxError(keywords,"Use = for assignment");
}

std::string CCP4base::getString(std::istringstream& input_stream)
{
  if (!tokenIs(input_stream,2,NAME,NUMBER))
  throw SyntaxError(keywords,"Character string not present or not valid");
  return string_value;
}

floatType CCP4base::get1num(std::istringstream& input_stream)
{
  if (!tokenIs(input_stream,1,NUMBER))
  throw SyntaxError(keywords,"Number not present or not valid");
  return number_value;
}

DVect3 CCP4base::get3nums(std::istringstream& input_stream)
{
  DVect3 temp;
  temp[0] = get1num(input_stream);
  temp[1] = get1num(input_stream);
  temp[2] = get1num(input_stream);
  return temp;
}

DMat33 CCP4base::get9nums(std::istringstream& input_stream)
{
  DMat33 temp;
  temp(0,0) = get1num(input_stream);
  temp(0,1) = get1num(input_stream);
  temp(0,2) = get1num(input_stream);
  temp(1,0) = get1num(input_stream);
  temp(1,1) = get1num(input_stream);
  temp(1,2) = get1num(input_stream);
  temp(2,0) = get1num(input_stream);
  temp(2,1) = get1num(input_stream);
  temp(2,2) = get1num(input_stream);
  return temp;
}

bool CCP4base::keyIs(std::string key)
{
  //use first four characters of key for comparison
  int len(4);
  //except in special cases
  if (key == "NUM" || key == "SIG" || key == "ROT" || key == "TRA"
      || key == "MIN" || key == "MAX" || key == "SEQ" || key == "PHI") len = 3;
  if (key == "F") len = 1;
  std::string upper_string_value(stoup(string_value));
  if (int(string_value.size()) >= len) {
    for (size_t i = 0; i < unsigned(len) && i < key.size(); i++)
      if (upper_string_value[i] != key[i])
        return false;
    return true;
  }
  return (upper_string_value == key);
}

bool CCP4base::tokenIs(std::istringstream& input_stream,int va_len, ...)
{
  get_token(input_stream);
  bool tokenFound(false);
  va_list tokens;
  va_start(tokens,va_len);
  //va_arg return type Token_value, but the compiler wants this promoted to int
  for (int i = 0; i < va_len; i++)
    if (curr_tok == va_arg(tokens,int))
      tokenFound = true;
  va_end(tokens);
  return tokenFound;
}

bool CCP4base::tokenIs(int va_len, ...)
{
  bool tokenFound(false);
  va_list tokens;
  va_start(tokens,va_len);
  //va_arg return type Token_value, but the compiler wants this promoted to int
  for (int i = 0; i < va_len; i++)
    if (curr_tok == va_arg(tokens,int))
      tokenFound = true;
  va_end(tokens);
  return tokenFound;
}

bool CCP4base::compulsoryKey(int va_len, ...)
{
  bool keyFound(false);
  va_list keys;
  va_start(keys,va_len);
  if (tokenIs(1,NAME))
  {
//va_arg will not accept std::string
    for (int i = 0; i < va_len; i++)
      if (keyIs(std::string(va_arg(keys,char*))))
        keyFound = true;
  }
  else  
    keyFound = false;
  va_end(keys);

  if (!(keyFound))
  {
    va_start(keys,va_len);
    std::string useMessage = "Use ";
    if (va_len > 1)
    {
      for (int i = 0; i < va_len-1; i++)
        useMessage += std::string(va_arg(keys,char*)) + " ";
      useMessage += "or " + std::string(va_arg(keys,char*));
    }
    else
      useMessage += std::string(va_arg(keys,char*));
     
    va_end(keys);
    throw SyntaxError(keywords,useMessage);
  }
  return keyFound;
}

bool CCP4base::compulsoryKey(std::istringstream& input_stream,int va_len, ...)
{
  get_token(input_stream);
  bool keyFound(false);
  va_list keys;
  va_start(keys,va_len);
  if (tokenIs(1,NAME))
  {
//va_arg will not accept std::string
    for (int i = 0; i < va_len; i++)
      if (keyIs(std::string(va_arg(keys,char*))))
        keyFound = true;
  }
  else  
    keyFound = false;
  va_end(keys);

  if (!(keyFound))
  {
    va_start(keys,va_len);
    std::string useMessage = "Use ";
    if (va_len > 1)
    {
      for (int i = 0; i < va_len-1; i++)
        useMessage += std::string(va_arg(keys,char*)) + " ";
      useMessage += "or " + std::string(va_arg(keys,char*));
    }
    else
      useMessage += std::string(va_arg(keys,char*));
     
    va_end(keys);
    throw SyntaxError(keywords,useMessage);
  }
  return keyFound;
}

bool CCP4base::optionalKey(std::istringstream& input_stream,int va_len, ...)
{
  get_token(input_stream);
  bool keyFound(false);
  va_list keys;
  va_start(keys,va_len);
  if (tokenIs(1,NAME))
  {
//va_arg will not accept std::string
    for (int i = 0; i < va_len; i++)
      if (keyIs(std::string(va_arg(keys,char*))))
        keyFound = true;
  }
  va_end(keys);

  //we expect ENDLINE or END if none of the keywords are found
  if ((tokenIs(1,NAME) && !keyFound) || !tokenIs(3,NAME,ENDLINE,END))
  {
    va_start(keys,va_len);
    std::string useMessage = "Use (optionally) ";
    if (va_len > 1)
    {
      for (int i = 0; i < va_len - 1; i++)
        useMessage += std::string(va_arg(keys,char*)) + " ";
      useMessage += "or " + std::string(va_arg(keys,char*));
    }
    else
      useMessage += std::string(va_arg(keys,char*));
     
    va_end(keys);
    throw SyntaxError(keywords,useMessage);
  }
  return keyFound;
}

bool CCP4base::optionalKey(int va_len, ...)
{
  bool keyFound(false);
  va_list keys;
  va_start(keys,va_len);
  if (tokenIs(1,NAME))
  {
    for (int i = 0; i < va_len; i++)
    {
//va_arg will not accept std::string
      std::string lookup = std::string(va_arg(keys,char*));
      if (keyIs(lookup)) keyFound = true;
    }
  }
  else
    keyFound = false;
  va_end(keys);

  //we expect ENDLINE or END if none of the keywords are found
  if ((tokenIs(1,NAME) && !keyFound) || !tokenIs(3,NAME,ENDLINE,END))
  {
    va_start(keys,va_len);
    std::string useMessage = "Use (optionally) ";
    if (va_len > 1)
    {
      for (int i = 0; i < va_len - 1; i++)
        useMessage += std::string(va_arg(keys,char*)) + " ";
      useMessage += "or " + std::string(va_arg(keys,char*));
    }
    else
      useMessage += std::string(va_arg(keys,char*));
     
    va_end(keys);
    throw SyntaxError(keywords,useMessage);
  }

  return keyFound;
}

void CCP4base::parseCCP4(Preprocessor capture)
{
  std::string capturestring = capture.Echo();
  //std::istringstream input_stream(const_cast<char*>(capturestring.c_str()),capturestring.size());
  std::istringstream input_stream(const_cast<char*>(capturestring.c_str()));

  // Private initial values
  curr_tok = ENDLINE;

  // Do the parsing
  input_stream.seekg(0);
  PROGRAM_ASSERT(possible_keys.size() == possible_fns.size());
  while (input_stream)
  {
    get_key(input_stream);
    if (curr_tok == END) break;
    if (curr_tok == ENDLINE) continue;
  }
}


floatType CCP4base::isperc(floatType percent)
{
  if (percent < 0) throw InputError(keywords,"Value not a percentage");
  if (percent > 1.0 && percent <= 100.0) percent /= 100;
  if (percent > 1) throw InputError(keywords,"Value not a percentage");
  return percent;
}


std::string CCP4base::Keywords()
{ return keywords; }

}//phaser_io

