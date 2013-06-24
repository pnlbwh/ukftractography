/**
 * \file ExpressionParser.h
 * \brief Contains implementation of class Expression Parser
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#include "ExpressionParser.h"
#include <iostream>

ExpressionParser::ExpressionParser()
{
}

bool ExpressionParser::Run()
{

  std::string tInfix = _sInfix;

  formatString(tInfix);
  if( !expresssionValid(tInfix) )
    {
    return 1;
    }
  parseStringForNot(tInfix);

  std::vector<char> stack;
  _sPostfix = "";

  char * i;
  i = &tInfix[0];
  for( unsigned int j = 0; j < tInfix.size(); j++ )
    {
    while( *i == ' ' ||  *i == '\t' )
      {
      i++;
      }

    if( isOperand(*i) )
      {
      _sPostfix.push_back(*(i++) );
      }

    if( *i == '(' )
      {
      stack.push_back(*(i++) );
      }

    if( *i == ')' )
      {
      while( stack.back() != '(' )
        {
        _sPostfix.push_back(stack.back() );
        stack.pop_back();
        }

      i++;
      }

    if( isOperator(*i) )
      {
      if( stack.empty() )
        {
        stack.push_back(*i);
        }
      else
        {
        while( priority(stack.back() ) >= priority(*i) )
          {
          _sPostfix.push_back(stack.back() );
          stack.pop_back();
          }

        stack.push_back(*i);
        }
      i++;
      }
    }

  while( !stack.empty() )
    {
    if( stack.back() != '(' )
      {
      _sPostfix.push_back(stack.back() );
      }
    stack.pop_back();
    }

  if( _bVerbose )
    {
    std::cout << "-Postfix: ";
    for( unsigned int i = 0; i < _sPostfix.size(); i++ )
      {
      std::cout << _sPostfix[i];
      }
    std::cout << std::endl;
    }

  return 0;
}

int ExpressionParser::priority(char op)
{
  switch( op )
    {
    case '&':
      {
      return 2;
      }
    case '|':
      {
      return 1;
      }
    default:
      return -1;
    }
}

void ExpressionParser::formatString(std::string & str)
{
  size_t position = str.find( "and" );

  while( position != std::string::npos )
    {
    str.replace( position, 3, "&" );
    position = str.find( "and");
    }

  position = str.find( "&&" );
  while( position != std::string::npos )
    {
    str.replace( position, 2, "&" );
    position = str.find( "&&", position + 1 );
    }

  position = str.find( "or" );
  while( position != std::string::npos )
    {
    str.replace( position, 2, "|" );
    position = str.find( "or", position + 1 );
    }

  position = str.find( "||" );
  while( position != std::string::npos )
    {
    str.replace( position, 2, "|" );
    position = str.find( "||", position + 1 );
    }

  position = str.find( "not " );
  while( position != std::string::npos )
    {
    str.replace( position, 4, "!" );
    position = str.find( "not ", position + 1 );
    }

  position = str.find( "not" );
  while( position != std::string::npos )
    {
    str.replace( position, 3, "!" );
    position = str.find( "not", position + 1 );
    }
  for( unsigned int i = 0; i < str.size(); i++ )
    {
    str[i] = toupper(str[i]);
    }
}

bool ExpressionParser::isOperand(char c)
{
  bool out = c == 'A' || c == 'B' || c == 'C' || c == 'D' || c == 'a' || c == 'b' || c == 'c' || c == 'd' ||
    c == 'E' || c == 'F' || c == 'G' || c == 'H' || c == 'e' || c == 'f' || c == 'g' || c == 'h';

  return out;
}

bool ExpressionParser::isOperator(char c)
{
  return c == '&' || c == '|' || c == '!';
}

bool ExpressionParser::expresssionValid(const std::string & str)
{

  int open = 0;
  int close = 0;

  for( unsigned int i = 0; i < str.size(); i++ )
    {
    if( str[i] == '(' )
      {
      open++;
      }
    else if( str[i] == ')' )
      {
      close++;
      }
    else if( !isOperand(str[i]) && !isOperator(str[i]) && str[i] != ' ' )
      {
      return 0;
      }
    }

  if( open != close )
    {
    return 0;
    }

  return 1;
}

void ExpressionParser::parseStringForNot(std::string & str)
{
  size_t position = str.find( "E" ); // find first space

  while( position != std::string::npos )
    {
    std::string repstring;
    repstring.push_back(str[position - 1] + 4);
    str.replace( position - 1, 2, repstring );
    position = str.find( "E", position + 1 );
    }

  position = str.find( "!" ); // find first space
  while( position != std::string::npos )
    {
    std::string repstring;
    repstring.push_back(str[position + 1] + 32); // to lowercase
    str.replace( position, 2, repstring );
    position = str.find( "!", position + 1 );
    }
}
