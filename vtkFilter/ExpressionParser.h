/**
 * \file ExpressionParser.h
 * \brief Contains class ExpressionParser for converting logical Infix expressions to Postfix expressions.
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#ifndef EXPRESSIONPARSER_H_
#define EXPRESSIONPARSER_H_

#include <string>
#include <vector>

/**
 * \class ExpressionParser
 * \brief Converts a logical expression in infix to postfix
*/
class ExpressionParser
{
public:

  /** Constructor */
  ExpressionParser();

  /** Destructor */
  ~ExpressionParser();

  /** Run the expression parser when all parameters are set */
  bool Run();

  /** Set input infix notation */
  void SetInput(std::string s)
  {
    _sInfix = s;
  }

  /** Set whether to output execution details in console */
  void SetVerbose(bool b)
  {
    _bVerbose = b;
  }

  /** Return the postfix Notation after running the expression parser */
  std::string GetPostfix()
  {
    return _sPostfix;
  }

protected:

  /** The infix notation */
  std::string _sInfix;

  /** The postfix notation */
  std::string _sPostfix;

  /** 'Look-up table' for priority of operators */
  int priority(char op);

  /** Clean up the string a little bit before parsing it */
  void formatString(std::string & str);

  /** Check if the current char is an operand */
  bool isOperand(char c);

  /** Check if the curent char in an operator */
  bool isOperator(char c);

  /** Check if it is a valid expression */
  bool expresssionValid(const std::string & str);

  /** Parse the string for 'not' or '!' */
  void parseStringForNot(std::string & str);

  /** Toggles command line output */
  bool _bVerbose;

};

#endif  // EXPRESSIONPARSER_H_
