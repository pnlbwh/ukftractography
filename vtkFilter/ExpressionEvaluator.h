/**
 * \file ExpressionEvaluator.h
 * \brief Contains class ExpressionEvaluator for evaluating logical expressions on fibers
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#ifndef EXPRESSIONEVALUATOR_H_
#define EXPRESSIONEVALUATOR_H_

#include <vector>
#include <string>
#include "fiber.h"
#include "FiberFilter.h"
#include "FiberLogic.h"

/**
 * \class ExpressionEvaluator
 * Evaluates a logical expression involving fibers
*/
class ExpressionEvaluator
{
public:

  /** Constructor */
  ExpressionEvaluator();

  /** Destructor */
  ~ExpressionEvaluator();

  /** Run the expression evaluator once all the parameters are set */
  bool Run();

  /** Set region A */
  void SetRegionA(const Region & reg);

  /** Set region B */
  void SetRegionB(const Region & reg);

  /** Set region C */
  void SetRegionC(const Region & reg);

  /** Set region D */
  void SetRegionD(const Region & reg);

  /** Chose whether to display execution infos in the console */
  void SetVerbose(const bool & b);

  /** Set pointer to input fibers */
  void SetInputFibers(const std::vector<Fiber> & fibers);

  /** Set pointer to output fibers */
  void SetOutputFibers(std::vector<Fiber> & fibers);

  /** Set the posfix expression */
  void SetPostfixExpr(const std::string & str);

  /** Pointer to the output Fibers */
  std::vector<Fiber> * _outFibers;
protected:

  /** Pointer to region A */
  const Region * _regionA;

  /** Pointer to region B */
  const Region * _regionB;

  /** Pointer to region C */
  const Region * _regionC;

  /** Pointer to region D */
  const Region * _regionD;

  /** Pointer to input fibers */
  const std::vector<Fiber> * _inFibers;

  /** Was Region A set? */
  bool _bASet;

  /** Was Region B set? */
  bool _bBSet;

  /** Was Region C set? */
  bool _bCSet;

  /** Was Region D set? */
  bool _bDSet;

  /** Verbose or silent execution? */
  bool _bVerbose;

  /** String holding the Postfix expression */
  std::string _sPostfixExpr;

  /** Stack of operands */
  std::vector<std::vector<Fiber> > _operandStack;

  /** A pointer to the fiber filter */
  FiberFilter * _filter;

  /** A pointer to class containing the logic functions for fibers */
  FiberLogic * _logic;

  // Note: For actually evaluating the logical expression we use a slightly different alphabet
  // Lower Case letters mean 'not'
  // Letters shifted by 4 ascii character mean 'ending in'

  /** Returns true if c is E,F,G, or H lower or upper case */
  bool isInEtoH(char c);

  /** Returns true if c is A,B,C, or D lower or upper case */
  bool isInAtoD(char c);

  /** Returns true if c is a operand */
  bool isOperand(char c);

  /** Returns true if c is a operator */
  bool isOperator(char c);

  /** Pops to elements from the opperand stack and uses the latest opertor to evaluate the elements */
  std::vector<Fiber> evaluateOperand(char operand);

};

#endif // EXPRESSIONEVALUATOR_H_
