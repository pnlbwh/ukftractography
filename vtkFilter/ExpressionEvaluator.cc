/**
 * \file ExpressionEvaluator.cc
 * \brief Implements class ExpressionEvaluator
*/

#include <iostream>
#include <ctype.h>
#include "ExpressionEvaluator.h"

ExpressionEvaluator::ExpressionEvaluator()
{
  _filter = new FiberFilter();
  _logic = new FiberLogic();
  _bVerbose = false;
}

ExpressionEvaluator::~ExpressionEvaluator()
{
  delete _filter;
}

bool ExpressionEvaluator::Run()
{

  for (unsigned int i = 0; i < _sPostfixExpr.size(); ++i) {
    if (isOperand(_sPostfixExpr[i])) {
      if (_bVerbose) std::cout << "-filtering " << _sPostfixExpr[i] << "...\n";
      _operandStack.push_back(evaluateOperand(_sPostfixExpr[i]));
    } else if (isOperator(_sPostfixExpr[i])) {
      std::vector<Fiber> A = _operandStack.back();
      _operandStack.pop_back();

      std::vector<Fiber> B = _operandStack.back();
      _operandStack.pop_back();

      if (_bVerbose) std::cout << "-doing '" << _sPostfixExpr[i] << "'-logic...\n";
      if (_sPostfixExpr[i] == '|') {
        _operandStack.push_back(_logic->doOr(A, B));
      } else if (_sPostfixExpr[i] == '&') {
        _operandStack.push_back(_logic->doAnd(A, B));
      }
    }
  }

  *_outFibers = _operandStack[0];

  return 0;

}

// SETTERS
void ExpressionEvaluator::SetRegionA(const Region & reg)
{
  _regionA = & reg;
}

void ExpressionEvaluator::SetRegionB(const Region & reg)
{
  _regionB = & reg;
}

void ExpressionEvaluator::SetRegionC(const Region & reg)
{
  _regionC = & reg;
}

void ExpressionEvaluator::SetRegionD(const Region & reg)
{
  _regionD = & reg;
}

void ExpressionEvaluator::SetInputFibers(const std::vector< Fiber > & fibers)
{
  _inFibers = & fibers;
}

void ExpressionEvaluator::SetOutputFibers(std::vector< Fiber > & fibers)
{
  _outFibers = & fibers;
}

void ExpressionEvaluator::SetPostfixExpr(const std::string & str)
{
  _sPostfixExpr = str;
}

void ExpressionEvaluator::SetVerbose(const bool & b)
{
  _bVerbose = b;
}

std::vector<Fiber> ExpressionEvaluator::evaluateOperand(char operand)
{

  std::vector<Fiber> out_fibers;

  _filter->SetInputFibers(*(_inFibers)); // <-- all fibers
  _filter->SetOutputFibers(out_fibers); // <-- new fibers

  if (isupper(operand))
    _filter->SetCalcInverse(false);
  else {
    _filter->SetCalcInverse(true);
    operand -= 32;
  }

  if (isInAtoD(operand)) { // means passing not ending;
    _filter->SetConnectionMode(FiberFilter::PASS);
  } else if (isInEtoH(operand)) {
    _filter->SetConnectionMode(FiberFilter::END);
    operand -= 4;
  }

  switch(operand) {
  case 'A':
    _filter->SetRegion(*_regionA);
    break;
  case 'B':
    _filter->SetRegion(*_regionB);
    break;
  case 'C':
    _filter->SetRegion(*_regionC);
    break;
  case 'D':
    _filter->SetRegion(*_regionD);
    break;
  default :
    std::cout << "error!\n";
    exit(1);
  }

  _filter->Run();

  return out_fibers;

}


bool ExpressionEvaluator::isInEtoH(char c)
{
  if ((c >= 69 && c <= 72) || (c >= 69 + 32 && c <= 72 + 32))
    return true;
  return false;
}

bool ExpressionEvaluator::isInAtoD(char c)
{
  if ((c >= 69 - 4 && c <= 72 - 4) || (c >= 69 + 32 - 4 && c <= 72 + 32 - 4))
    return true;
  return false;
}

bool ExpressionEvaluator::isOperand(char c)
{
  bool out = c == 'A' || c == 'B' || c == 'C' || c == 'D' || c == 'a' || c == 'b' || c == 'c' || c == 'd' ||
             c == 'E' || c == 'F' || c == 'G' || c == 'H' || c == 'e' || c == 'f' || c == 'g' || c == 'h';
  return out;
}

bool ExpressionEvaluator::isOperator(char c)
{
  return  c == '&' || c == '|' || c == '!';
}

