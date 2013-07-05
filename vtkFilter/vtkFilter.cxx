/**
 * \file vtkFilter.cxx
 * \brief Contains main functionality of the vtkFilter
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#include <iostream>
#include <string>
#include "vtkReader.h"
#include "vtkWriter.h"
#include "fiber.h"
#include "Region.h"
#include "ExpressionParser.h"
#include "ExpressionEvaluator.h"
#include "vtkFilterCLP.h"

int main(int argc, char* argv[]) {


  // PARSE COMMAND LINE OR SLICER INPUT ////////////////////////////////////
  // Defined in GenerateCLP
  PARSE_ARGS ;

  // Check input TODO

  // ALLOCATE THE FIBERS ///////////////////////////////////////////////////
  std::vector<Fiber> in_fibers;
  std::vector<Fiber> out_fibers;

  // LOAD THE REGIONS //////////////////////////////////////////////////////
  Region regionA(LabelFileA, LabelA);
  Region regionB(LabelFileB, LabelB);
  Region regionC(LabelFileC, LabelC);
  Region regionD(LabelFileD, LabelD);

  // SET THINGS UP AND RUN ////////////////////////////////////////////////

  // Read the VTK Fiber
  if (Verbose) std::cout << "** Reading VTK file...\n";
  vtkReader<Fiber> * reader = new vtkReader<Fiber>();
  reader->SetInputPath(FiberInFile);
  reader->SetOutputFibers(in_fibers);
  reader->SetReadFieldData(CopyFields);
  reader->SetVerbose(Verbose);
  reader->Run();
  delete reader;
  if (Verbose) std::cout << "-Number of fibers in the input: " << in_fibers.size() << std::endl;

  // Parse end convert the input logic expression to postifx
  ExpressionParser * expParser = new ExpressionParser();
  expParser->SetInput(LogicFunction);
  expParser->SetVerbose(Verbose);
  expParser->Run();
  std::string sPostfixExpr = expParser->GetPostfix();

  // Filter out the fibers that don't match the conditions
  if (Verbose) std::cout << "** Filtering out Fibers...\n";

  ExpressionEvaluator * evaluator = new ExpressionEvaluator();
  if (!regionA.IsEmpty())
    evaluator->SetRegionA(regionA);
  if (!regionB.IsEmpty())
    evaluator->SetRegionB(regionB);
  if (!regionC.IsEmpty())
    evaluator->SetRegionC(regionC);
  if (!regionD.IsEmpty())
    evaluator->SetRegionD(regionD);
  evaluator->SetVerbose(Verbose);
  evaluator->SetPostfixExpr(sPostfixExpr);
  evaluator->SetInputFibers(in_fibers);
  evaluator->SetOutputFibers(out_fibers);
  evaluator->Run();

  if (Verbose) std::cout << "-Number of fibers in the output: " << out_fibers.size() << std::endl;

  if (out_fibers.size() == 0) {
    std::cout << "** There are no Fibers that match the expression: " << LogicFunction << ". Not Writing the output file.\n";
    return EXIT_SUCCESS;
  }

  // Write output files
  if (Verbose) std::cout << "** Writing VTK File...\n";
  vtkWriter * writer = new vtkWriter();
  writer->SetOutputPath(FiberOutFile);
  writer->SetInputFibers(out_fibers);
  writer->Run();
  delete writer;

  if (Verbose) std::cout << "** done\n";

  return EXIT_SUCCESS;

}
