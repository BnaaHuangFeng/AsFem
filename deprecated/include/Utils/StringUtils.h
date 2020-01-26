#ifndef ASFEM_STRINGUTILS_H
#define ASFEM_STRINGUTILS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>

using namespace std;

string RemoveStrSpace(string &instr);
string RemoveSymbolFromStr(string &instr,char symbol);
vector<string> SplitStr(string &instr,char symbol);

bool IsUniqueStrVec(vector<string> &strvec);

bool IsCommentLine(string &instr);

bool IsBracketMatch(ifstream &in,const int &linenum0);
bool IsBracketMatch(ifstream &in,const int &linenum0,int &lastend_linenum);

// for string to number convert
vector<double> SplitStrNum(string &instr);
vector<double> SplitStrNumAfter(string instr,int pos);

void GoToLine(ifstream &in,const int &linenum);

// for time dependent dirichlet bc
bool IsValidExpression(string str);

#endif //ASFEM_STRINGUTILS_H