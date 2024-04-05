#ifndef __SDMMannerH__
#define __SDMMannerH__
#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;

unsigned char** occupancyData;
int occupancyHeight;
int occupancyWidth;
bool                QPr5;
int                 OorGorA;
ofstream          extraFeatures;
map<string, string> xyd_features;
map<string, int>    frontModeFlag; 

void occupancyDciInit(int height, int width, int frameNum, string occupancyName);
void checkData(int frameNum, int height, int width);

void imgClassicate(string name);
#endif