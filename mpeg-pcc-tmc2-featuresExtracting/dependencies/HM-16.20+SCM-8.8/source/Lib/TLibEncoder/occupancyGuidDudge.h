#ifndef __occupancyDciDudgeH__
#define __occupancyDciDudgeH__
#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;

unsigned char** occupancyData;
int             occupancyHeight;
int             occupancyWidth;
int             occupancyFrameNumber;
bool            QPr5;
int             OorGorA;

void occupancyDciInit( int height, int width, int frameNum, string occupancyName );
void checkData( int frameNum, int height, int width );
#endif