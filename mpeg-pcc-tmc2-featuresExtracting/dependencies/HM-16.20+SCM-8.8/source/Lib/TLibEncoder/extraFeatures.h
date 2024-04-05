#ifndef __extraFeaturesH__
#define __extraFeaturesH__
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
using namespace std;

ofstream            extraGeoMergeFeatures;
ofstream            extraAttriMergeFeatures;
ofstream            extraGeoInterFeatures;
ofstream            extraAttriInterFeatures;
map<string, string> xyd_GeoMergeFeatures;
map<string, string> xyd_AttriMergeFeatures;
map<string, string> xyd_GeoInterFeatures;
map<string, string> xyd_AttriInterFeatures;
map<string, int> frontModeFlag;

void initExtraFeatures(int QP);
void destroyExtraFeatures();
#endif