#ifndef __extraFeaturesCPP__
#define __extraFeaturesCPP__ 
#include "extraFeatures.h"

extern ofstream                   extraGeoMergeFeatures;
extern ofstream                   extraAttriMergeFeatures;
extern ofstream                   extraGeoInterFeatures;
extern ofstream                   extraAttriInterFeatures;
extern map<string, string>        xyd_GeoMergeFeatures;
extern map<string, string>        xyd_AttriMergeFeatures;
extern map<string, string>        xyd_GeoInterFeatures;
extern map<string, string>        xyd_AttriInterFeatures;
extern int                        OorGorA;

void initExtraFeatures(int QP) {
  xyd_GeoMergeFeatures.clear();
  xyd_AttriMergeFeatures.clear();
  xyd_GeoInterFeatures.clear();
  xyd_AttriInterFeatures.clear();
  frontModeFlag.clear();
  stringstream sm_gpath, sm_apath, si_gpath, si_apath, sp_gpath, sp_apath;
  sm_gpath << "../__extraFeatures/oriP_extraFeatures_Geo.csv";
  sm_apath << "../__extraFeatures/oriP_extraFeatures_Att.csv";
  si_gpath << "../__extraFeatures/oriI_extraFeatures_Geo.csv";
  si_apath << "../__extraFeatures/oriI_extraFeatures_Att.csv";
  string gmpath, ampath, gipath, aipath;
  sm_gpath >> gmpath;
  sm_apath >> ampath;
  si_gpath >> gipath;
  si_apath >> aipath;
  extraGeoMergeFeatures.open( gmpath, ios::app );
  extraAttriMergeFeatures.open( ampath, ios::app );
  extraGeoInterFeatures.open( gipath, ios::app );
  extraAttriInterFeatures.open( aipath, ios::app );
}

void destroyExtraFeatures() {
  xyd_GeoMergeFeatures.clear();
  xyd_AttriMergeFeatures.clear();
  xyd_GeoInterFeatures.clear();
  xyd_AttriInterFeatures.clear();
  frontModeFlag.clear();
  extraGeoMergeFeatures.close();
  extraAttriMergeFeatures.close();
  extraGeoInterFeatures.close();
  extraAttriInterFeatures.close();
}
#endif