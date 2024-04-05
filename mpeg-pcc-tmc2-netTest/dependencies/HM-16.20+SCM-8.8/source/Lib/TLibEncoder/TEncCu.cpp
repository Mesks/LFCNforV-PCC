/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2017, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TEncCu.cpp
    \brief    Coding Unit (CU) encoder class
*/

#include <stdio.h>
#include "TEncTop.h"
#include "TEncCu.h"
#include "TEncAnalyze.h"
#include "TLibCommon/Debug.h"

#include <cmath>
#include <algorithm>

#ifdef SDMTEST            // MesksCode
#define P_GEO_INPUT 3     
#define P_ATT_INPUT 3  
#define P_GEO_OUTPUT 1    
#define P_ATT_OUTPUT 1  
#define I_GEO_INPUT 3  
#define I_ATT_INPUT 3   
#define I_GEO_OUTPUT 1    
#define I_ATT_OUTPUT 1  


#define P_GEO_CONV_INPUT_SIZE 16     // PGeometry Frames conv input size
#define P_ATT_CONV_INPUT_SIZE 16     // PAttribute Frames conv input size
#define I_GEO_CONV_INPUT_SIZE 16     // IGeometry Frames conv input size
#define I_ATT_CONV_INPUT_SIZE 16     // IAttribute Frames conv input size
#define P_GEO_CONV_KERNEL_SIZE1 3     // PGeometry Frames conv kernel size
#define P_ATT_CONV_KERNEL_SIZE1 3     // PAttribute Frames conv kernel size
#define I_GEO_CONV_KERNEL_SIZE1 3     // IGeometry Frames conv kernel size
#define I_ATT_CONV_KERNEL_SIZE1 3     // IAttribute Frames conv kernel size
#define P_GEO_CONV_POOL_SIZE1  4      // PGeometry Frames pool size
#define P_ATT_CONV_POOL_SIZE1  4      // PAttribute Frames pool size
#define I_GEO_CONV_POOL_SIZE1  4      // IGeometry Frames pool size
#define I_ATT_CONV_POOL_SIZE1  4      // IAttribute Frames pool size
#define P_GEO_CONV_KERNEL_SIZE2 3      // PGeometry Frames dense of conv module size
#define P_ATT_CONV_KERNEL_SIZE2 3      // PAttribute Frames dense of conv module size
#define I_GEO_CONV_KERNEL_SIZE2 3      // IGeometry Frames dense of conv module size
#define I_ATT_CONV_KERNEL_SIZE2 3      // IAttribute Frames dense of conv module size
#define P_GEO_CONV_POOL_SIZE2 4        // PGeometry Frames pool size
#define P_ATT_CONV_POOL_SIZE2 4        // PAttribute Frames pool size
#define I_GEO_CONV_POOL_SIZE2 4        // IGeometry Frames pool size
#define I_ATT_CONV_POOL_SIZE2 4        // IAttribute Frames pool size

#define P_GEO_SLP_INPUT  5           // PGeometry Frames SLP input number
#define P_ATT_SLP_INPUT  5           // PAttribute Frames SLP input number
#define I_GEO_SLP_INPUT  5           // IGeometry Frames SLP input number
#define I_ATT_SLP_INPUT  5           // IAttribute Frames SLP input number
#define P_GEO_SLP_OUTPUT 5           // PGeometry Frames SLP output number
#define P_ATT_SLP_OUTPUT 5           // PAttribute Frames SLP output number
#define I_GEO_SLP_OUTPUT 5           // IGeometry Frames SLP output number
#define I_ATT_SLP_OUTPUT 5           // PAttribute Frames SLP output number

#define P_GEO_CONCATENATE_INPUT ( 1 + P_GEO_SLP_OUTPUT )  // PGeometry Frames CONCATENATE input number
#define P_ATT_CONCATENATE_INPUT ( 1 + P_ATT_SLP_OUTPUT )  // PAttribute Frames CONCATENATE input number
#define I_GEO_CONCATENATE_INPUT ( 1 + I_GEO_SLP_OUTPUT )  // IGeometry Frames CONCATENATE input number
#define I_ATT_CONCATENATE_INPUT ( 1 + I_ATT_SLP_OUTPUT )  // IAttribute Frames CONCATENATE input number
#define P_GEO_CONCATENATE_DENSE_SIZE1 10          // PGeometry Frames dense of CONCATENATE module size
#define P_ATT_CONCATENATE_DENSE_SIZE1 10          // PAttribute Frames dense of CONCATENATE module size
#define I_GEO_CONCATENATE_DENSE_SIZE1 10          // IGeometry Frames dense of CONCATENATE module size
#define I_ATT_CONCATENATE_DENSE_SIZE1 10          // IAttribute Frames dense of CONCATENATE module size
#define P_GEO_CONCATENATE_DENSE_SIZE2 5          // PGeometry Frames dense of CONCATENATE module size
#define P_ATT_CONCATENATE_DENSE_SIZE2 5          // PAttribute Frames dense of CONCATENATE module size
#define I_GEO_CONCATENATE_DENSE_SIZE2 5          // IGeometry Frames dense of CONCATENATE module size
#define I_ATT_CONCATENATE_DENSE_SIZE2 5          // IAttribute Frames dense of CONCATENATE module size
#define P_GEO_CONCATENATE_OUTPUT 1               // PGeometry Frames CONCATENATE output number
#define P_ATT_CONCATENATE_OUTPUT 1               // PAttribute Frames CONCATENATE output number
#define I_GEO_CONCATENATE_OUTPUT 1               // IGeometry Frames CONCATENATE output number
#define I_ATT_CONCATENATE_OUTPUT 1               // IAttribute Frames CONCATENATE output number

extern unsigned char** occupancyData;
extern int             occupancyHeight;
extern int             occupancyWidth;
extern bool            QPr5;
extern int             OorGorA;

int CUClassify( int Width, int Height, int Y, int X, int nowPOC ) {
  int numberOfOne = 0; 
  int shrink      = 0; 
  if ( QPr5 )
    shrink = 2;
  else
    shrink = 4;

  int oCUWidth  = Width / shrink;
  int oCUHeight = Height / shrink;
  int oCUY      = Y / shrink;
  int oCUX      = X / shrink;

  for ( int i = 0; i < oCUHeight; i++ )
    for ( int j = 0; j < oCUWidth; j++ )
      numberOfOne += ( occupancyData[nowPOC / 2][( oCUY + i ) * occupancyWidth + oCUX + j] - 0 );

  // checkCUoinfo(oCUWidth, oCUHeight, oCUY, oCUX,nowPOC);

  if ( numberOfOne == 0 )
    return 0;  // unoccupancy block
  else if ( numberOfOne == oCUWidth * oCUHeight )
    return 1;  // fill block
  else
    return 2;  // boundary block
}

extern map<string, int> frontModeFlag;
typedef double          input_dtype;
const unsigned int      P_G_hidden1_num = 10;
const unsigned int      P_A_hidden1_num = 10;
const unsigned int      P_G_hidden2_num = 5;
const unsigned int      P_A_hidden2_num = 5;
const unsigned int      P_G_input_node  = P_GEO_INPUT;
const unsigned int      P_A_input_node  = P_ATT_INPUT;
const unsigned int      P_G_output_node = P_GEO_OUTPUT;
const unsigned int      P_A_output_node = P_ATT_OUTPUT;
const unsigned int      I_G_hidden1_num = 10;
const unsigned int      I_A_hidden1_num = 10;
const unsigned int      I_G_hidden2_num = 5;
const unsigned int      I_A_hidden2_num = 5;
const unsigned int      I_G_input_node  = I_GEO_INPUT;
const unsigned int      I_A_input_node  = I_ATT_INPUT;
const unsigned int      I_G_output_node = I_GEO_OUTPUT;
const unsigned int      I_A_output_node = I_ATT_OUTPUT;

// PDCNNwithMLP
// P Geo
double P_G_convConvWeight1[P_GEO_CONV_KERNEL_SIZE1 * P_GEO_CONV_KERNEL_SIZE1]{
    2.4758632, 1.1896214, -0.91053516, -2.7898495, -4.7970176, -9.42424, 8.630292, 7.809705, 6.170876 };
double P_G_convConvBias1[1]{ -1.3561134 };
double P_G_convConvWeight2[P_GEO_CONV_KERNEL_SIZE2 * P_GEO_CONV_KERNEL_SIZE2]{
    -1.2488679, 1.7871522, 1.0782553, 1.1152016, -0.009114418, -0.22153492, -1.9935546, 0.029436084, -1.290348 };
double P_G_convConvBias2[1]{ -0.28994897 };
double P_G_concatenateWeight1[P_GEO_CONCATENATE_INPUT * P_GEO_CONCATENATE_DENSE_SIZE1]{
    -1.0759901,  -0.0492474,    -1.2397894,  1.0118258,   -1.4125019,  0.90340674,  -0.97320753, 0.80374414,
    -0.7616137,  -0.0018930855, 0.7799791,   -4.5741467,  -1.477468,   -0.89610386, -1.7796128,  -0.98001516,
    1.1095102,   0.303772,      1.0178773,   -4.06721,    0.027867977, 0.003425187, -1.2981775,  -0.15392792,
    0.707992,    -0.16404249,   0.5127644,   1.0735879,   0.3138606,   -0.77745146, 1.4843153,   0.47959185,
    2.1357474,   -1.118577,     -2.8583927,  -0.91797566, 0.5864384,   -1.3364332,  0.6348829,   0.5910162,
    1.7955865,   -2.1240559,    2.4041002,   -1.8683606,  2.833893,    -1.9442354,  1.9541028,   -2.9001398,
    1.9458418,   -2.2012632,    -1.2247157,  0.33230487,  -1.0574832,  0.74369615,  6.872281,    0.54041696,
    -0.47260016, -0.79331535,   -0.29382518, 0.48194712 };
double P_G_concatenateBias1[P_GEO_CONCATENATE_DENSE_SIZE1 * 1]{ -0.88153243, 0.4483849, -0.3846405,  0.83345264,
                                                                -0.512444,   0.8605206, -0.88916564, 0.29526424,
                                                                -1.0142366,  0.37927988 };
double P_G_concatenateWeight2[P_GEO_CONCATENATE_DENSE_SIZE2 * P_GEO_CONCATENATE_DENSE_SIZE1]{
    -1.1227528, -1.1436865, -1.1002206, 1.0814103, -1.1722418, 2.0769336, 2.1096964, 2.1676087, -1.9802585, 2.1004987,
    -1.2727028, -1.2326187, -1.2642188, 1.2624966, -1.2391735, 1.5183864, 1.4587961, 1.4360263, -1.5041107, 1.399252,
    -1.4833606, -1.3539084, -1.7025777, 1.5065544, -1.7255229, 1.5443676, 1.6061351, 1.5328381, -1.5579032, 1.5044162,
    -1.3242137, -1.2287745, -1.3379658, 1.2263037, -1.3333625, 2.094658,  1.9817677, 2.0148566, -1.9849694, 2.0861876,
    -1.1415348, -1.1879857, -1.1672837, 1.1826041, -1.2406297, 2.08304,   2.0813377, 2.09045,   -2.1304822, 2.0076292 };
double P_G_concatenateBias2[P_GEO_CONCATENATE_DENSE_SIZE2 * 1]{ -0.120900124, -0.07236237, -0.11709604, 0.088577695,
                                                                -0.04462634 };
double P_G_concatenateOutputWeight[P_GEO_CONCATENATE_DENSE_SIZE2 * P_GEO_CONCATENATE_OUTPUT]{
    -2.6094694, -2.6136334, -2.677842, 4.449585, -2.6705377 };
double P_G_concatenateOutputBias[P_GEO_CONCATENATE_OUTPUT * 1]{ 1.1744573 };

// P Att
double P_A_convConvWeight1[P_ATT_CONV_KERNEL_SIZE1 * P_ATT_CONV_KERNEL_SIZE1]{
    -16.964169, -11.51067, -7.325712, -4.9276814, -3.9299402, 0.3791383, -4.196525, -3.0311632, -2.1498566 };
double P_A_convConvBias1[1]{ -0.65761745 };
double P_A_convConvWeight2[P_ATT_CONV_KERNEL_SIZE2 * P_ATT_CONV_KERNEL_SIZE2]{
    0.99691933, 0.77970564, 0.85397446, 0.84983695, -0.71553457, -0.018739762, 0.5121835, -1.0602562, 0.35058513 };
double P_A_convConvBias2[1]{ -0.19663292 };
double P_A_concatenateWeight1[P_ATT_CONCATENATE_INPUT * P_ATT_CONCATENATE_DENSE_SIZE1]{
    -0.032710288, 1.4889127,   0.89384097, 1.4529092,   0.36860648,  0.91277385,  1.4129859,  0.16975632,  -0.013084309,
    0.9351847,    0.6767992,   -1.9461156, 2.4529018,   -2.646398,   0.97556126,  -1.3086708, -2.1214013,  8.445762,
    0.6555587,    -1.2670227,  -2.2544613, 0.5143039,   -1.7396946,  1.8503525,   -2.7503145, -1.5198746,  0.9946552,
    0.70373005,   -2.0999002,  -2.0333483, 0.44351536,  -2.3205783,  0.42681605,  1.4711133,  0.30958605,  -0.51860356,
    -2.0150208,   0.3188992,   0.3757428,  -0.56542516, -6.423156,   -0.2559604,  -7.3162575, 0.8165374,   -7.5795054,
    -0.6928014,   0.009441176, -3.6579366, -6.167086,   -0.5831627,  -0.17875731, 0.9002587,  -0.32242608, -0.17989476,
    -0.034986496, -0.4107071,  0.5728541,  -0.4161747,  -0.35338545, -0.6789865 };
double P_A_concatenateBias1[P_ATT_CONCATENATE_DENSE_SIZE1 * 1]{ -0.21181835, 1.4126166,  0.32600805, 1.1606296,
                                                                0.05444333,  0.49696076, 1.2324318,  1.2322121,
                                                                -0.18756117, 0.42534268 };
double P_A_concatenateWeight2[P_ATT_CONCATENATE_DENSE_SIZE2 * P_ATT_CONCATENATE_DENSE_SIZE1]{
    3.8335218,   -0.7105041,  1.5309472,  -1.2056398,  -1.4289868,  0.12845626,  -0.7915686,  2.0034022,  -1.7352916,
    -2.0423007,  7.099755,    -2.3044403, 0.47168055,  -0.63115233, -0.37021837, -0.37957335, -0.9435457, 2.1733468,
    -2.02117,    -2.2040381,  6.6456256,  -2.5257983,  1.2062473,   -1.3545538,  -1.3360053,  -0.3834816, -0.12382437,
    1.1053572,   -1.0038338,  -1.1702371, -0.18173823, -0.6022258,  2.0151198,   -1.7862194,  -1.9549915, -1.9637787,
    3.210268,    -3.2620373,  3.2225335,  3.2472556,   3.36856,     -0.3028352,  1.4985173,   -1.2248062, -1.4759295,
    -0.42794046, -0.15899517, 1.1180614,  -1.018939,   -1.1762383 };
double P_A_concatenateBias2[P_ATT_CONCATENATE_DENSE_SIZE2 * 1]{ -1.4279286, 1.2463415, -1.6314437, 1.4380723,
                                                                1.6176305 };
double P_A_concatenateOutputWeight[P_ATT_CONCATENATE_DENSE_SIZE2 * P_ATT_CONCATENATE_OUTPUT]{
    -3.8938975, 2.752624, -3.1005256, 4.6912317, 5.4514894 };
double P_A_concatenateOutputBias[P_ATT_CONCATENATE_OUTPUT * 1]{ -0.42919534 };

// I Geo
double I_G_convConvWeight1[I_GEO_CONV_KERNEL_SIZE1 * I_GEO_CONV_KERNEL_SIZE1]{
    -0.22954915, 1.9033159, -0.87728554, 2.4076807, 0.7516228, -2.1240177, -2.012418, -2.6249645, -1.8492402 };
double I_G_convConvBias1[1]{ -0.3274676 };
double I_G_convConvWeight2[I_GEO_CONV_KERNEL_SIZE2 * I_GEO_CONV_KERNEL_SIZE2]{
    0.640245, 0.19110294, 0.49701044, 0.048497852, -0.18459806, -0.0085841715, 1.0944448, -0.0833085, 0.30345392 };
double I_G_convConvBias2[1]{ -0.11067175 };
double I_G_concatenateWeight1[I_GEO_CONCATENATE_INPUT * I_GEO_CONCATENATE_DENSE_SIZE1]{
    -1.2156739, 0.9088411,   0.467494,   -0.8541831,  -0.8574022, -0.9639003, 0.87062454, -0.68620884, 1.4931751,
    -1.0842396, 1.8133969,   -1.3639419, -0.34885728, 1.2514706,  2.2248757,  1.398127,   1.6449664,   0.13868055,
    0.9444075,  1.5305994,   0.80286354, -1.0688899,  -1.1937977, 1.0168968,  1.2393475,  0.9242286,   -0.91741884,
    1.1878415,  -0.23807989, 0.7233633,  -2.145281,   -1.1852452, -3.1334703, 1.6308534,  -4.7782054,  1.4236984,
    -3.5798743, 2.5103958,   -1.3004229, 1.0359595,   1.3746276,  -1.2753546, -3.0728462, 0.9430954,   1.6196722,
    0.9992485,  -3.3693173,  0.9805231,  -0.3394561,  0.8614934,  -1.0240932, 1.0597516,  0.5085611,   -1.6215363,
    -1.3717827, -1.3680128,  0.43000492, -2.51508,    -5.9255157, -1.2641972 };
double I_G_concatenateBias1[I_GEO_CONCATENATE_DENSE_SIZE1 * 1]{ -1.2725621,  0.690456,    0.46663705, -0.6376881,
                                                                -0.66698545, -0.91276634, 0.8194415,  -0.46650434,
                                                                1.4056009,   -0.9653065 };
double I_G_concatenateWeight2[I_GEO_CONCATENATE_DENSE_SIZE2 * I_GEO_CONCATENATE_DENSE_SIZE1]{
    1.3823704,   -1.5274259, 1.2925509,  -1.4554502,  1.5435754,  -1.6804659, 1.7447766,  -1.6175056, 1.6671771,
    -1.7093256,  -1.8144429, 1.9928136,  -1.8260574,  2.0176878,  -2.0889704, 0.83046174, -1.0522404, 0.87358505,
    -0.96108586, 1.0603346,  1.8925407,  -1.9532624,  1.8821065,  -2.0153613, 1.9749352,  0.98415536, -0.937241,
    0.86796993,  -1.0128359, 1.1358864,  -2.0785148,  2.1628094,  -2.1083164, 2.2654207,  -2.3611178, 1.202473,
    -1.3255696,  1.2439586,  -1.2378174, 1.3642526,   -3.1626537, 3.392433,   -3.1374614, 3.3665216,  -3.5371337,
    0.8740252,   -0.9518637, 0.82592607, -0.92370886, 1.0701635 };
double I_G_concatenateBias2[I_GEO_CONCATENATE_DENSE_SIZE2 * 1]{ -0.3086615, 0.25321433, -0.30044767, 0.21724622,
                                                                -0.18322258 };
double I_G_concatenateOutputWeight[I_GEO_CONCATENATE_DENSE_SIZE2 * I_GEO_CONCATENATE_OUTPUT]{
    2.004059, -3.8482218, 2.0326178, -3.8882325, 2.1491752 };
double I_G_concatenateOutputBias[I_GEO_CONCATENATE_OUTPUT * 1]{ -0.8240644 };

// I Att
double I_A_convConvWeight1[I_ATT_CONV_KERNEL_SIZE1 * I_ATT_CONV_KERNEL_SIZE1]{
    -2.313751, -2.5657394, -3.5435393, -1.625864, -2.3104625, -2.9658153, -2.0316648, -1.6398494, -0.40192693 };
double I_A_convConvBias1[1]{ 0.6452012 };
double I_A_convConvWeight2[I_ATT_CONV_KERNEL_SIZE2 * I_ATT_CONV_KERNEL_SIZE2]{
    -8.84542, 1.0680593, 0.98437625, -8.097758, 0.42062566, 0.28371465, -7.4370794, 0.8762599, 0.8828099 };
double I_A_convConvBias2[1]{ -1.0178541 };
double I_A_concatenateWeight1[I_ATT_CONCATENATE_INPUT * I_ATT_CONCATENATE_DENSE_SIZE1]{
    1.4155316,   0.85760653,  1.0845468,   -0.39423513, -0.029985564, -0.714534,   -1.3492212,  -0.9887018, 2.150993,
    -0.23534842, -0.24626194, -0.34997857, -2.5576985,  -0.47634545,  -6.8082843,  0.9872692,   19.351383,  -0.834732,
    -1.2627342,  -0.9576514,  -0.18307836, -0.53646487, -0.064039014, 0.43902996,  -1.5825542,  0.6422879,  0.7573643,
    0.19846967,  0.2170225,   0.4754407,   -2.664847,   -0.9921201,   -0.99691856, 0.4219246,   -1.0807207, 1.9806353,
    1.3125556,   0.22260486,  7.4621468,   0.38077715,  -3.2867785,   -4.05741,    -0.95738095, 2.3385267,  2.0226688,
    3.277886,    -0.10141508, 4.974981,    -2.8144388,  2.8423483,    -0.37775132, -1.5271894,  0.9825234,  -1.2062823,
    0.40825826,  -0.7190107,  2.8415391,   3.6243372,   -0.14635897,  -1.2593734 };
double I_A_concatenateBias1[I_ATT_CONCATENATE_DENSE_SIZE1 * 1]{ 0.50649834, 0.045459267, 0.99696255,  -0.7059469,
                                                                0.16624726, -0.3276451,  -0.70123947, -0.49645388,
                                                                1.0248702,  -0.585051 };
double I_A_concatenateWeight2[I_ATT_CONCATENATE_DENSE_SIZE2 * I_ATT_CONCATENATE_DENSE_SIZE1]{
    -1.9213756, 2.0575862,  -1.8859855, -1.8630377, 1.8717602,  -1.7241359, 1.7612948,  -1.7566996, -1.808578,
    1.6362387,  -1.9870633, 2.227016,   -1.9141161, -2.0921948, 2.0285618,  0.8594412,  -1.0688462, 0.6294947,
    0.6794265,  -0.9672363, -1.7986736, 1.0391009,  -3.3266435, -3.3898597, 1.5129923,  0.95932764, -1.1734631,
    0.6726749,  0.6653629,  -1.0132568, 1.9808166,  -2.260547,  1.512272,   1.5542772,  -2.050884,  1.2245796,
    -1.2805103, 0.9252321,  0.9778746,  -1.2439013, -1.7717177, 2.0110571,  -1.9628607, -2.0220978, 1.8477818,
    0.9945002,  -1.2875917, 0.81730056, 0.7782702,  -1.0862621 };
double I_A_concatenateBias2[I_ATT_CONCATENATE_DENSE_SIZE2 * 1]{ -0.16705863, 0.06383626, -0.21596476, -0.112949155,
                                                                0.10939691 };
double I_A_concatenateOutputWeight[I_ATT_CONCATENATE_DENSE_SIZE2 * I_ATT_CONCATENATE_OUTPUT]{
    2.3603811, -4.5765676, 2.5549712, 2.5744483, -4.0866103 };
double I_A_concatenateOutputBias[I_ATT_CONCATENATE_OUTPUT * 1]{ -0.6012961 };

// LFCN
double P_G_weight1[P_G_hidden1_num * P_G_input_node]{
    -0.1798137,   0.015715983,  -0.09532983,  -0.12949441,  -0.016367447, 0.56853455,   -0.13140851,  -0.12590483,
    0.5730921,    0.07826559,   -0.119140275, 0.010980937,  -0.041487183, -0.015314037, -0.019670215, 0.4799034,
    -0.040534467, -0.050889056, 0.57576144,   -0.029060448, -0.5683938,   0.013555974,  0.042025037,  -0.7310104,
    -0.004070325, 1.1931382,    -0.6871277,   -0.026506832, 1.1419412,    -0.056911822 };
double P_G_bias1[P_G_hidden1_num]{ 0.8574694,  -0.03832592, -0.029819846, 0.8761555,   -0.05974367,
                                   -0.2603977, 0.86793953,  -0.04469238,  -0.18879512, -0.078819275 };
double P_G_weight2[P_G_hidden1_num * P_G_hidden2_num]{
    -1.5059104,  -1.3450825,  -1.43051,    1.477028,     1.480278,     -0.0066613876, -0.010042782, -0.10854612,
    0.043490637, 0.025427084, 0.010881705, -0.031549543, -0.010803781, -0.035900652,  -0.08626275,  -1.8950084,
    -1.8535957,  -1.8098618,  1.8292805,   1.7869675,    0.008852309,  -0.011014264,  0.023642553,  0.037875876,
    0.035256337, 1.2992209,   1.3258094,   1.350649,     -1.3084512,   -1.28103,      -1.8162,      -1.7824688,
    -1.7318273,  1.8208646,   1.7572726,   0.06479486,   0.057287093,  0.059295353,   0.07287708,   -0.05606584,
    1.1962177,   1.0893365,   1.1233077,   -1.1840847,   -1.1880049,   0.071519755,   0.0057525644, -0.06867385,
    -0.06903744, -0.036188617 };
double P_G_bias2[P_G_hidden2_num]{ -0.18987358, -0.08684765, -0.29997537, 0.23931256, 0.29157344 };
double P_G_weight3[P_G_hidden2_num * P_G_output_node]{ 2.3271866, 2.291772, 2.3157434, -2.4887452, -2.4943295 };
double P_G_bias3[P_G_output_node]{ -0.3431725 };

double P_A_weight1[P_A_hidden1_num * P_A_input_node]{
    0.38778436,   -0.007314995,  0.33233833,   0.36591738,  0.4208757,   0.3915192,   0.022857603, -0.015012168,
    0.0021270008, -0.0051038982, 0.059003677,  0.065315135, 0.07027739,  0.41112804,  0.0588191,   0.02490963,
    0.04832394,   0.008314324,   -0.03013271,  0.024819389, 0.92511517,  -0.09174364, 0.9075054,   0.8722612,
    1.0292256,    1.1075729,     -0.039412506, 0.039779317, 0.009898708, 0.036773242 };
double P_A_bias1[P_A_hidden1_num]{ -0.09741153, -0.05833369, -0.09079688,  -0.6771248,  -0.11397923,
                                   -0.6344252,  -0.06651615, -0.029567974, -0.05385367, -0.04727762 };
double P_A_weight2[P_A_hidden1_num * P_A_hidden2_num]{
    -0.9188661,   -1.0477335,   -0.9269158,   0.8647607,    1.0614995,   -0.08677276, 0.050005186, -0.0061682737,
    0.10102158,   -0.047413867, -0.8767638,   -1.1435426,   -0.9832114,  0.9889068,   1.1729834,   -2.9734714,
    -2.757544,    -2.9321315,   2.9902015,    2.7917132,    -1.0080203,  -1.1630646,  -1.0566517,  0.93173665,
    1.1421729,    -3.598869,    -3.1084049,   -3.3781476,   3.4908228,   3.0106144,   -0.09670318, -0.04485282,
    -0.022119856, -0.050046027, -0.060104396, -0.010417917, 0.03764891,  0.026927702, 0.016772207, -0.025491703,
    0.01552383,   -0.0738545,   0.030917257,  -0.0945312,   0.037042968, 0.004715443, -0.06587834, 0.0051364466,
    -0.041220848, -0.021615531 };
double P_A_bias2[P_A_hidden2_num]{ 1.7809409, 1.437488, 1.700708, -1.6271535, -1.4059961 };
double P_A_weight3[P_A_hidden2_num * P_A_output_node]{ -2.4459467, -2.511491, -2.4142091, 2.8143487, 2.7357085 };
double P_A_bias3[P_A_output_node]{ 0.48611394 };

double I_G_weight1[I_G_hidden1_num * I_G_input_node]{
    0.6891816,    -0.032851584, -0.07149959, 0.1327268,    -0.031111201, -0.013056578, -0.17099088, -0.13880853,
    -2.5727057,   -0.5805173,   0.13167034,  -0.054828633, -0.47254786,  0.6590843,    0.014007426, -0.032431256,
    -0.6763731,   -0.43721578,  -0.11417717, -0.59266067,  1.7472414,    0.04121032,   -0.75080377, 1.8364912,
    -0.101214945, -0.004063612, -0.5388431,  -0.76438826,  0.29453498,   -0.35285744 };
double I_G_bias1[I_G_hidden1_num]{ 0.053268004,  -0.15858771, 0.9373012, -0.0056827045, -0.060157437,
                                   -0.011107092, 1.0216832,   0.9162116, 0.7341579,     0.9352663 };
double I_G_weight2[I_G_hidden1_num * I_G_hidden2_num]{
    0.730095,     -1.5902374,   -1.4260284,  1.8777653,   2.1345243,      0.009160237,  0.0020384693, -0.024743749,
    0.029532894,  -0.025439087, -2.1620114,  1.9589972,   -0.50337994,    -2.2293231,   -2.6179252,   0.57881624,
    -1.7947954,   -0.67919695,  1.8029268,   2.536627,    -0.00031412483, -0.029971926, 0.057302352,  0.03557874,
    -0.014848113, 0.072738044,  -0.08819018, 0.013006314, -0.1172115,     0.086894706,  -3.095519,    1.5223482,
    -0.6347482,   -1.7568599,   -2.2616343,  -2.2669876,  2.0884976,      -0.42797062,  -2.1458473,   -2.5100088,
    -3.8781316,   1.2264118,    4.9437575,   -0.46441874, 1.3488111,      -3.7643235,   1.1565139,    -0.026370969,
    -1.3415911,   -1.7141935 };
double I_G_bias2[I_G_hidden2_num]{ -0.40222603, 0.038892575, -0.48762688, -0.06750394, 0.04712745 };
double I_G_weight3[I_G_hidden2_num * I_G_output_node]{ 4.585394, -4.635746, -5.3902903, 1.4180381, 1.6833978 };
double I_G_bias3[I_G_output_node]{ -1.2666923 };

double I_A_weight1[I_A_hidden1_num * I_A_input_node]{
    -0.3399334,  -4.2346177,  -0.09606053, -1.9594648,  -0.083944865, 0.7173832,    0.5212202,   -0.029234419,
    -0.42640632, 0.035182115, -0.28761864, -0.2260642,  -0.030421784, -0.1238716,   -0.05297501, 0.84043366,
    0.7088765,   -0.05069234, -0.20634308, -0.04941558, -1.027468,    0.0034996956, 0.015534068, 1.1255032,
    0.024796126, 1.6061252,   1.6438757,   0.08391107,  -1.0014163,   -0.10561229 };
double I_A_bias1[I_A_hidden1_num]{ 1.0356098,  0.4776824,  -0.1347511,  0.11729125, -0.052083444,
                                   -0.2584819, -0.1723621, -0.09161258, 1.0184947,  -0.03605778 };
double I_A_weight2[I_A_hidden1_num * I_A_hidden2_num]{
    -2.437292,    1.610624,   -1.5597881,  1.6304326,   -1.4551742,   1.157308,       6.680512,     -2.5165,
    6.8099475,    -5.1118183, 0.043091673, 0.05795094,  -0.016092842, 0.0010595184,   0.042020276,  -3.456434,
    0.05806105,   -2.4778214, -0.06356614, -1.4688472,  -0.09521533,  -0.00036730454, -0.089662425, -0.059034307,
    -0.010170054, 0.8435426,  -1.6657287,  1.1937068,   -1.6841154,   1.5154018,      0.6038032,    -1.8238503,
    0.95434284,   -1.814911,  1.2709001,   0.044319823, 0.031158386,  0.038542103,    -0.024861578, -0.06378681,
    -2.3413324,   1.5640033,  -1.2784584,  1.5854845,   -1.3741953,   -0.022077857,   0.013171807,  -0.014771372,
    0.022822618,  0.03274649 };
double I_A_bias2[I_A_hidden2_num]{ -0.88686234, 0.00551265, -0.46492913, -0.017391888, -0.17816842 };
double I_A_weight3[I_A_hidden2_num * I_A_output_node]{ 2.8143318, -3.5012615, 1.5834863, -3.524372, 1.5158705 };
double I_A_bias3[I_A_output_node]{ -0.7981659 };


double softmax( double* y, int cateNum ) {
  double sum   = 0;
  double p_max = 0, cate = -1;
  for ( int i = 0; i < cateNum; i++ ) sum += exp( y[i] );
  for ( int i = 0; i < cateNum; i++ ) {
    double p = exp( y[i] ) / sum;
    if ( p > p_max ) {
      p_max = p;
      cate  = i;
    }
  }
  return cate;
}
double optimizer( double u, string activation ) {
  if ( activation == "sigmoid" )
    return 1.0 / ( 1.0 + exp( -u ) );  // sigmoid
  else if ( activation == "relu" )     // relu
    if ( 0.0 >= u )
      return 0.0;
    else
      return u;
  else
    return 0.0;
}
double vanillayNN( string whichMode, int GorA, const vector<input_dtype>& x, string activation ) {
  int     i, j;
  double  u1;
  double  y0, y;
  int     hidden1_num, hidden2_num;
  int     input_node, output_node;
  double *weight_h1, *weight_h2, *bias_h1, *bias_h2, *weight_h3, *bias_h3;
  if ( whichMode == "PModule" ) {
    if ( GorA == 0 ) {
      input_node  = P_G_input_node;
      hidden1_num = P_G_hidden1_num;
      weight_h1   = P_G_weight1;
      weight_h2   = P_G_weight2;
      bias_h1     = P_G_bias1;
      bias_h2     = P_G_bias2;
      hidden2_num = P_G_hidden2_num;
      weight_h3   = P_G_weight3;
      bias_h3     = P_G_bias3;
    } else {
      input_node  = P_A_input_node;
      hidden1_num = P_A_hidden1_num;
      weight_h1   = P_A_weight1;
      weight_h2   = P_A_weight2;
      bias_h1     = P_A_bias1;
      bias_h2     = P_A_bias2;
      hidden2_num = P_A_hidden2_num;
      weight_h3   = P_A_weight3;
      bias_h3     = P_A_bias3;
    }
    double *x_hidden1, *x_hidden2;
    x_hidden1 = new double[hidden1_num];
    x_hidden2 = new double[hidden2_num];
    // layer1: input layer --> hidden layer1
    for ( i = 0; i < hidden1_num; i++ ) {
      u1 = 0;
      for ( j = 0; j < input_node; j++ ) u1 += (double)x[j] * weight_h1[j * hidden1_num + i];
      u1 += bias_h1[i];  // bias of layer 1
      x_hidden1[i] = optimizer( u1, "relu" );
    }
    // layer2: hidden layer1 --> hidden layer2
    for ( i = 0; i < hidden2_num; i++ ) {
      u1 = 0;
      for ( j = 0; j < hidden1_num; j++ ) u1 += (double)x_hidden1[j] * weight_h2[j * hidden2_num + i];
      u1 += bias_h2[i];  // bias of layer 1
      x_hidden2[i] = optimizer( u1, activation );
    }
    // layer3: hidden layer3 --> output layer
    y = 0;
    for ( i = 0; i < hidden2_num; ++i ) y += x_hidden2[i] * weight_h3[i];
    y += bias_h3[0];  // bias of layer 2
    delete[] x_hidden2;
    y0 = optimizer( y, "sigmoid" );
    delete[] x_hidden1;
  } 
  else if ( whichMode == "IModule" ) {
    if ( GorA == 0 ) {
      input_node  = I_G_input_node;
      hidden1_num = I_G_hidden1_num;
      weight_h1   = I_G_weight1;
      weight_h2   = I_G_weight2;
      bias_h1     = I_G_bias1;
      bias_h2     = I_G_bias2;
      hidden2_num = I_G_hidden2_num;
      weight_h3   = I_G_weight3;
      bias_h3     = I_G_bias3;
    } else {
      input_node  = I_A_input_node;
      hidden1_num = I_A_hidden1_num;
      weight_h1   = I_A_weight1;
      weight_h2   = I_A_weight2;
      bias_h1     = I_A_bias1;
      bias_h2     = I_A_bias2;
      hidden2_num = I_A_hidden2_num;
      weight_h3   = I_A_weight3;
      bias_h3     = I_A_bias3;
    }
    double *x_hidden1, *x_hidden2;
    x_hidden1 = new double[hidden1_num];
    x_hidden2 = new double[hidden2_num];
    // layer1: input layer --> hidden layer1
    for ( i = 0; i < hidden1_num; i++ ) {
      u1 = 0;
      for ( j = 0; j < input_node; j++ ) u1 += (double)x[j] * weight_h1[j * hidden1_num + i];
      u1 += bias_h1[i];  // bias of layer 1
      x_hidden1[i] = optimizer( u1, "relu" );
    }
    // layer2: hidden layer1 --> hidden layer2
    for ( i = 0; i < hidden2_num; i++ ) {
      u1 = 0;
      for ( j = 0; j < hidden1_num; j++ ) u1 += (double)x_hidden1[j] * weight_h2[j * hidden2_num + i];
      u1 += bias_h2[i];  // bias of layer 1
      x_hidden2[i] = optimizer( u1, activation );
    }
    // layer3: hidden layer3 --> output layer
    y = 0;
    for ( i = 0; i < hidden2_num; ++i ) y += x_hidden2[i] * weight_h3[i];
    y += bias_h3[0];  // bias of layer 2
    delete[] x_hidden2;
    y0 = optimizer( y, "sigmoid" );
    delete[] x_hidden1;
  }
  return y0;
}

// Same conv layer
double** convLayer( double** inputMatrix,
                    int      inputSize,
                    double*  kernel,
                    int      kernelSize,
                    double*  bias,
                    string   activation ) {
  double** outputMatrix = new double*[inputSize];
  for ( int i = 0; i < inputSize; i++ ) outputMatrix[i] = new double[inputSize];
  int      paddingSize  = kernelSize / 2;
  double**    matrixPadding = new double*[inputSize + 2 * paddingSize];
  for ( int i = 0; i < inputSize + 2 * paddingSize; i++ ) matrixPadding[i] = new double[inputSize + 2 * paddingSize];
  for ( int i = 0; i < inputSize + 2 * paddingSize; i++ ) {
    for ( int j = 0; j < inputSize + 2 * paddingSize; j++ ) {
      if ( i < paddingSize || i >= inputSize + paddingSize || j < paddingSize || j >= inputSize + paddingSize )
        matrixPadding[i][j] = 0;
      else
        matrixPadding[i][j] = inputMatrix[i - paddingSize][j - paddingSize];
    }
  }
  for ( int i = paddingSize; i < inputSize + paddingSize; i++ ) {
    for ( int j = paddingSize; j < inputSize + paddingSize; j++ ) {
      double sum = 0;
      for ( int kerneli = 0; kerneli < kernelSize; kerneli++ )
        for ( int kernelj = 0; kernelj < kernelSize; kernelj++ )
          sum += matrixPadding[i - paddingSize + kerneli][j - paddingSize + kernelj] * kernel[kerneli*kernelSize+kernelj];
      sum += bias[0];
      outputMatrix[i - paddingSize][j - paddingSize] = optimizer( sum, activation );
    }
  }
  //cout << "conv output: " << endl;
  //for ( int i = 0; i < inputSize; i++ )
  //  for ( int j = 0; j < inputSize; j++ ) cout << outputMatrix[i][j] << ",";
  //cout << endl;
  return outputMatrix;
}

// Pooling layer.
double** poolingLayer( double** inputMatrix, int inputSize, int poolingSize ){
  int      outputSize   = inputSize / poolingSize;
  double** outputMatrix = new double*[outputSize];
  for ( int i = 0; i < outputSize; i++ ) outputMatrix[i] = new double[outputSize];
  for ( int i = 0; i < inputSize; i += poolingSize ) {
    for ( int j = 0; j < inputSize; j += poolingSize ) {
      double maxEle = -1;
      for ( int poolingi = 0; poolingi < poolingSize; poolingi++ )
        for ( int poolingj = 0; poolingj < poolingSize; poolingj++ )
          if ( maxEle < inputMatrix[i + poolingi][j + poolingj] ) maxEle = inputMatrix[i + poolingi][j + poolingj];
      outputMatrix[int(i / poolingSize)][int(j / poolingSize)] = maxEle;
    }
  }
  //cout << "pooling output: "<<endl;
  //for ( int i = 0; i < outputSize; i++ )
  //  for ( int j = 0; j < outputSize; j++ ) cout << outputMatrix[i][j] << ",";
  //cout << endl;
  return outputMatrix;
}


double SLPADDCNN( double** inputMatrix, const vector<input_dtype>& x, int GorA, string whichModule ) {
  int convInputSize = 16, convKernelSize1 = 3, convOutputSize1 = 16, convPoolingSize1 = 4, 
                          convKernelSize2 = 2, convOutputSize2 = 4, convPoolingSize2 = 4;
  double **convOutput1;
  double *convWeight1, *convBias1, *convWeight2, *convBias2, *convWeight3, *convBias3, *convWeight4, *convBias4;

  int concatenateInputSize = 1 + 5, concatenateDenseSize1 = 10, concatenateDenseSize2 = 5, concatenateOutputSize = 1;
  double *concatenateInput, *concatenateDenseWeight1, *concatenateDenseBias1, *concatenateDense1Output, 
                            *concatenateDenseWeight2, *concatenateDenseBias2, *concatenateDense2Output, 
                            *concatenateOutputWeight, *concatenateOutputBias, *concatenateOutput;
  double y = 0;

  if ( whichModule == "PModule" ) {
    convInputSize         = ( GorA == 0 ) ? P_GEO_CONV_INPUT_SIZE : P_ATT_CONV_INPUT_SIZE;
    convKernelSize1        = ( GorA == 0 ) ? P_GEO_CONV_KERNEL_SIZE1 : P_ATT_CONV_KERNEL_SIZE1;
    convWeight1        = ( GorA == 0 ) ? P_G_convConvWeight1 : P_A_convConvWeight1;
    convBias1          = ( GorA == 0 ) ? P_G_convConvBias1 : P_A_convConvBias1;
    convOutputSize1      = convInputSize;
    convPoolingSize1        = ( GorA == 0 ) ? P_GEO_CONV_POOL_SIZE1 : P_ATT_CONV_POOL_SIZE1;
    convKernelSize2         = ( GorA == 0 ) ? P_GEO_CONV_KERNEL_SIZE2 : P_ATT_CONV_KERNEL_SIZE2;
    convWeight2             = ( GorA == 0 ) ? P_G_convConvWeight2 : P_A_convConvWeight2;
    convBias2               = ( GorA == 0 ) ? P_G_convConvBias2 : P_A_convConvBias2;
    convOutputSize2         = convOutputSize1 / convPoolingSize1;
    convPoolingSize2        = ( GorA == 0 ) ? P_GEO_CONV_POOL_SIZE2 : P_ATT_CONV_POOL_SIZE2;

    concatenateInputSize  = ( GorA == 0 ) ? P_GEO_CONCATENATE_INPUT : P_ATT_CONCATENATE_INPUT;
    concatenateDenseSize1  = ( GorA == 0 ) ? P_GEO_CONCATENATE_DENSE_SIZE1 : P_ATT_CONCATENATE_DENSE_SIZE1;
    concatenateDenseSize2   = ( GorA == 0 ) ? P_GEO_CONCATENATE_DENSE_SIZE2 : P_ATT_CONCATENATE_DENSE_SIZE2;
    concatenateDenseWeight1 = ( GorA == 0 ) ? P_G_concatenateWeight1 : P_A_concatenateWeight1;
    concatenateDenseBias1   = ( GorA == 0 ) ? P_G_concatenateBias1 : P_A_concatenateBias1;
    concatenateDenseWeight2 = ( GorA == 0 ) ? P_G_concatenateWeight2 : P_A_concatenateWeight2;
    concatenateDenseBias2   = ( GorA == 0 ) ? P_G_concatenateBias2 : P_A_concatenateBias2;
    concatenateOutputWeight = ( GorA == 0 ) ? P_G_concatenateOutputWeight : P_A_concatenateOutputWeight;
    concatenateOutputBias   = ( GorA == 0 ) ? P_G_concatenateOutputBias : P_A_concatenateOutputBias;
  } else if ( whichModule == "IModule" ) {
    convInputSize      = ( GorA == 0 ) ? I_GEO_CONV_INPUT_SIZE : I_ATT_CONV_INPUT_SIZE;
    convKernelSize1    = ( GorA == 0 ) ? I_GEO_CONV_KERNEL_SIZE1 : I_ATT_CONV_KERNEL_SIZE1;
    convWeight1        = ( GorA == 0 ) ? I_G_convConvWeight1 : I_A_convConvWeight1;
    convBias1          = ( GorA == 0 ) ? I_G_convConvBias1 : I_A_convConvBias1;
    convOutputSize1    = convInputSize;
    convPoolingSize1   = ( GorA == 0 ) ? I_GEO_CONV_POOL_SIZE1 : I_ATT_CONV_POOL_SIZE1;
    convKernelSize2    = ( GorA == 0 ) ? I_GEO_CONV_KERNEL_SIZE2 : I_ATT_CONV_KERNEL_SIZE2;
    convWeight2        = ( GorA == 0 ) ? I_G_convConvWeight2 : I_A_convConvWeight2;
    convBias2          = ( GorA == 0 ) ? I_G_convConvBias2 : I_A_convConvBias2;
    convOutputSize2    = convOutputSize1 / convPoolingSize1;
    convPoolingSize2   = ( GorA == 0 ) ? I_GEO_CONV_POOL_SIZE2 : I_ATT_CONV_POOL_SIZE2;

    concatenateInputSize    = ( GorA == 0 ) ? I_GEO_CONCATENATE_INPUT : I_ATT_CONCATENATE_INPUT;
    concatenateDenseSize1   = ( GorA == 0 ) ? I_GEO_CONCATENATE_DENSE_SIZE1 : I_ATT_CONCATENATE_DENSE_SIZE1;
    concatenateDenseSize2   = ( GorA == 0 ) ? I_GEO_CONCATENATE_DENSE_SIZE2 : I_ATT_CONCATENATE_DENSE_SIZE2;
    concatenateDenseWeight1 = ( GorA == 0 ) ? I_G_concatenateWeight1 : I_A_concatenateWeight1;
    concatenateDenseBias1   = ( GorA == 0 ) ? I_G_concatenateBias1 : I_A_concatenateBias1;
    concatenateDenseWeight2 = ( GorA == 0 ) ? I_G_concatenateWeight2 : I_A_concatenateWeight2;
    concatenateDenseBias2   = ( GorA == 0 ) ? I_G_concatenateBias2 : I_A_concatenateBias2;
    concatenateOutputWeight = ( GorA == 0 ) ? I_G_concatenateOutputWeight : I_A_concatenateOutputWeight;
    concatenateOutputBias   = ( GorA == 0 ) ? I_G_concatenateOutputBias : I_A_concatenateOutputBias;
  }

  // initialization
  convOutput1 = new double*[1];
  convOutput1[0] = new double[1];
  
  convOutput1 = poolingLayer(
      convLayer(
          poolingLayer( convLayer( inputMatrix, convInputSize, convWeight1, convKernelSize1, convBias1, "sigmoid" ),
                        convOutputSize1, convPoolingSize1 ),
          convOutputSize1 / convPoolingSize1, convWeight2, convKernelSize2, convBias2, "sigmoid" ),
      convOutputSize2, convPoolingSize2 );

  // concatenate.
  concatenateInput       = new double[concatenateInputSize];
  concatenateDense1Output = new double[concatenateDenseSize1];
  concatenateDense2Output = new double[concatenateDenseSize2];
  concatenateOutput      = new double[concatenateOutputSize];
  for ( int i = 0; i < concatenateInputSize; i++ ) {
    if ( i == 0)
      concatenateInput[i] = convOutput1[0][0];
    else
      concatenateInput[i] = x[i - 1];
  }
  for ( int i = 0; i < concatenateDenseSize1; i++ ) {
    double sum = 0;
    for ( int j = 0; j < concatenateInputSize; j++ )
      sum += concatenateInput[j] * concatenateDenseWeight1[j * concatenateDenseSize1 + i];
    sum += concatenateDenseBias1[i];
    concatenateDense1Output[i] = optimizer( sum, "sigmoid" );
  }
  for ( int i = 0; i < concatenateDenseSize2; i++ ) {
    double sum = 0;
    for ( int j = 0; j < concatenateDenseSize1; j++ )
      sum += concatenateDense1Output[j] * concatenateDenseWeight2[j * concatenateDenseSize2 + i];
    sum += concatenateDenseBias2[i];
    concatenateDense2Output[i] = optimizer( sum, "sigmoid" );
  }
  for ( int i = 0; i < concatenateOutputSize; i++ ) {
    double sum = 0;
    for ( int j = 0; j < concatenateDenseSize2; j++ )
      sum += concatenateDense2Output[j] * concatenateOutputWeight[j * concatenateOutputSize + i];
    sum += concatenateOutputBias[i];
    concatenateOutput[i] = optimizer( sum, "sigmoid" );
  }

  for ( int i = 0; i < concatenateOutputSize; i++ ) y += concatenateOutput[i];
  return y;
}

double autoTH(double QP,int GorA, string IorP) { 
    double x = (QP / 51 - 1.0 / 2.0)*2.0; 
    //if ( GorA == 0 && IorP == "IModule" )
    //  return ( ( 1.0 / 10.0 * ( exp( -x ) - exp( x ) ) / ( exp( -x ) + exp( x ) ) ) + 9.0 / 10.0 );  // I Frame
    //else if ( GorA > 0 && IorP == "IModule" )
    //  return ( ( 1.0 / 10.0 * ( exp( -x ) - exp( x ) ) / ( exp( -x ) + exp( x ) ) ) + 9.3 / 10.0 );  // I Frame
    //else if ( GorA == 0 && IorP == "PModule" )
    //  return ( ( 1.0 / 10.0 * ( exp( -x ) - exp( x ) ) / ( exp( -x ) + exp( x ) ) ) + 9.3 / 10.0 );  // P Frame
    //else if ( GorA > 0 && IorP == "PModule" )
    //    return ((1.0 / 10.0 * ( exp( -x ) - exp( x ) ) / ( exp( -x ) + exp( x ) )) + 9.5 / 10.0);    // P Frame

    if ( IorP == "IModule" )
      return ( ( 1.0 / 5.0 * ( exp( -x ) - exp( x ) ) / ( exp( -x ) + exp( x ) ) ) + 4.0 / 10.0 );
    else if ( IorP == "PModule" )
      return ( ( 1.0 / 5.0 * ( exp( -x ) - exp( x ) ) / ( exp( -x ) + exp( x ) ) ) + 6.0 / 10.0 );
}

int LCU[64][64]   = { 0 };
int P_LCU[64][64] = { 0 };

double ave_8[8][8]     = { 0 };
double D_sum_8[8][8]   = { 0 };
double P_ave_8[8][8]   = { 0 };
double P_D_sum_8[8][8] = { 0 };

double ave_16[4][4]     = { 0 };
double D_sum_16[4][4]   = { 0 };
double P_ave_16[4][4]   = { 0 };
double P_D_sum_16[4][4] = { 0 };

double ave_32[2][2]     = { 0 };
double D_sum_32[2][2]   = { 0 };
double P_ave_32[2][2]   = { 0 };
double P_D_sum_32[2][2] = { 0 };

double cacAverage( int x_begin, int y_begin, int x_end, int y_end, int**& pixels ) {
  double sum = 0;
  for ( int y = y_begin; y < y_end; y++ ) {
    for ( int x = x_begin; x < x_end; x++ ) { sum += pixels[y][x]; }
  }
  return sum / ( ( (double)x_end - x_begin ) * ( (double)y_end - y_begin ) );
}

double cacVariance( int x_begin, int y_begin, int x_end, int y_end, int**& pixels, double average )
{
  double sum              = 0;
  int*   tmp_edge         = new int[y_end - y_begin];
  int    tmp_edge_pointer = 0;
  int    pixels_count     = 0;
  for ( int y = y_begin; y < y_end; y++ ) {
    for ( int x = x_begin; x < x_end; x++ ) { sum += pow( pixels[y][x] - average, 2 ); }
  }
  pixels_count = ( ( x_end - x_begin ) * ( y_end - y_begin ) );

  return sum / pixels_count;
}

double twoDecimalDouble( const double& dbNum ) {
  stringstream strCode;
  strCode << std::setiosflags( std::ios::fixed ) << std::setprecision( 2 ) << dbNum;
  double douCode;
  strCode >> douCode;
  return douCode;
}
#endif  // SDMTEST

namespace pcc_hm {
using namespace std;


//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

/**
 \param    uhTotalDepth  total number of allowable depth
 \param    uiMaxWidth    largest CU width
 \param    uiMaxHeight   largest CU height
 \param    chromaFormat  chroma format
 */
Void TEncCu::create(UChar uhTotalDepth, UInt uiMaxWidth, UInt uiMaxHeight, ChromaFormat chromaFormat, UInt paletteMaxSize, UInt paletteMaxPredSize )
{
  Int i;

  m_uhTotalDepth   = uhTotalDepth + 1;
  m_ppcBestCU      = new TComDataCU*[m_uhTotalDepth-1];
  m_ppcTempCU      = new TComDataCU*[m_uhTotalDepth-1];

  m_ppcPredYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcResiYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcRecoYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcPredYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcResiYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcRecoYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcOrigYuv     = new TComYuv*[m_uhTotalDepth-1];
  m_ppcNoCorrYuv   = new TComYuv*[m_uhTotalDepth-1];
#if PCC_RDO_EXT
  m_ppcOccupancyYuv = new TComYuv*[m_uhTotalDepth - 1];
#endif

  UInt uiNumPartitions;
  for( i=0 ; i<m_uhTotalDepth-1 ; i++)
  {
    uiNumPartitions = 1<<( ( m_uhTotalDepth - i - 1 )<<1 );
    UInt uiWidth  = uiMaxWidth  >> i;
    UInt uiHeight = uiMaxHeight >> i;

    m_ppcBestCU[i] = new TComDataCU; m_ppcBestCU[i]->create( chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1), paletteMaxSize, paletteMaxPredSize );
    m_ppcTempCU[i] = new TComDataCU; m_ppcTempCU[i]->create( chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1), paletteMaxSize, paletteMaxPredSize );

    m_ppcPredYuvBest[i] = new TComYuv; m_ppcPredYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcResiYuvBest[i] = new TComYuv; m_ppcResiYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcRecoYuvBest[i] = new TComYuv; m_ppcRecoYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);

    m_ppcPredYuvTemp[i] = new TComYuv; m_ppcPredYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcResiYuvTemp[i] = new TComYuv; m_ppcResiYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcRecoYuvTemp[i] = new TComYuv; m_ppcRecoYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);

    m_ppcOrigYuv    [i] = new TComYuv; m_ppcOrigYuv    [i]->create(uiWidth, uiHeight, chromaFormat);
    m_ppcNoCorrYuv  [i] = new TComYuv; m_ppcNoCorrYuv  [i]->create(uiWidth, uiHeight, chromaFormat);
#if PCC_RDO_EXT
    m_ppcOccupancyYuv[i] = new TComYuv; m_ppcOccupancyYuv[i]->create(uiWidth, uiHeight, chromaFormat);
#endif
  }

  m_bEncodeDQP                     = false;
  m_stillToCodeChromaQpOffsetFlag  = false;
  m_cuChromaQpOffsetIdxPlus1       = 0;
  m_bFastDeltaQP                   = false;

  // initialize partition order.
  UInt* piTmp = &g_auiZscanToRaster[0];
  initZscanToRaster( m_uhTotalDepth, 1, 0, piTmp);
  initRasterToZscan( uiMaxWidth, uiMaxHeight, m_uhTotalDepth );

  // initialize conversion matrix from partition index to pel
  initRasterToPelXY( uiMaxWidth, uiMaxHeight, m_uhTotalDepth );
}

Void TEncCu::destroy()
{
  Int i;

  for( i=0 ; i<m_uhTotalDepth-1 ; i++)
  {
    if(m_ppcBestCU[i])
    {
      m_ppcBestCU[i]->destroy();      delete m_ppcBestCU[i];      m_ppcBestCU[i] = NULL;
    }
    if(m_ppcTempCU[i])
    {
      m_ppcTempCU[i]->destroy();      delete m_ppcTempCU[i];      m_ppcTempCU[i] = NULL;
    }
    if(m_ppcPredYuvBest[i])
    {
      m_ppcPredYuvBest[i]->destroy(); delete m_ppcPredYuvBest[i]; m_ppcPredYuvBest[i] = NULL;
    }
    if(m_ppcResiYuvBest[i])
    {
      m_ppcResiYuvBest[i]->destroy(); delete m_ppcResiYuvBest[i]; m_ppcResiYuvBest[i] = NULL;
    }
    if(m_ppcRecoYuvBest[i])
    {
      m_ppcRecoYuvBest[i]->destroy(); delete m_ppcRecoYuvBest[i]; m_ppcRecoYuvBest[i] = NULL;
    }
    if(m_ppcPredYuvTemp[i])
    {
      m_ppcPredYuvTemp[i]->destroy(); delete m_ppcPredYuvTemp[i]; m_ppcPredYuvTemp[i] = NULL;
    }
    if(m_ppcResiYuvTemp[i])
    {
      m_ppcResiYuvTemp[i]->destroy(); delete m_ppcResiYuvTemp[i]; m_ppcResiYuvTemp[i] = NULL;
    }
    if(m_ppcRecoYuvTemp[i])
    {
      m_ppcRecoYuvTemp[i]->destroy(); delete m_ppcRecoYuvTemp[i]; m_ppcRecoYuvTemp[i] = NULL;
    }
    if(m_ppcOrigYuv[i])
    {
      m_ppcOrigYuv[i]->destroy();     delete m_ppcOrigYuv[i];     m_ppcOrigYuv[i] = NULL;
    }
    if(m_ppcNoCorrYuv[i])
    {
      m_ppcNoCorrYuv[i]->destroy();   delete m_ppcNoCorrYuv[i];   m_ppcNoCorrYuv[i] = NULL;
    }
#if PCC_RDO_EXT
    if (m_ppcOccupancyYuv[i])
    {
      m_ppcOccupancyYuv[i]->destroy(); delete m_ppcOccupancyYuv[i]; m_ppcOccupancyYuv[i] = NULL;
    }
#endif
  }
  if(m_ppcBestCU)
  {
    delete [] m_ppcBestCU;
    m_ppcBestCU = NULL;
  }
  if(m_ppcTempCU)
  {
    delete [] m_ppcTempCU;
    m_ppcTempCU = NULL;
  }

  if(m_ppcPredYuvBest)
  {
    delete [] m_ppcPredYuvBest;
    m_ppcPredYuvBest = NULL;
  }
  if(m_ppcResiYuvBest)
  {
    delete [] m_ppcResiYuvBest;
    m_ppcResiYuvBest = NULL;
  }
  if(m_ppcRecoYuvBest)
  {
    delete [] m_ppcRecoYuvBest;
    m_ppcRecoYuvBest = NULL;
  }
  if(m_ppcPredYuvTemp)
  {
    delete [] m_ppcPredYuvTemp;
    m_ppcPredYuvTemp = NULL;
  }
  if(m_ppcResiYuvTemp)
  {
    delete [] m_ppcResiYuvTemp;
    m_ppcResiYuvTemp = NULL;
  }
  if(m_ppcRecoYuvTemp)
  {
    delete [] m_ppcRecoYuvTemp;
    m_ppcRecoYuvTemp = NULL;
  }
  if(m_ppcOrigYuv)
  {
    delete [] m_ppcOrigYuv;
    m_ppcOrigYuv = NULL;
  }
  if(m_ppcNoCorrYuv)
  {
    delete [] m_ppcNoCorrYuv;
    m_ppcNoCorrYuv = NULL;
  }
#if PCC_RDO_EXT
  if (m_ppcOccupancyYuv)
  {
    delete[] m_ppcOccupancyYuv;
    m_ppcOccupancyYuv = NULL;
  }
#endif
}

/** \param    pcEncTop      pointer of encoder class
 */
Void TEncCu::init( TEncTop* pcEncTop )
{
  m_pcEncCfg           = pcEncTop;
  m_pcPredSearch       = pcEncTop->getPredSearch();
  m_pcTrQuant          = pcEncTop->getTrQuant();
  m_pcRdCost           = pcEncTop->getRdCost();

  m_pcEntropyCoder     = pcEncTop->getEntropyCoder();
  m_pcBinCABAC         = pcEncTop->getBinCABAC();

  m_pppcRDSbacCoder    = pcEncTop->getRDSbacCoder();
  m_pcRDGoOnSbacCoder  = pcEncTop->getRDGoOnSbacCoder();

  m_pcRateCtrl         = pcEncTop->getRateCtrl();
  m_lumaQPOffset       = 0;
  initLumaDeltaQpLUT();
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

/** 
 \param  pCtu pointer of CU data class
 */
Void TEncCu::compressCtu( TComDataCU* pCtu, UChar* lastPaletteSize, Pel lastPalette[][MAX_PALETTE_PRED_SIZE] )
{
  // initialize CU data
  m_ppcBestCU[0]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );
  m_ppcTempCU[0]->initCtu( pCtu->getPic(), pCtu->getCtuRsAddr() );

  for (UChar comp = 0; comp < (pCtu->getSlice()->getSPS()->getChromaFormatIdc() == CHROMA_400 ? 1 : 3); comp++)
  {
    m_ppcBestCU[0]->setLastPaletteInLcuSizeFinal(comp, lastPaletteSize[comp]);
    m_ppcTempCU[0]->setLastPaletteInLcuSizeFinal(comp, lastPaletteSize[comp]);
    for (UInt idx = 0; idx < pCtu->getSlice()->getSPS()->getSpsScreenExtension().getPaletteMaxPredSize(); idx++)
    {
      m_ppcBestCU[0]->setLastPaletteInLcuFinal(comp, lastPalette[comp][idx], idx);
      m_ppcTempCU[0]->setLastPaletteInLcuFinal(comp, lastPalette[comp][idx], idx);
    }
  }

  // analysis of CU
  DEBUG_STRING_NEW(sDebug)

  xCompressCU( m_ppcBestCU[0], m_ppcTempCU[0], 0 DEBUG_STRING_PASS_INTO(sDebug) );
  DEBUG_STRING_OUTPUT(std::cout, sDebug)

#if ADAPTIVE_QP_SELECTION
  if( m_pcEncCfg->getUseAdaptQpSelect() )
  {
    if(pCtu->getSlice()->getSliceType()!=I_SLICE) //IIII
    {
      xCtuCollectARLStats( pCtu );
    }
  }
#endif
}
/** \param  pCtu  pointer of CU data class
 */
Void TEncCu::encodeCtu ( TComDataCU* pCtu )
{
  if ( pCtu->getSlice()->getPPS()->getUseDQP() )
  {
    setdQPFlag(true);
  }

  if ( pCtu->getSlice()->getUseChromaQpAdj() )
  {
    setCodeChromaQpAdjFlag(true);
  }

  // Encode CU data
  xEncodeCU( pCtu, 0, 0 );
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

Void TEncCu::initLumaDeltaQpLUT()
{
  const LumaLevelToDeltaQPMapping &mapping=m_pcEncCfg->getLumaLevelToDeltaQPMapping();

  if ( !mapping.isEnabled() )
  {
    return;
  }

  // map the sparse LumaLevelToDeltaQPMapping.mapping to a fully populated linear table.

  Int         lastDeltaQPValue=0;
  std::size_t nextSparseIndex=0;
  for(Int index=0; index<LUMA_LEVEL_TO_DQP_LUT_MAXSIZE; index++)
  {
    while (nextSparseIndex < mapping.mapping.size() && index>=mapping.mapping[nextSparseIndex].first)
    {
      lastDeltaQPValue=mapping.mapping[nextSparseIndex].second;
      nextSparseIndex++;
    }
    m_lumaLevelToDeltaQPLUT[index]=lastDeltaQPValue;
  }
}

Int TEncCu::calculateLumaDQP(TComDataCU *pCU, const UInt absPartIdx, const TComYuv * pOrgYuv)
{
  const Pel *pY = pOrgYuv->getAddr(COMPONENT_Y, absPartIdx);
  const Int stride  = pOrgYuv->getStride(COMPONENT_Y);
  Int width = pCU->getWidth(absPartIdx);
  Int height = pCU->getHeight(absPartIdx);
  Double avg = 0;

  // limit the block by picture size
  const TComSPS* pSPS = pCU->getSlice()->getSPS();
  if ( pCU->getCUPelX() + width > pSPS->getPicWidthInLumaSamples() )
  {
    width = pSPS->getPicWidthInLumaSamples() - pCU->getCUPelX();
  }
  if ( pCU->getCUPelY() + height > pSPS->getPicHeightInLumaSamples() )
  {
    height = pSPS->getPicHeightInLumaSamples() - pCU->getCUPelY();
  }

  // Get QP offset derived from Luma level
  if ( m_pcEncCfg->getLumaLevelToDeltaQPMapping().mode == LUMALVL_TO_DQP_AVG_METHOD )
  {
    // Use avg method
    Int sum = 0;
    for (Int y = 0; y < height; y++)
    {
      for (Int x = 0; x < width; x++)
      {
        sum += pY[x];
      }
      pY += stride;
    }
    avg = (Double)sum/(width*height);
  }
  else
  {
    // Use maximum luma value
    Int maxVal = 0;
    for (Int y = 0; y < height; y++)
    {
      for (Int x = 0; x < width; x++)
      {
        if (pY[x] > maxVal)
        {
          maxVal = pY[x];
        }
      }
      pY += stride;
    }
    // use a percentage of the maxVal
    avg = (Double)maxVal * m_pcEncCfg->getLumaLevelToDeltaQPMapping().maxMethodWeight;
  }

  Int lumaIdx = Clip3<Int>(0, Int(LUMA_LEVEL_TO_DQP_LUT_MAXSIZE)-1, Int(avg+0.5) );
  Int QP = m_lumaLevelToDeltaQPLUT[lumaIdx];
  return QP;
}

//! Derive small set of test modes for AMP encoder speed-up
#if AMP_ENC_SPEEDUP
#if AMP_MRG
Void TEncCu::deriveTestModeAMP (TComDataCU *pcBestCU, PartSize eParentPartSize, Bool &bTestAMP_Hor, Bool &bTestAMP_Ver, Bool &bTestMergeAMP_Hor, Bool &bTestMergeAMP_Ver)
#else
Void TEncCu::deriveTestModeAMP (TComDataCU *pcBestCU, PartSize eParentPartSize, Bool &bTestAMP_Hor, Bool &bTestAMP_Ver)
#endif
{
  if ( pcBestCU->getPartitionSize(0) == SIZE_2NxN )
  {
    bTestAMP_Hor = true;
  }
  else if ( pcBestCU->getPartitionSize(0) == SIZE_Nx2N )
  {
    bTestAMP_Ver = true;
  }
  else if ( pcBestCU->getPartitionSize(0) == SIZE_2Nx2N && pcBestCU->getMergeFlag(0) == false && pcBestCU->isSkipped(0) == false )
  {
    bTestAMP_Hor = true;
    bTestAMP_Ver = true;
  }

#if AMP_MRG
  //! Utilizing the partition size of parent PU
  if ( eParentPartSize >= SIZE_2NxnU && eParentPartSize <= SIZE_nRx2N )
  {
    bTestMergeAMP_Hor = true;
    bTestMergeAMP_Ver = true;
  }

  if ( eParentPartSize == NUMBER_OF_PART_SIZES ) //! if parent is intra
  {
    if ( pcBestCU->getPartitionSize(0) == SIZE_2NxN )
    {
      bTestMergeAMP_Hor = true;
    }
    else if ( pcBestCU->getPartitionSize(0) == SIZE_Nx2N )
    {
      bTestMergeAMP_Ver = true;
    }
  }

  if ( pcBestCU->getPartitionSize(0) == SIZE_2Nx2N && pcBestCU->isSkipped(0) == false )
  {
    bTestMergeAMP_Hor = true;
    bTestMergeAMP_Ver = true;
  }

  if ( pcBestCU->getWidth(0) == 64 )
  {
    bTestAMP_Hor = false;
    bTestAMP_Ver = false;
  }
#else
  //! Utilizing the partition size of parent PU
  if ( eParentPartSize >= SIZE_2NxnU && eParentPartSize <= SIZE_nRx2N )
  {
    bTestAMP_Hor = true;
    bTestAMP_Ver = true;
  }

  if ( eParentPartSize == SIZE_2Nx2N )
  {
    bTestAMP_Hor = false;
    bTestAMP_Ver = false;
  }
#endif
}
#endif


static Intermediate_Int CalculateMinimumHVLumaActivity(TComDataCU *pcCU, const UInt uiAbsPartIdx, const TComYuv * const * const ppcOrigYuv);

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================
/** Compress a CU block recursively with enabling sub-CTU-level delta QP
 *  - for loop of QP value to compress the current CU with all possible QP
*/
#if AMP_ENC_SPEEDUP
Void TEncCu::xCompressCU( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth DEBUG_STRING_FN_DECLARE(sDebug_), PartSize eParentPartSize )
#else
Void TEncCu::xCompressCU( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth )
#endif
{
  TComMv lastIntraBCMv[2];
  for(Int i=0; i<2; i++)
  {
    lastIntraBCMv[i] = rpcBestCU->getLastIntraBCMv(i);
  }

  UChar lastPaletteSize[3];
  UInt numValidComp = rpcBestCU->getPic()->getNumberValidComponents();
  Pel*  lastPalette[MAX_NUM_COMPONENT];
  for (UInt ch = 0; ch < numValidComp; ch++)
  {
    lastPalette[ch] = (Pel*)xMalloc(Pel, rpcBestCU->getSlice()->getSPS()->getSpsScreenExtension().getPaletteMaxPredSize());
  }
  for (UInt ch = 0; ch < numValidComp; ch++)
  {
    lastPaletteSize[ch] = rpcBestCU->getLastPaletteInLcuSizeFinal(ch);
    for (UInt i = 0; i < rpcBestCU->getSlice()->getSPS()->getSpsScreenExtension().getPaletteMaxPredSize(); i++)
    {
      lastPalette[ch][i] = rpcBestCU->getLastPaletteInLcuFinal(ch, i);
    }
  }

  Bool isPerfectMatch = false;
  Bool terminateAllFurtherRDO = false;

  TComMv iMVCandList[4][10];
  memset( iMVCandList, 0, sizeof( TComMv )*4*10 );

  TComPic* pcPic = rpcBestCU->getPic();
  DEBUG_STRING_NEW(sDebug)
  const TComPPS &pps=*(rpcTempCU->getSlice()->getPPS());
  const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
  
  // These are only used if getFastDeltaQp() is true
  const UInt fastDeltaQPCuMaxSize    = Clip3(sps.getMaxCUHeight()>>sps.getLog2DiffMaxMinCodingBlockSize(), sps.getMaxCUHeight(), 32u);

  // get Original YUV data from picture
  m_ppcOrigYuv[uiDepth]->copyFromPicYuv( pcPic->getPicYuvOrg(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu() );
#if PCC_RDO_EXT
  m_ppcOccupancyYuv[uiDepth]->copyFromPicYuv(pcPic->getOccupancyMapYuv(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu());
#endif

  // variable for Cbf fast mode PU decision
  Bool    doNotBlockPu = true;
  Bool    earlyDetectionSkipMode = false;

  const UInt uiLPelX   = rpcBestCU->getCUPelX();
  const UInt uiRPelX   = uiLPelX + rpcBestCU->getWidth(0)  - 1;
  const UInt uiTPelY   = rpcBestCU->getCUPelY();
  const UInt uiBPelY   = uiTPelY + rpcBestCU->getHeight(0) - 1;
  const UInt uiWidth   = rpcBestCU->getWidth(0);

  Int iBaseQP = xComputeQP( rpcBestCU, uiDepth );
  Int iMinQP;
  Int iMaxQP;
  Bool isAddLowestQP = false;

  const UInt numberValidComponents = rpcBestCU->getPic()->getNumberValidComponents();

  if( uiDepth <= pps.getMaxCuDQPDepth() )
  {
    Int idQP = m_pcEncCfg->getMaxDeltaQP();
    iMinQP = Clip3( -sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP-idQP );
    iMaxQP = Clip3( -sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP+idQP );
  }
  else
  {
    iMinQP = rpcTempCU->getQP(0);
    iMaxQP = rpcTempCU->getQP(0);
  }

  if ( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() )
  {
    if ( uiDepth <= pps.getMaxCuDQPDepth() )
    {
      // keep using the same m_QP_LUMA_OFFSET in the same CTU
      m_lumaQPOffset = calculateLumaDQP(rpcTempCU, 0, m_ppcOrigYuv[uiDepth]);
    }
    iMinQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP - m_lumaQPOffset);
    iMaxQP = iMinQP; // force encode choose the modified QO
  }

  if ( m_pcEncCfg->getUseRateCtrl() )
  {
    iMinQP = m_pcRateCtrl->getRCQP();
    iMaxQP = m_pcRateCtrl->getRCQP();
  }

  // transquant-bypass (TQB) processing loop variable initialisation ---

  const Int lowestQP = iMinQP; // For TQB, use this QP which is the lowest non TQB QP tested (rather than QP'=0) - that way delta QPs are smaller, and TQB can be tested at all CU levels.

  if ( (pps.getTransquantBypassEnabledFlag()) )
  {
    isAddLowestQP = true; // mark that the first iteration is to cost TQB mode.
    iMinQP = iMinQP - 1;  // increase loop variable range by 1, to allow testing of TQB mode along with other QPs
    if ( m_pcEncCfg->getCUTransquantBypassFlagForceValue() )
    {
      iMaxQP = iMinQP;
    }
  }

  TComSlice * pcSlice = rpcTempCU->getPic()->getSlice(rpcTempCU->getPic()->getCurrSliceIdx());

  if(rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans())
  {
    m_ppcBestCU[uiDepth]->setTmpIntraBCRDCost(MAX_DOUBLE);
    m_ppcTempCU[uiDepth]->setTmpIntraBCRDCost(MAX_DOUBLE);
    m_ppcBestCU[uiDepth]->setIntraBCCSCEnabled(true);
    m_ppcTempCU[uiDepth]->setIntraBCCSCEnabled(true);

    m_ppcBestCU[uiDepth]->setTmpInterRDCost(MAX_DOUBLE);
    m_ppcTempCU[uiDepth]->setTmpInterRDCost(MAX_DOUBLE);
    m_ppcBestCU[uiDepth]->setInterCSCEnabled(true);
    m_ppcTempCU[uiDepth]->setInterCSCEnabled(true);
    setEnableIntraTUACT(uiDepth, pcSlice);
    setEnableIBCTUACT(uiDepth, pcSlice);
    setEnableInterTUACT(uiDepth, pcSlice);
  }

#ifdef SDMTEST           // MesksCode
  bool DONTSKIP = true;
  int  POC      = rpcBestCU->getSlice()->getPOC();
  int  QP = rpcBestCU->getQP( 0 );
  int  CUcate   = -1;

  double b_max_variance = 0;
  double b_min_variance = 999999;
  double q_max_variance = 0;
  double q_min_variance = 999999;

  double IGeoTH = 0.3;
  double IAttTH = 0.3;
  double PGeoTH = 0.6;
  double PAttTH = 0.6;

  if ( OorGorA >= 0 )
    CUcate = CUClassify( uiWidth, uiWidth, uiTPelY, uiLPelX, POC );  // =0 unoccupancy block£¬=1 fill block£¬=2 boundary block

  double LFCNSWITCH = rpcBestCU->getSlice()->getSliceType() == P_SLICE ? true : false;
  //LFCNSWITCH        = true;
#endif


  const Bool bBoundary = !( uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples() );

  if ( !bBoundary )
  {
    for (Int iQP=iMinQP; iQP<=iMaxQP; iQP++)
    {
      const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP);

      if (bIsLosslessMode)
      {
        iQP = lowestQP;
      }
      if ( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() && uiDepth <= pps.getMaxCuDQPDepth() )
      {
        getSliceEncoder()->updateLambda(pcSlice, iQP);
      }

      m_cuChromaQpOffsetIdxPlus1 = 0;
      if (pcSlice->getUseChromaQpAdj())
      {
        /* Pre-estimation of chroma QP based on input block activity may be performed
         * here, using for example m_ppcOrigYuv[uiDepth] */
        /* To exercise the current code, the index used for adjustment is based on
         * block position
         */
        Int lgMinCuSize = sps.getLog2MinCodingBlockSize() +
                          std::max<Int>(0, sps.getLog2DiffMaxMinCodingBlockSize()-Int(pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth()));
        m_cuChromaQpOffsetIdxPlus1 = ((uiLPelX >> lgMinCuSize) + (uiTPelY >> lgMinCuSize)) % (pps.getPpsRangeExtension().getChromaQpOffsetListLen() + 1);
      }

      rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

      // do inter modes, SKIP and 2Nx2N
      if ( ( !rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && rpcBestCU->getSlice()->getSliceType() != I_SLICE ) ||
           ( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && !rpcBestCU->getSlice()->isOnlyCurrentPictureAsReference() ) )
      {
        if ( m_pcEncCfg->getUseHashBasedME() )
        {
          xCheckRDCostHashInter( rpcBestCU, rpcTempCU, isPerfectMatch DEBUG_STRING_PASS_INTO(sDebug) );
          rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

          if ( isPerfectMatch )
          {
            if ( uiDepth == 0 )
            {
              terminateAllFurtherRDO = true;
            }
          }
        }

        // 2Nx2N
        if(m_pcEncCfg->getUseEarlySkipDetection() && !terminateAllFurtherRDO)
        {
          xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug) );
          rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );//by Competition for inter_2Nx2N
        }
        // SKIP
        xCheckRDCostMerge2Nx2N( rpcBestCU, rpcTempCU DEBUG_STRING_PASS_INTO(sDebug), &earlyDetectionSkipMode, terminateAllFurtherRDO );//by Merge for inter_2Nx2N
        rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

        if(!m_pcEncCfg->getUseEarlySkipDetection() && !terminateAllFurtherRDO)
        {
          // 2Nx2N, NxN
          xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug) );
          rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
          if(m_pcEncCfg->getUseCbfFastMode())
          {
            doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
          }
        }
      }

      if (bIsLosslessMode) // Restore loop variable if lossless mode was searched.
      {
        iQP = iMinQP;
      }
    }

    if(!earlyDetectionSkipMode && !terminateAllFurtherRDO)
    {
      for (Int iQP=iMinQP; iQP<=iMaxQP; iQP++)
      {
        const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP); // If lossless, then iQP is irrelevant for subsequent modules.

        if (bIsLosslessMode)
        {
          iQP = lowestQP;
        }

        rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

        // do inter modes, NxN, 2NxN, and Nx2N
        if ( ( !rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && rpcBestCU->getSlice()->getSliceType() != I_SLICE ) ||
             ( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && !rpcBestCU->getSlice()->isOnlyCurrentPictureAsReference() ) )
        {
          // 2Nx2N, NxN

          if(!( (rpcBestCU->getWidth(0)==8) && (rpcBestCU->getHeight(0)==8) ))
          {
            if( uiDepth == sps.getLog2DiffMaxMinCodingBlockSize() && doNotBlockPu)
            {
              xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug)   );
              rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            }
          }

          if(doNotBlockPu)
          {
            xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_Nx2N DEBUG_STRING_PASS_INTO(sDebug), false, iMVCandList[SIZE_Nx2N] );
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_Nx2N )
            {
              doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
            }
          }
          if(doNotBlockPu)
          {
            xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxN DEBUG_STRING_PASS_INTO(sDebug), false, iMVCandList[SIZE_2NxN] );
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxN)
            {
              doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
            }
          }

          //! Try AMP (SIZE_2NxnU, SIZE_2NxnD, SIZE_nLx2N, SIZE_nRx2N)
          if(sps.getUseAMP() && uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() )
          {
#if AMP_ENC_SPEEDUP
            Bool bTestAMP_Hor = false, bTestAMP_Ver = false;

#if AMP_MRG
            Bool bTestMergeAMP_Hor = false, bTestMergeAMP_Ver = false;

            deriveTestModeAMP (rpcBestCU, eParentPartSize, bTestAMP_Hor, bTestAMP_Ver, bTestMergeAMP_Hor, bTestMergeAMP_Ver);
#else
            deriveTestModeAMP (rpcBestCU, eParentPartSize, bTestAMP_Hor, bTestAMP_Ver);
#endif

            //! Do horizontal AMP
            if ( bTestAMP_Hor )
            {
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnU DEBUG_STRING_PASS_INTO(sDebug) );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnU )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnD DEBUG_STRING_PASS_INTO(sDebug) );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnD )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
            }
#if AMP_MRG
            else if ( bTestMergeAMP_Hor )
            {
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnU DEBUG_STRING_PASS_INTO(sDebug), true );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnU )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnD DEBUG_STRING_PASS_INTO(sDebug), true );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnD )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
            }
#endif

            //! Do horizontal AMP
            if ( bTestAMP_Ver )
            {
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nLx2N DEBUG_STRING_PASS_INTO(sDebug) );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_nLx2N )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nRx2N DEBUG_STRING_PASS_INTO(sDebug) );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
              }
            }
#if AMP_MRG
            else if ( bTestMergeAMP_Ver )
            {
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nLx2N DEBUG_STRING_PASS_INTO(sDebug), true );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_nLx2N )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nRx2N DEBUG_STRING_PASS_INTO(sDebug), true );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
              }
            }
#endif

#else
            xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnU );
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnD );
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nLx2N );
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

            xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nRx2N );
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

#endif
          }
        }

        // do normal intra modes
        // speedup for inter frames
        Double intraCost = MAX_DOUBLE;
        Double dIntraBcCostPred = 0.0;
        if ( ( !rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && rpcBestCU->getSlice()->getSliceType() == I_SLICE ) ||
             ( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && rpcBestCU->getSlice()->isOnlyCurrentPictureAsReference() ) ||
             !rpcBestCU->isSkipped(0) ) // avoid very complex intra if it is unlikely
        {
          if (m_pcEncCfg->getUseIntraBlockCopyFastSearch() && rpcTempCU->getWidth(0) <= SCM_S0067_MAX_CAND_SIZE )
          {
            xCheckRDCostIntraBC( rpcBestCU, rpcTempCU, false, SIZE_2Nx2N, dIntraBcCostPred, true DEBUG_STRING_PASS_INTO(sDebug));
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
          }

#if MCTS_ENC_CHECK
          if ( m_pcEncCfg->getTMCTSSEITileConstraint() || (rpcBestCU->getPredictionMode(0) == NUMBER_OF_PREDICTION_MODES) ||
              ((!m_pcEncCfg->getDisableIntraPUsInInterSlices()) && (
              (rpcBestCU->getCbf(0, COMPONENT_Y) != 0) ||
              ((rpcBestCU->getCbf(0, COMPONENT_Cb) != 0) && (numberValidComponents > COMPONENT_Cb)) ||
              ((rpcBestCU->getCbf(0, COMPONENT_Cr) != 0) && (numberValidComponents > COMPONENT_Cr))  // avoid very complex intra if it is unlikely
              )))
          {
#else
          if( rpcBestCU->getPredictionMode(0) == NUMBER_OF_PREDICTION_MODES ||
              ((!m_pcEncCfg->getDisableIntraPUsInInterSlices()) && (
               (rpcBestCU->getCbf( 0, COMPONENT_Y  ) != 0)                                            ||
              ((rpcBestCU->getCbf( 0, COMPONENT_Cb ) != 0) && (numberValidComponents > COMPONENT_Cb)) ||
              ((rpcBestCU->getCbf( 0, COMPONENT_Cr ) != 0) && (numberValidComponents > COMPONENT_Cr))  ) // avoid very complex intra if it is unlikely
            ))
          {
#endif 
            if(rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() && ( !bIsLosslessMode || (sps.getBitDepth( CHANNEL_TYPE_LUMA ) == sps.getBitDepth( CHANNEL_TYPE_CHROMA ))))
            {
              Double tempIntraCost = MAX_DOUBLE;
              TComDataCU *pStoredBestCU = rpcBestCU;
              if(m_pcEncCfg->getRGBFormatFlag())
              {
                xCheckRDCostIntraCSC( rpcBestCU, rpcTempCU, tempIntraCost, SIZE_2Nx2N, ACT_TRAN_CLR, false DEBUG_STRING_PASS_INTO(sDebug) );
              }
              else
              {
                if( uiDepth < 3 || bIsLosslessMode )
                {
                  xCheckRDCostIntra( rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
                }
                else
                {
                  xCheckRDCostIntraCSC( rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N, ACT_ORG_CLR, false DEBUG_STRING_PASS_INTO(sDebug) );
                }
              }
              Bool bSkipRGBCoding = pStoredBestCU == rpcBestCU ? !rpcTempCU->getQtRootCbf( 0 ): !rpcBestCU->getQtRootCbf( 0 );
              if(!bSkipRGBCoding)
              {
                rpcTempCU->initRQTData( uiDepth, rpcBestCU, (pStoredBestCU != rpcBestCU), false, true );
                if(m_pcEncCfg->getRGBFormatFlag())
                {
                  if( !getEnableIntraTUACT() )
                  {
                    xCheckRDCostIntra( rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug), true );
                  }
                  else
                  {
                    xCheckRDCostIntraCSC( rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N, ACT_TWO_CLR, true DEBUG_STRING_PASS_INTO(sDebug) );
                  }
                }
                else
                {
                  xCheckRDCostIntraCSC( rpcBestCU, rpcTempCU, tempIntraCost, SIZE_2Nx2N, ACT_TRAN_CLR, true DEBUG_STRING_PASS_INTO(sDebug) );
                }
                intraCost = std::min(intraCost, tempIntraCost);
              }
              else
                intraCost = tempIntraCost;
            }
            else
            {
              xCheckRDCostIntra( rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug) );
            }
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            if( uiDepth == sps.getLog2DiffMaxMinCodingBlockSize() )
            {
              if( rpcTempCU->getWidth(0) > ( 1 << sps.getQuadtreeTULog2MinSize() ) )
              {
                Double tmpIntraCost;
                if(rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() && ( !bIsLosslessMode || (sps.getBitDepth( CHANNEL_TYPE_LUMA ) == sps.getBitDepth( CHANNEL_TYPE_CHROMA ))))
                {
                  Double      tempIntraCost = MAX_DOUBLE;
                  TComDataCU *pStoredBestCU = rpcBestCU;
                  if(m_pcEncCfg->getRGBFormatFlag())
                  {
                    xCheckRDCostIntraCSC( rpcBestCU, rpcTempCU, tempIntraCost, SIZE_NxN, ACT_TRAN_CLR, false DEBUG_STRING_PASS_INTO(sDebug) );
                  }
                  else
                  {
                    if( uiDepth < 3 || bIsLosslessMode )
                    {
                      xCheckRDCostIntra( rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug) );
                    }
                    else
                    {
                      xCheckRDCostIntraCSC( rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN, ACT_ORG_CLR, false DEBUG_STRING_PASS_INTO(sDebug) );
                    }
                  }
                  Bool bSkipRGBCoding = pStoredBestCU == rpcBestCU ? !rpcTempCU->getQtRootCbf( 0 ) : !rpcBestCU->getQtRootCbf( 0 );
                  if(!bSkipRGBCoding)
                  {
                    rpcTempCU->initRQTData( uiDepth, rpcBestCU, (pStoredBestCU != rpcBestCU), false, true );
                    if(m_pcEncCfg->getRGBFormatFlag())
                    {
                      if( !getEnableIntraTUACT() )
                      {
                        xCheckRDCostIntra( rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug), true );
                      }
                      else
                      {
                        xCheckRDCostIntraCSC( rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN, ACT_TWO_CLR, true DEBUG_STRING_PASS_INTO(sDebug) );
                      }
                    }
                    else
                    {
                      xCheckRDCostIntraCSC( rpcBestCU, rpcTempCU, tempIntraCost, SIZE_NxN, ACT_TRAN_CLR, true DEBUG_STRING_PASS_INTO(sDebug) );
                    }
                    tmpIntraCost = std::min(tmpIntraCost, tempIntraCost);
                  }
                  else
                    tmpIntraCost = tempIntraCost;
                }
                else
                {
                  xCheckRDCostIntra( rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug)   );
                }
                intraCost = std::min(intraCost, tmpIntraCost);
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
              }
            }
          }

          if( rpcBestCU->isIntraBC( 0 ) )
          {
            intraCost = dIntraBcCostPred;
          }
        }

        intraCost = MAX_DOUBLE;

        // test PCM
        if(sps.getUsePCM()
          && rpcTempCU->getWidth(0) <= (1<<sps.getPCMLog2MaxSize())
          && rpcTempCU->getWidth(0) >= (1<<sps.getPCMLog2MinSize()) )
        {
          UInt uiRawBits = getTotalBits(rpcBestCU->getWidth(0), rpcBestCU->getHeight(0), rpcBestCU->getPic()->getChromaFormat(), sps.getBitDepths().recon);
          UInt uiBestBits = rpcBestCU->getTotalBits();
          if((uiBestBits > uiRawBits) || (rpcBestCU->getTotalCost() > m_pcRdCost->calcRdCost(uiRawBits, 0)))
          {
            xCheckIntraPCM (rpcBestCU, rpcTempCU);
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
          }
        }

        Bool bUse1DSearchFor8x8        = false;
        Bool bSkipIntraBlockCopySearch = false;

        if (rpcTempCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy())
        {
          xCheckRDCostIntraBCMerge2Nx2N( rpcBestCU, rpcTempCU );
          rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
        }

        if( !rpcBestCU->isSkipped(0) ) // avoid very complex intra if it is unlikely
        {
          if (m_pcEncCfg->getUseIntraBlockCopyFastSearch())
          {
            bSkipIntraBlockCopySearch = ((rpcTempCU->getWidth(0) > 16) || (intraCost < std::max(32*m_pcRdCost->getLambda(), 48.0)));

            if (rpcTempCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() &&
                !bSkipIntraBlockCopySearch &&
                rpcTempCU->getWidth(0) == 8 &&
                !m_ppcBestCU[uiDepth -1]->isIntraBC(0) )
            {
              bUse1DSearchFor8x8 = (CalculateMinimumHVLumaActivity(rpcTempCU, 0, m_ppcOrigYuv) < (168 << (sps.getBitDepth( CHANNEL_TYPE_LUMA ) - 8)));
            }
          }

          if (rpcTempCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy())
          {
            if (!bSkipIntraBlockCopySearch)
            {
              Double adIntraBcCost[NUMBER_OF_PART_SIZES] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
              xCheckRDCostIntraBC( rpcBestCU, rpcTempCU, bUse1DSearchFor8x8, SIZE_2Nx2N, adIntraBcCost[SIZE_2Nx2N] DEBUG_STRING_PASS_INTO(sDebug));
              rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

              if (m_pcEncCfg->getUseIntraBlockCopyFastSearch())
              {
                if( uiDepth == sps.getLog2DiffMaxMinCodingBlockSize() ) // Only additionally check Nx2N, 2NxN, NxN at the bottom level in fast search
                {
                  intraCost = std::min( intraCost, adIntraBcCost[SIZE_2Nx2N] );

                  Double dTH2 = std::max( 60 * m_pcRdCost->getLambda(),  56.0 );
                  Double dTH3 = std::max( 66 * m_pcRdCost->getLambda(), 800.0 );
                  if( intraCost >= dTH2 ) // only check Nx2N depending on best intraCost (and intraBCcost so far)
                  {
                    xCheckRDCostIntraBC( rpcBestCU, rpcTempCU, true, SIZE_Nx2N, adIntraBcCost[SIZE_Nx2N], false, (iMVCandList[SIZE_Nx2N]+8) DEBUG_STRING_PASS_INTO(sDebug));
                    rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                    intraCost = std::min( intraCost, adIntraBcCost[SIZE_Nx2N] );

                    if ( ( !rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && rpcBestCU->getSlice()->getSliceType() != I_SLICE ) ||
                         ( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && !rpcBestCU->getSlice()->isOnlyCurrentPictureAsReference() ) )
                    {
                      xCheckRDCostIntraBCMixed( rpcBestCU, rpcTempCU, SIZE_Nx2N, adIntraBcCost[SIZE_Nx2N] DEBUG_STRING_PASS_INTO(sDebug), iMVCandList[SIZE_Nx2N]);
                      rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                      intraCost = std::min( intraCost, adIntraBcCost[SIZE_Nx2N] );
                    }
                  }

                  if( intraCost >= dTH2 && !bIsLosslessMode ) // only check 2NxN depending on best intraCost (and intraBCcost so far) and if it is not lossless
                  {
                    xCheckRDCostIntraBC( rpcBestCU, rpcTempCU, ( bUse1DSearchFor8x8 || intraCost < dTH3 ), SIZE_2NxN, adIntraBcCost[SIZE_2NxN], false, (iMVCandList[SIZE_2NxN]+8) DEBUG_STRING_PASS_INTO(sDebug));
                    rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                    intraCost = std::min( intraCost, adIntraBcCost[SIZE_2NxN] );

                    if ( ( !rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && rpcBestCU->getSlice()->getSliceType() != I_SLICE ) ||
                         ( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseIntraBlockCopy() && !rpcBestCU->getSlice()->isOnlyCurrentPictureAsReference() ) )
                    {
                      xCheckRDCostIntraBCMixed( rpcBestCU, rpcTempCU, SIZE_2NxN, adIntraBcCost[SIZE_2NxN] DEBUG_STRING_PASS_INTO(sDebug), iMVCandList[SIZE_2NxN]);
                      rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                      intraCost = std::min( intraCost, adIntraBcCost[SIZE_2NxN] );
                    }
                }
              }
            }
            else
            {
              // full search (bUse1DSearchFor8x8 will be false but is kept here for consistency).

              xCheckRDCostIntraBC( rpcBestCU, rpcTempCU, bUse1DSearchFor8x8, SIZE_Nx2N, adIntraBcCost[SIZE_Nx2N] DEBUG_STRING_PASS_INTO(sDebug));
              rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
              xCheckRDCostIntraBC( rpcBestCU, rpcTempCU, bUse1DSearchFor8x8, SIZE_2NxN, adIntraBcCost[SIZE_2NxN] DEBUG_STRING_PASS_INTO(sDebug));
              rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            }
          }
        }

          if (bIsLosslessMode) // Restore loop variable if lossless mode was searched.
          {
            iQP = iMinQP;
          }
          if ( rpcBestCU->getSlice()->getSPS()->getSpsScreenExtension().getUsePaletteMode() )
          {
            //Change Palette QP dependent error limit
            Int iQP_Palette=Int(rpcBestCU->getQP(0));
            Int iQPrem = iQP_Palette % 6;
            Int iQPper = iQP_Palette / 6;
            Double quantiserScale = g_quantScales[iQPrem];
            Int quantiserRightShift = QUANT_SHIFT + iQPper;

            Double dQP=((Double)(1<<quantiserRightShift))/quantiserScale;

            UInt paletteQP;
            paletteQP=(UInt)(2.0*dQP/3.0+0.5);
            m_pcPredSearch->setPaletteErrLimit(paletteQP);

            Bool forcePalettePrediction = false;
            for( UChar ch = 0; ch < numValidComp; ch++ )
            {
              forcePalettePrediction = forcePalettePrediction || ( rpcTempCU->getLastPaletteInLcuSizeFinal( ch ) > 0 );
            }

            UInt iterNumber=0, paletteSize[2] = {MAX_PALETTE_SIZE, MAX_PALETTE_SIZE}, testedModes[4];

            if( rpcTempCU->getWidth(0) != 64)
            {
              iterNumber = 0;
              testedModes[iterNumber]=xCheckPaletteMode( rpcBestCU, rpcTempCU, false, iterNumber, paletteSize);
              rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

              if( !bIsLosslessMode )
              {
                if (paletteSize[0]>2 && testedModes[0]>0)
                {
                  iterNumber = 2;
                  testedModes[iterNumber]=xCheckPaletteMode( rpcBestCU, rpcTempCU, false, iterNumber, paletteSize);
                  rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                }

                if( forcePalettePrediction)
                {
                  iterNumber = 1;
                  testedModes[iterNumber]=xCheckPaletteMode( rpcBestCU, rpcTempCU, true, iterNumber, paletteSize+1);
                  rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                }

                if (forcePalettePrediction && paletteSize[1]>2 && testedModes[1]>0)
                {
                  iterNumber = 3;
                  testedModes[iterNumber]=xCheckPaletteMode( rpcBestCU, rpcTempCU, true, iterNumber, paletteSize+1);
                  rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                }
              }
            }
          }
        }
      }
    }

    // If Intra BC keep last coded Mv
    if(rpcBestCU->isInter(0))
    {
      Int iRefIdxFirst = rpcBestCU->getCUMvField( REF_PIC_LIST_0 )->getRefIdx( 0 );
      Int iRefIdxLast = rpcBestCU->getCUMvField( REF_PIC_LIST_0 )->getRefIdx( rpcBestCU->getTotalNumPart() - 1 );
      Bool isIntraBCFirst = ( iRefIdxFirst >= 0 ) ? rpcBestCU->getSlice()->getRefPic( REF_PIC_LIST_0, iRefIdxFirst )->getPOC() == rpcBestCU->getSlice()->getPOC() : false;
      Bool isIntraBCLast = ( iRefIdxLast >= 0 ) ? rpcBestCU->getSlice()->getRefPic( REF_PIC_LIST_0, iRefIdxLast )->getPOC() == rpcBestCU->getSlice()->getPOC() : false;

      if (  isIntraBCFirst || isIntraBCLast  )
      {
        if( rpcBestCU->getPartitionSize( 0 ) == SIZE_2Nx2N)
        {
          if( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( rpcBestCU->getTotalNumPart() - 1 ) != rpcBestCU->getLastIntraBCMv(0))
          {
            rpcBestCU->setLastIntraBCMv( rpcBestCU->getLastIntraBCMv(0), 1 );
            rpcBestCU->setLastIntraBCMv( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( rpcBestCU->getTotalNumPart() - 1 ) );
          }
        }
        else if(rpcBestCU->getPartitionSize( 0 ) == SIZE_2NxN || rpcBestCU->getPartitionSize( 0 ) == SIZE_Nx2N)
        {
          // mixed PU, only one partition is IntraBC coded
          if( isIntraBCFirst != isIntraBCLast )
          {
            if(isIntraBCFirst)
            {
              // Part 0
              if( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( 0 ) != rpcBestCU->getLastIntraBCMv())
              {
                rpcBestCU->setLastIntraBCMv( rpcBestCU->getLastIntraBCMv(), 1);
                rpcBestCU->setLastIntraBCMv( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( 0 ) );
              }
            }
            else if(isIntraBCLast)
            {
              // Part 1
              if( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( rpcBestCU->getTotalNumPart() - 1 ) != rpcBestCU->getLastIntraBCMv())
              {
                rpcBestCU->setLastIntraBCMv( rpcBestCU->getLastIntraBCMv(), 1);
                rpcBestCU->setLastIntraBCMv( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( rpcBestCU->getTotalNumPart() - 1 ) );
              }
            }
          }
          else // normal IntraBC CU
          {
            // Part 0
            if( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( 0 ) != rpcBestCU->getLastIntraBCMv())
            {
              rpcBestCU->setLastIntraBCMv( rpcBestCU->getLastIntraBCMv(), 1);
              rpcBestCU->setLastIntraBCMv( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( 0 ) );
            }
            // Part 1
            if( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( rpcBestCU->getTotalNumPart() - 1 ) != rpcBestCU->getLastIntraBCMv())
            {
              rpcBestCU->setLastIntraBCMv( rpcBestCU->getLastIntraBCMv(), 1);
              rpcBestCU->setLastIntraBCMv( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( rpcBestCU->getTotalNumPart() - 1 ) );
            }
          }
        }
        else
        {
          // NxN
          for(Int part=0; part<4; part++)
          {
            if( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( rpcBestCU->getTotalNumPart() - 4 + part ) != rpcBestCU->getLastIntraBCMv())
            {
              rpcBestCU->setLastIntraBCMv( rpcBestCU->getLastIntraBCMv(), 1 );
              rpcBestCU->setLastIntraBCMv( rpcBestCU->getCUMvField(REF_PIC_LIST_0)->getMv( rpcBestCU->getTotalNumPart() - 4 + part ) );
            }
          }
        }
      }
    } // is inter
    if (rpcBestCU->getPaletteModeFlag(0))
    {
      rpcBestCU->saveLastPaletteInLcuFinal( rpcBestCU, 0, rpcBestCU->getSlice()->getSPS()->getChromaFormatIdc() == CHROMA_400 ? 1 : 3 );
    }

    if( rpcBestCU->getTotalCost()!=MAX_DOUBLE )
    {
      m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);
      m_pcEntropyCoder->resetBits();
      m_pcEntropyCoder->encodeSplitFlag( rpcBestCU, 0, uiDepth, true );
      rpcBestCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
      rpcBestCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
      rpcBestCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcBestCU->getTotalBits(), rpcBestCU->getTotalDistortion() );
      m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);
    }
  }

  // copy original YUV samples to PCM buffer
  if( rpcBestCU->getTotalCost()!=MAX_DOUBLE && rpcBestCU->isLosslessCoded(0) && (rpcBestCU->getIPCMFlag(0) == false))
  {
    xFillPCMBuffer(rpcBestCU, m_ppcOrigYuv[uiDepth]);
  }

  if( uiDepth == pps.getMaxCuDQPDepth() )
  {
    Int idQP = m_pcEncCfg->getMaxDeltaQP();
    iMinQP = Clip3( -sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP-idQP );
    iMaxQP = Clip3( -sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP+idQP );
  }
  else if( uiDepth < pps.getMaxCuDQPDepth() )
  {
    iMinQP = iBaseQP;
    iMaxQP = iBaseQP;
  }
  else
  {
    const Int iStartQP = rpcTempCU->getQP(0);
    iMinQP = iStartQP;
    iMaxQP = iStartQP;
  }

  if ( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() )
  {
    iMinQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP - m_lumaQPOffset);
    iMaxQP = iMinQP;
  }

  if ( m_pcEncCfg->getUseRateCtrl() )
  {
    iMinQP = m_pcRateCtrl->getRCQP();
    iMaxQP = m_pcRateCtrl->getRCQP();
  }

  if ( m_pcEncCfg->getCUTransquantBypassFlagForceValue() )
  {
    iMaxQP = iMinQP; // If all TUs are forced into using transquant bypass, do not loop here.
  }

  Bool bSubBranch = bBoundary || !( m_pcEncCfg->getUseEarlyCU() && rpcBestCU->getTotalCost() != MAX_DOUBLE &&
                                    rpcBestCU->isSkipped( 0 ) );

  if ( rpcBestCU->isIntraBC( 0 ) && rpcBestCU->getQtRootCbf( 0 ) == 0 ) { bSubBranch = false; }
   
  #ifdef SDMTEST  // MesksCode
  // Compute residuals of split.
  if ( bSubBranch && OorGorA >= 0 && uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() && LFCNSWITCH) {
      if ( rpcBestCU->getSlice()->getSliceType() == P_SLICE ) {
        Pel* pReco = m_ppcRecoYuvBest[uiDepth]->getAddr( COMPONENT_Y );
        Pel* pPred = m_ppcPredYuvBest[uiDepth]->getAddr( COMPONENT_Y );
        Pel* pOri  = m_ppcOrigYuv[uiDepth]->getAddr( COMPONENT_Y );
        Pel* pOriU = m_ppcOrigYuv[uiDepth]->getAddr( COMPONENT_Cb );

        int   cuHeight = rpcBestCU->getHeight( uiDepth );
        int   cuWidth  = rpcBestCU->getWidth( uiDepth );
        int** pixels;
        int   x_begin = 0;
        int   x_end   = cuWidth;
        int   y_begin = 0;
        int   y_end   = cuHeight;

        pixels = new int*[cuHeight];
        for ( int i = 0; i < cuHeight; i++ ) { pixels[i] = new int[cuWidth]; }

        // processing geometry and attribute separately
        if ( OorGorA == 0 ) {
          for ( int y = 0; y < cuHeight; y++ ) {
            for ( int x = 0; x < cuWidth; x++ ) {
              // computational prediction distortion
              pixels[y][x] = 24 * abs( ( *pOri ) - ( *pPred ) ) / QP;
              pPred++;
              pOri++;
            }
          }
        } else if ( OorGorA > 0 ) {
          for ( int y = 0; y < cuHeight; y++ ) {
            for ( int x = 0; x < cuWidth; x++ ) {
              // computational prediction distortion
              pixels[y][x] = 32 * abs( ( *pOri ) - ( *pPred ) ) / QP;
              pPred++;
              pOri++;
            }
          }
        }

        // normalization processing
        double normalizingFactor = 0;
        if ( OorGorA == 0 )
          normalizingFactor = 20;
        else if ( OorGorA > 0 )
          normalizingFactor = 200;

        int kernelSize = 0;
        if ( uiDepth == 0 )
          kernelSize = 4;
        else if ( uiDepth == 1 )
          kernelSize = 2;
        else if ( uiDepth == 2 )
          kernelSize = 1;

        double** pdm = new double*[16];
        for ( int i = 0; i < 16; i++ ) pdm[i] = new double[16];
        for ( int i = 0; i < cuHeight; i += kernelSize ) {
          for ( int j = 0; j < cuWidth; j += kernelSize ) {
            double maxDistortion = -1;
            for ( int k = i; k < i + kernelSize; k++ ) {
              for ( int l = j; l < j + kernelSize; l++ ) {
                if ( pixels[k][l] > maxDistortion ) maxDistortion = pixels[k][l];
              }
            }
            pdm[i / kernelSize][j / kernelSize] =
                ( maxDistortion / normalizingFactor ) <= 1 ? ( maxDistortion / normalizingFactor ) : 1;
          }
        }

        double    AVERAGE         = 0;
        const int SPLIT           = 2;
        double    overallVariance = 0;
        AVERAGE                   = cacAverage( 0, 0, cuWidth, cuHeight, pixels );
        overallVariance           = cacVariance( 0, 0, cuWidth, cuHeight, pixels, AVERAGE );
        for ( int i = 0; i < SPLIT; i++ ) {
          for ( int j = 0; j < SPLIT; j++ ) {
            AVERAGE            = cacAverage( j * cuWidth / SPLIT, i * cuHeight / SPLIT, ( j + 1 ) * cuWidth / SPLIT,
                                  ( i + 1 ) * cuHeight / SPLIT, pixels );
            double subVariance = cacVariance( j * cuWidth / SPLIT, i * cuHeight / SPLIT, ( j + 1 ) * cuWidth / SPLIT,
                                              ( i + 1 ) * cuHeight / SPLIT, pixels, AVERAGE );
            if ( subVariance > overallVariance ) { overallVariance = subVariance; }
          }
        }

        double statisticVarMax = twoDecimalDouble( overallVariance / normalizingFactor );  // DV
        if ( statisticVarMax > 1 ) statisticVarMax = 1;
        int    statisticCBF    = ( rpcBestCU->getQtRootCbf( 0 ) == 0 ) ? 0 : 1;  // CBF
        double statisticDepth  = uiDepth == 0 ? 2 : ( uiDepth == 1 ? 1 : ( uiDepth == 2 ? 0 : -1 ) );
        statisticDepth         = twoDecimalDouble( statisticDepth / 2.000 );  // CD
        double statisticQP     = twoDecimalDouble( ( 51 - QP ) / 51.000 );    // QP
        double statisticCUcate = CUcate / 2.000;                              // CUC
        // int    statisticPOC    = ( POC % 2 == 0 ) ? 1 : 0;                    // POC
        // double statisticPredictionMode = rpcBestCU->getPredictionMode( 0 );
        // gain front mode
        // map<string, int>::iterator iter;
        // string                     frontKey;
        // double                     statisticFrontMode = -1;
        //// =0 dont split, =1 further split, =-2 means error, =-1 means have not front mode
        // if ( uiDepth == 0 )
        //  statisticFrontMode = -1;
        // else {
        //  int          nowX = int( uiLPelX / pow( 2, 6 - uiDepth + 1 ) ) * pow( 2, 6 - uiDepth + 1 );
        //  int          nowY = int( uiTPelY / pow( 2, 6 - uiDepth + 1 ) ) * pow( 2, 6 - uiDepth + 1 );
        //  stringstream ss_frontKey;
        //  ss_frontKey << POC << "_" << nowX << "_" << nowY << "_" << uiDepth - 1;
        //  ss_frontKey >> frontKey;
        //  iter = frontModeFlag.find( frontKey );
        //  if ( iter != frontModeFlag.end() )
        //    statisticFrontMode = iter->second;
        //  else {
        //    statisticFrontMode = -2;
        //  }
        //}

        // statisticFrontMode =
        //    ( statisticFrontMode == 0 )
        //        ? 0
        //        : ( ( statisticFrontMode == -1 ) ? 0.5 : ( ( statisticFrontMode == 1 ) ? 0.7 : 1 ) );  // FLM

        std::vector<double> x( ( OorGorA == 0 ) ? P_GEO_INPUT : P_ATT_INPUT );
        if ( OorGorA == 0 ) {
          x[0] = statisticVarMax;
          x[1] = statisticDepth;
          x[2] = statisticQP;
          //x[0] = statisticVarMax;
          //x[1] = statisticCBF;
          //x[2] = statisticDepth;
          //x[3] = statisticQP;
          //x[4] = statisticCUcate;
        } else if ( OorGorA > 0 ) {
          x[0] = statisticVarMax;
          x[1] = statisticDepth;
          x[2] = statisticQP;
          //x[0] = statisticVarMax;
          //x[1] = statisticCBF;
          //x[2] = statisticDepth;
          //x[3] = statisticQP;
          //x[4] = statisticCUcate;
        }

        double y = vanillayNN( "PModule", OorGorA, x, "sigmoid" );
        //double y = SLPADDCNN( pdm, x, OorGorA, "PModule" );
        if ( OorGorA == 0 && y < PGeoTH ) 
            bSubBranch = false;
        else if ( OorGorA > 0 && y < PAttTH )
          bSubBranch = false;
        //cout << "P Frame Geometry: " << endl << "  y: " << y << endl << "  autoTH: " << autoTH( QP, 1 ) << endl;

        //cout << "P Frame: --" <<OorGorA<< endl << "pdm: " << endl;
        //for ( int i = 0; i < 16; i++ )
        //  for ( int j = 0; j < 16; j++ ) cout << pdm[i][j] << ",";
        //cout << endl << "\tfx: " << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << ","<<x[4]<< endl << "y: " << y << endl;
      }
      if ( rpcBestCU->getSlice()->getSliceType() == I_SLICE ) {
        Pel* pReco = m_ppcRecoYuvBest[uiDepth]->getAddr( COMPONENT_Y );
        Pel* pPred = m_ppcPredYuvBest[uiDepth]->getAddr( COMPONENT_Y );
        Pel* pOri  = m_ppcOrigYuv[uiDepth]->getAddr( COMPONENT_Y );
        Pel* pOriU = m_ppcOrigYuv[uiDepth]->getAddr( COMPONENT_Cb );

        int   cuHeight = rpcBestCU->getHeight( uiDepth );
        int   cuWidth  = rpcBestCU->getWidth( uiDepth );
        int** pixels;
        int   x_begin = 0;
        int   x_end   = cuWidth;
        int   y_begin = 0;
        int   y_end   = cuHeight;

        pixels = new int*[cuHeight];
        for ( int i = 0; i < cuHeight; i++ ) { pixels[i] = new int[cuWidth]; }

        // processing geometry and attribute separately
        if ( OorGorA == 0 ) {
          for ( int y = 0; y < cuHeight; y++ ) {
            for ( int x = 0; x < cuWidth; x++ ) {
              // computational prediction distortion
              pixels[y][x] = 24 * abs( ( *pOri ) - ( *pPred ) ) / QP;
              pPred++;
              pOri++;
            }
          }
        } else if ( OorGorA > 0 ) {
          for ( int y = 0; y < cuHeight; y++ ) {
            for ( int x = 0; x < cuWidth; x++ ) {
              // computational prediction distortion
              pixels[y][x] = 32 * abs( ( *pOri ) - ( *pPred ) ) / QP;
              pPred++;
              pOri++;
            }
          }
        }

        // normalization processing
        double normalizingFactor = 0;
        if ( OorGorA == 0 )
          normalizingFactor = 20;
        else if ( OorGorA > 0 )
          normalizingFactor = 200;

        int kernelSize = 0;
        if ( uiDepth == 0 )
          kernelSize = 4;
        else if ( uiDepth == 1 )
          kernelSize = 2;
        else if ( uiDepth == 2 )
          kernelSize = 1;

        double** pdm = new double*[16];
        for ( int i = 0; i < 16; i++ ) pdm[i] = new double[16];
        for ( int i = 0; i < cuHeight; i += kernelSize ) {
          for ( int j = 0; j < cuWidth; j += kernelSize ) {
            double maxDistortion = -1;
            for ( int k = i; k < i + kernelSize; k++ ) {
              for ( int l = j; l < j + kernelSize; l++ ) {
                if ( pixels[k][l] > maxDistortion ) maxDistortion = pixels[k][l];
              }
            }
            pdm[i / kernelSize][j / kernelSize] =
                ( maxDistortion / normalizingFactor ) <= 1 ? ( maxDistortion / normalizingFactor ) : 1;
          }
        }

        double    AVERAGE         = 0;
        const int SPLIT           = 2;
        double    overallVariance = 0;
        AVERAGE                   = cacAverage( 0, 0, cuWidth, cuHeight, pixels );
        overallVariance           = cacVariance( 0, 0, cuWidth, cuHeight, pixels, AVERAGE );
        for ( int i = 0; i < SPLIT; i++ ) {
          for ( int j = 0; j < SPLIT; j++ ) {
            AVERAGE            = cacAverage( j * cuWidth / SPLIT, i * cuHeight / SPLIT, ( j + 1 ) * cuWidth / SPLIT,
                                  ( i + 1 ) * cuHeight / SPLIT, pixels );
            double subVariance = cacVariance( j * cuWidth / SPLIT, i * cuHeight / SPLIT, ( j + 1 ) * cuWidth / SPLIT,
                                              ( i + 1 ) * cuHeight / SPLIT, pixels, AVERAGE );
            if ( subVariance > overallVariance ) { overallVariance = subVariance; }
          }
        }

        double statisticVarMax = twoDecimalDouble( overallVariance / normalizingFactor );  // DV
        if ( statisticVarMax > 1 ) statisticVarMax = 1;
        int    statisticCBF    = ( rpcBestCU->getQtRootCbf( 0 ) == 0 ) ? 0 : 1;  // CBF
        double statisticDepth  = uiDepth == 0 ? 2 : ( uiDepth == 1 ? 1 : ( uiDepth == 2 ? 0 : -1 ) );
        statisticDepth         = twoDecimalDouble( statisticDepth / 2.000 );  // CD
        double statisticQP     = twoDecimalDouble( ( 51 - QP ) / 51.000 );    // QP
        double statisticCUcate = CUcate / 2.000;                              // CUC
        // int    statisticPOC    = ( POC % 2 == 0 ) ? 1 : 0;                    // POC
        // double statisticPredictionMode = rpcBestCU->getPredictionMode( 0 );
        //// gain front mode
        // map<string, int>::iterator iter;
        // string                     frontKey;
        // double                     statisticFrontMode = -1;
        //// =0 dont split, =1 further split, =-2 means error, =-1 means have not front mode
        // if ( uiDepth == 0 )
        //  statisticFrontMode = -1;
        // else {
        //  int          nowX = int( uiLPelX / pow( 2, 6 - uiDepth + 1 ) ) * pow( 2, 6 - uiDepth + 1 );
        //  int          nowY = int( uiTPelY / pow( 2, 6 - uiDepth + 1 ) ) * pow( 2, 6 - uiDepth + 1 );
        //  stringstream ss_frontKey;
        //  ss_frontKey << POC << "_" << nowX << "_" << nowY << "_" << uiDepth - 1;
        //  ss_frontKey >> frontKey;
        //  iter = frontModeFlag.find( frontKey );
        //  if ( iter != frontModeFlag.end() )
        //    statisticFrontMode = iter->second;
        //  else {
        //    statisticFrontMode = -2;
        //  }
        //}

        // statisticFrontMode =
        //    ( statisticFrontMode == 0 )
        //        ? 0
        //        : ( ( statisticFrontMode == -1 ) ? 0.5 : ( ( statisticFrontMode == 1 ) ? 0.7 : 1 ) );  // FLM
        std::vector<double> x( ( OorGorA == 0 ) ? I_GEO_INPUT : I_ATT_INPUT );
        if ( OorGorA == 0 ) {
          x[0] = statisticVarMax;
          x[1] = statisticDepth;
          x[2] = statisticQP;
          //x[0] = statisticVarMax;
          //x[1] = statisticCBF;
          //x[2] = statisticDepth;
          //x[3] = statisticQP;
          //x[4] = statisticCUcate;
        } else if ( OorGorA > 0 ) {
          x[0] = statisticVarMax;
          x[1] = statisticDepth;
          x[2] = statisticQP;
          //x[0] = statisticVarMax;
          //x[1] = statisticCBF;
          //x[2] = statisticDepth;
          //x[3] = statisticQP;
          //x[4] = statisticCUcate;
        }

         double y = vanillayNN( "IModule", OorGorA, x, "sigmoid" );
        //double y = SLPADDCNN( pdm, x, OorGorA, "IModule" );
        if ( OorGorA == 0 && y < IGeoTH )
          bSubBranch = false;
        else if ( OorGorA > 0 && y < IAttTH )
          bSubBranch = false;

        //cout << "I Frame: (OorGorA: " << OorGorA << " )" << endl
        //     << "  y: " << y << endl
        //     << "  autoTH: " << autoTH( QP, 0 ) << endl;
        //cout << "I Frame: --" << OorGorA << endl << "pdm: " << endl;
        //for(int i=0;i<16;i++)
        //  for (int j=0;j<16;j++)
        //    cout << pdm[i][j] << ",";
        //cout << endl << "\tfx: " << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << ","<<x[4]<< endl << "y: " << y << endl;

        //cout << "MesksCode: " << endl << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] << endl << "y: "<<y<<endl;
      }
    }
#endif  // SDMTEST

  if( !terminateAllFurtherRDO && bSubBranch && uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() && (!getFastDeltaQp() || uiWidth > fastDeltaQPCuMaxSize || bBoundary))
  {
    PaletteInfoBuffer tempPalettePredictor;

    if( iMinQP != iMaxQP )
    {
      memcpy( tempPalettePredictor.lastPaletteSize, lastPaletteSize, sizeof( lastPaletteSize ) );
      memcpy( tempPalettePredictor.lastPalette[0],  lastPalette[0], sizeof( tempPalettePredictor.lastPalette[0] ) );
      memcpy( tempPalettePredictor.lastPalette[1],  lastPalette[1], sizeof( tempPalettePredictor.lastPalette[1] ) );
      memcpy( tempPalettePredictor.lastPalette[2],  lastPalette[2], sizeof( tempPalettePredictor.lastPalette[2] ) );
    }

    // further split
    Double splitTotalCost = 0;

    for (Int iQP=iMinQP; iQP<=iMaxQP; iQP++)
    {
      const Bool bIsLosslessMode = false; // False at this level. Next level down may set it to true.

      rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

      UChar       uhNextDepth         = uiDepth+1;
      TComDataCU* pcSubBestPartCU     = m_ppcBestCU[uhNextDepth];
      TComDataCU* pcSubTempPartCU     = m_ppcTempCU[uhNextDepth];
      DEBUG_STRING_NEW(sTempDebug)

      if( iMinQP != iMaxQP && iQP != iMinQP )
      {
        memcpy( lastPaletteSize, tempPalettePredictor.lastPaletteSize, sizeof( lastPaletteSize ) );
        memcpy( lastPalette[0],  tempPalettePredictor.lastPalette[0],  sizeof( tempPalettePredictor.lastPalette[0] ) );
        memcpy( lastPalette[1],  tempPalettePredictor.lastPalette[1],  sizeof( tempPalettePredictor.lastPalette[1] ) );
        memcpy( lastPalette[2],  tempPalettePredictor.lastPalette[2],  sizeof( tempPalettePredictor.lastPalette[2] ) );
      }

      for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++ )
      {
        pcSubBestPartCU->initSubCU( rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP );           // clear sub partition datas or init.
        pcSubTempPartCU->initSubCU( rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP );           // clear sub partition datas or init.

        for(Int i=0; i<2; i++)
        {
          pcSubBestPartCU->setLastIntraBCMv( lastIntraBCMv[i], i );
          pcSubTempPartCU->setLastIntraBCMv( lastIntraBCMv[i], i );
        }

        for (UInt ch = 0; ch < numValidComp; ch++)
        {
          pcSubBestPartCU->setLastPaletteInLcuSizeFinal(ch, lastPaletteSize[ch]);
          pcSubTempPartCU->setLastPaletteInLcuSizeFinal(ch, lastPaletteSize[ch]);
          for (UInt i = 0; i < pcSlice->getSPS()->getSpsScreenExtension().getPaletteMaxPredSize(); i++)
          {
            pcSubBestPartCU->setLastPaletteInLcuFinal(ch, lastPalette[ch][i], i);
            pcSubTempPartCU->setLastPaletteInLcuFinal(ch, lastPalette[ch][i], i);
          }
        }

        if( ( pcSubBestPartCU->getCUPelX() < sps.getPicWidthInLumaSamples() )&&
             ( pcSubBestPartCU->getCUPelY() < sps.getPicHeightInLumaSamples() ) )
        {
          if ( 0 == uiPartUnitIdx) //initialize RD with previous depth buffer
          {
            m_pppcRDSbacCoder[uhNextDepth][CI_CURR_BEST]->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
          }
          else
          {
            m_pppcRDSbacCoder[uhNextDepth][CI_CURR_BEST]->load(m_pppcRDSbacCoder[uhNextDepth][CI_NEXT_BEST]);
          }

#if AMP_ENC_SPEEDUP
          DEBUG_STRING_NEW(sChild)
          if ( !(rpcBestCU->getTotalCost()!=MAX_DOUBLE && rpcBestCU->isInter(0)) )
          {
            xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uhNextDepth DEBUG_STRING_PASS_INTO(sChild), NUMBER_OF_PART_SIZES );
          }
          else
          {

            xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uhNextDepth DEBUG_STRING_PASS_INTO(sChild), rpcBestCU->getPartitionSize(0) );
          }
          DEBUG_STRING_APPEND(sTempDebug, sChild)
#else
          xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uhNextDepth );
#endif

          // NOTE: (0,0) is used as an indicator that IntraBC has not been used within the CU.
          if( pcSubBestPartCU->getLastIntraBCMv().getHor() != 0 || pcSubBestPartCU->getLastIntraBCMv().getVer() != 0 )
          {
            for(Int i=0; i<2; i++)
            {
              lastIntraBCMv[i] = pcSubBestPartCU->getLastIntraBCMv(i);
            }
          }

          if (pcSubBestPartCU->getLastPaletteInLcuSizeFinal(COMPONENT_Y) != 0)
          {
            for (UInt ch = 0; ch < numValidComp; ch++)
            {
              lastPaletteSize[ch] = pcSubBestPartCU->getLastPaletteInLcuSizeFinal(ch);
              for (UInt i = 0; i < pcSlice->getSPS()->getSpsScreenExtension().getPaletteMaxPredSize(); i++)
              {

                lastPalette[ch][i] = pcSubBestPartCU->getLastPaletteInLcuFinal(ch, i);
              }
            }
          }

          rpcTempCU->copyPartFrom( pcSubBestPartCU, uiPartUnitIdx, uhNextDepth );         // Keep best part data to current temporary data.
          xCopyYuv2Tmp( pcSubBestPartCU->getTotalNumPart()*uiPartUnitIdx, uhNextDepth );
          if ( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() && pps.getMaxCuDQPDepth() >= 1 )
          {
            splitTotalCost += pcSubBestPartCU->getTotalCost();
          }
        }
        else
        {
          pcSubBestPartCU->copyToPic( uhNextDepth );
          rpcTempCU->copyPartFrom( pcSubBestPartCU, uiPartUnitIdx, uhNextDepth );
        }
      }

      m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uhNextDepth][CI_NEXT_BEST]);
      if( !bBoundary )
      {
        m_pcEntropyCoder->resetBits();
        m_pcEntropyCoder->encodeSplitFlag( rpcTempCU, 0, uiDepth, true );
        if ( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() && pps.getMaxCuDQPDepth() >= 1 )
        {
          Int splitBits = m_pcEntropyCoder->getNumberOfWrittenBits();
          Double splitBitCost = m_pcRdCost->calcRdCost( splitBits, 0 );
          splitTotalCost += splitBitCost;
        }

        rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
        rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
      }

      if ( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() && pps.getMaxCuDQPDepth() >= 1 )
      {
        rpcTempCU->getTotalCost() = splitTotalCost;
      }
      else
      {
        rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );
      }

      if( uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
      {
        Bool hasResidual = false;
        for( UInt uiBlkIdx = 0; uiBlkIdx < rpcTempCU->getTotalNumPart(); uiBlkIdx ++)
        {
          if( (     rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Y)
                || (rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Cb) && (numberValidComponents > COMPONENT_Cb))
                || (rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Cr) && (numberValidComponents > COMPONENT_Cr)) ) 
                || ( rpcTempCU->getPaletteModeFlag(uiBlkIdx) && rpcTempCU->getPaletteEscape(COMPONENT_Y, uiBlkIdx) ) )
          {
            hasResidual = true;
            break;
          }
        }

        if ( hasResidual )
        {
          m_pcEntropyCoder->resetBits();
          m_pcEntropyCoder->encodeQP( rpcTempCU, 0, false );
          rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
          rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
          rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

          Bool foundNonZeroCbf = false;
          rpcTempCU->setQPSubCUs( rpcTempCU->getRefQP( 0 ), 0, uiDepth, foundNonZeroCbf );
          assert( foundNonZeroCbf );
        }
        else
        {
          rpcTempCU->setQPSubParts( rpcTempCU->getRefQP( 0 ), 0, uiDepth ); // set QP to default QP
        }
      }

      m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);

      // If the configuration being tested exceeds the maximum number of bytes for a slice / slice-segment, then
      // a proper RD evaluation cannot be performed. Therefore, termination of the
      // slice/slice-segment must be made prior to this CTU.
      // This can be achieved by forcing the decision to be that of the rpcTempCU.
      // The exception is each slice / slice-segment must have at least one CTU.
      if (rpcBestCU->getTotalCost()!=MAX_DOUBLE)
      {
        const Bool isEndOfSlice        =    pcSlice->getSliceMode()==FIXED_NUMBER_OF_BYTES
                                         && ((pcSlice->getSliceBits()+rpcBestCU->getTotalBits())>pcSlice->getSliceArgument()<<3)
                                         && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceCurStartCtuTsAddr())
                                         && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
        const Bool isEndOfSliceSegment =    pcSlice->getSliceSegmentMode()==FIXED_NUMBER_OF_BYTES
                                         && ((pcSlice->getSliceSegmentBits()+rpcBestCU->getTotalBits()) > pcSlice->getSliceSegmentArgument()<<3)
                                         && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
                                             // Do not need to check slice condition for slice-segment since a slice-segment is a subset of a slice.
        if(isEndOfSlice||isEndOfSliceSegment)
        {
          rpcBestCU->getTotalCost()=MAX_DOUBLE;
        }
      }

      xCheckBestMode( rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTempDebug) DEBUG_STRING_PASS_INTO(false) ); // RD compare current larger prediction
                                                                                                                                                       // with sub partitioned prediction.
    }
  }

  DEBUG_STRING_APPEND(sDebug_, sDebug);

  rpcBestCU->copyToPic(uiDepth);                                                     // Copy Best data to Picture for next partition prediction.

  xCopyYuv2Pic( rpcBestCU->getPic(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu(), uiDepth, uiDepth, rpcBestCU );   // Copy Yuv data to picture Yuv
  for (UInt ch = 0; ch < numValidComp; ch++)
  {
    if (lastPalette[ch])
    {
      xFree(lastPalette[ch]);
      lastPalette[ch] = NULL;
    }
  }
  if (bBoundary)
  {
    return;
  }
  // Assert if Best prediction mode is NONE
  // Selected mode's RD-cost must be not MAX_DOUBLE.
  assert( rpcBestCU->getPartitionSize ( 0 ) != NUMBER_OF_PART_SIZES       );
  assert( rpcBestCU->getPredictionMode( 0 ) != NUMBER_OF_PREDICTION_MODES );
  assert( rpcBestCU->getTotalCost     (   ) != MAX_DOUBLE                 );
}

/** finish encoding a cu and handle end-of-slice conditions
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TEncCu::finishCU( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  TComPic* pcPic = pcCU->getPic();
  TComSlice * pcSlice = pcCU->getPic()->getSlice(pcCU->getPic()->getCurrSliceIdx());

  //Calculate end address
  const Int  currentCTUTsAddr = pcPic->getPicSym()->getCtuRsToTsAddrMap(pcCU->getCtuRsAddr());
  const Bool isLastSubCUOfCtu = pcCU->isLastSubCUOfCtu(uiAbsPartIdx);
  if ( isLastSubCUOfCtu )
  {
    // The 1-terminating bit is added to all streams, so don't add it here when it's 1.
    // i.e. when the slice segment CurEnd CTU address is the current CTU address+1.
    if (pcSlice->getSliceSegmentCurEndCtuTsAddr() != currentCTUTsAddr+1)
    {
      m_pcEntropyCoder->encodeTerminatingBit( 0 );
    }
  }
}

/** Compute QP for each CU
 * \param pcCU Target CU
 * \param uiDepth CU depth
 * \returns quantization parameter
 */
Int TEncCu::xComputeQP( TComDataCU* pcCU, UInt uiDepth )
{
  Int iBaseQp = pcCU->getSlice()->getSliceQp();
  Int iQpOffset = 0;
  if ( m_pcEncCfg->getUseAdaptiveQP() )
  {
    TEncPic* pcEPic = dynamic_cast<TEncPic*>( pcCU->getPic() );
    UInt uiAQDepth = min( uiDepth, pcEPic->getMaxAQDepth()-1 );
    TEncPicQPAdaptationLayer* pcAQLayer = pcEPic->getAQLayer( uiAQDepth );
    UInt uiAQUPosX = pcCU->getCUPelX() / pcAQLayer->getAQPartWidth();
    UInt uiAQUPosY = pcCU->getCUPelY() / pcAQLayer->getAQPartHeight();
    UInt uiAQUStride = pcAQLayer->getAQPartStride();
    TEncQPAdaptationUnit* acAQU = pcAQLayer->getQPAdaptationUnit();

    Double dMaxQScale = pow(2.0, m_pcEncCfg->getQPAdaptationRange()/6.0);
    Double dAvgAct = pcAQLayer->getAvgActivity();
    Double dCUAct = acAQU[uiAQUPosY * uiAQUStride + uiAQUPosX].getActivity();
    Double dNormAct = (dMaxQScale*dCUAct + dAvgAct) / (dCUAct + dMaxQScale*dAvgAct);
    Double dQpOffset = log(dNormAct) / log(2.0) * 6.0;
    iQpOffset = Int(floor( dQpOffset + 0.49999 ));
  }

  return Clip3(-pcCU->getSlice()->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQp+iQpOffset );
}

/** encode a CU block recursively
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TEncCu::xEncodeCU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
        TComPic   *const pcPic   = pcCU->getPic();
        TComSlice *const pcSlice = pcCU->getSlice();
  const TComSPS   &sps =*(pcSlice->getSPS());
  const TComPPS   &pps =*(pcSlice->getPPS());

  const UInt maxCUWidth  = sps.getMaxCUWidth();
  const UInt maxCUHeight = sps.getMaxCUHeight();

        Bool bBoundary = false;
        UInt uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
  const UInt uiRPelX   = uiLPelX + (maxCUWidth>>uiDepth)  - 1;
        UInt uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
  const UInt uiBPelY   = uiTPelY + (maxCUHeight>>uiDepth) - 1;

  if( ( uiRPelX < sps.getPicWidthInLumaSamples() ) && ( uiBPelY < sps.getPicHeightInLumaSamples() ) )
  {
    m_pcEntropyCoder->encodeSplitFlag( pcCU, uiAbsPartIdx, uiDepth );
  }
  else
  {
    bBoundary = true;
  }

  if( ( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) && ( uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() ) ) || bBoundary )
  {
    UInt uiQNumParts = ( pcPic->getNumPartitionsInCtu() >> (uiDepth<<1) )>>2;
    if( uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
    {
      setdQPFlag(true);
    }

    if( uiDepth == pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcSlice->getUseChromaQpAdj())
    {
      setCodeChromaQpAdjFlag(true);
    }

    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx+=uiQNumParts )
    {
      uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
      uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
      if( ( uiLPelX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY < sps.getPicHeightInLumaSamples() ) )
      {
        xEncodeCU( pcCU, uiAbsPartIdx, uiDepth+1 );
      }
    }
    return;
  }

  if( uiDepth <= pps.getMaxCuDQPDepth() && pps.getUseDQP())
  {
    setdQPFlag(true);
  }

  if( uiDepth <= pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcSlice->getUseChromaQpAdj())
  {
    setCodeChromaQpAdjFlag(true);
  }

  if (pps.getTransquantBypassEnabledFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( pcCU, uiAbsPartIdx );
  }

  if( !pcSlice->isIntra() )
  {
    m_pcEntropyCoder->encodeSkipFlag( pcCU, uiAbsPartIdx );
  }

  if( pcCU->isSkipped( uiAbsPartIdx ) )
  {
    m_pcEntropyCoder->encodeMergeIndex( pcCU, uiAbsPartIdx );
    finishCU(pcCU,uiAbsPartIdx);
    return;
  }

  m_pcEntropyCoder->encodePredMode( pcCU, uiAbsPartIdx );

  Bool bCodeDQP = getdQPFlag();
  Bool codeChromaQpAdj = getCodeChromaQpAdjFlag();

  m_pcEntropyCoder->encodePaletteModeInfo( pcCU, uiAbsPartIdx, false, &bCodeDQP, &codeChromaQpAdj );

  if ( pcCU->getPaletteModeFlag(uiAbsPartIdx) )
  {
    setCodeChromaQpAdjFlag( codeChromaQpAdj );
    setdQPFlag( bCodeDQP );
    finishCU(pcCU,uiAbsPartIdx);
    return;
  }

  m_pcEntropyCoder->encodePartSize( pcCU, uiAbsPartIdx, uiDepth );

  if (pcCU->isIntra( uiAbsPartIdx ) && pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
  {
    m_pcEntropyCoder->encodeIPCMInfo( pcCU, uiAbsPartIdx );

    if(pcCU->getIPCMFlag(uiAbsPartIdx))
    {
      // Encode slice finish
      finishCU(pcCU,uiAbsPartIdx);
      return;
    }
  }

  // prediction Info ( Intra : direction mode, Inter : Mv, reference idx )
  m_pcEntropyCoder->encodePredInfo( pcCU, uiAbsPartIdx );

  // Encode Coefficients
  m_pcEntropyCoder->encodeCoeff( pcCU, uiAbsPartIdx, uiDepth, bCodeDQP, codeChromaQpAdj );
  setCodeChromaQpAdjFlag( codeChromaQpAdj );
  setdQPFlag( bCodeDQP );

  // --- write terminating bit ---
  finishCU(pcCU,uiAbsPartIdx);
}

Int xCalcHADs8x8_ISlice(Pel *piOrg, Int iStrideOrg)
{
  Int k, i, j, jj;
  Int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;

  for( k = 0; k < 64; k += 8 )
  {
    diff[k+0] = piOrg[0] ;
    diff[k+1] = piOrg[1] ;
    diff[k+2] = piOrg[2] ;
    diff[k+3] = piOrg[3] ;
    diff[k+4] = piOrg[4] ;
    diff[k+5] = piOrg[5] ;
    diff[k+6] = piOrg[6] ;
    diff[k+7] = piOrg[7] ;

    piOrg += iStrideOrg;
  }

  //horizontal
  for (j=0; j < 8; j++)
  {
    jj = j << 3;
    m2[j][0] = diff[jj  ] + diff[jj+4];
    m2[j][1] = diff[jj+1] + diff[jj+5];
    m2[j][2] = diff[jj+2] + diff[jj+6];
    m2[j][3] = diff[jj+3] + diff[jj+7];
    m2[j][4] = diff[jj  ] - diff[jj+4];
    m2[j][5] = diff[jj+1] - diff[jj+5];
    m2[j][6] = diff[jj+2] - diff[jj+6];
    m2[j][7] = diff[jj+3] - diff[jj+7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  //vertical
  for (i=0; i < 8; i++)
  {
    m3[0][i] = m2[0][i] + m2[4][i];
    m3[1][i] = m2[1][i] + m2[5][i];
    m3[2][i] = m2[2][i] + m2[6][i];
    m3[3][i] = m2[3][i] + m2[7][i];
    m3[4][i] = m2[0][i] - m2[4][i];
    m3[5][i] = m2[1][i] - m2[5][i];
    m3[6][i] = m2[2][i] - m2[6][i];
    m3[7][i] = m2[3][i] - m2[7][i];

    m1[0][i] = m3[0][i] + m3[2][i];
    m1[1][i] = m3[1][i] + m3[3][i];
    m1[2][i] = m3[0][i] - m3[2][i];
    m1[3][i] = m3[1][i] - m3[3][i];
    m1[4][i] = m3[4][i] + m3[6][i];
    m1[5][i] = m3[5][i] + m3[7][i];
    m1[6][i] = m3[4][i] - m3[6][i];
    m1[7][i] = m3[5][i] - m3[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      iSumHad += abs(m2[i][j]);
    }
  }
  iSumHad -= abs(m2[0][0]);
  iSumHad =(iSumHad+2)>>2;
  return(iSumHad);
}

Int  TEncCu::updateCtuDataISlice(TComDataCU* pCtu, Int width, Int height)
{
  Int  xBl, yBl;
  const Int iBlkSize = 8;

  Pel* pOrgInit   = pCtu->getPic()->getPicYuvOrg()->getAddr(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);
  Int  iStrideOrig = pCtu->getPic()->getPicYuvOrg()->getStride(COMPONENT_Y);
  Pel  *pOrg;

  Int iSumHad = 0;
  for ( yBl=0; (yBl+iBlkSize)<=height; yBl+= iBlkSize)
  {
    for ( xBl=0; (xBl+iBlkSize)<=width; xBl+= iBlkSize)
    {
      pOrg = pOrgInit + iStrideOrig*yBl + xBl;
      iSumHad += xCalcHADs8x8_ISlice(pOrg, iStrideOrig);
    }
  }
  return(iSumHad);
}

/** check RD costs for a CU block encoded with merge
 * \param rpcBestCU
 * \param rpcTempCU
 * \param earlyDetectionSkipMode
 */
Void TEncCu::xCheckRDCostMerge2Nx2N( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU DEBUG_STRING_FN_DECLARE(sDebug), Bool *earlyDetectionSkipMode, Bool checkSkipOnly )
{
  assert( rpcTempCU->getSlice()->getSliceType() != I_SLICE );
  if(getFastDeltaQp())
  {
    return;   // never check merge in fast deltaqp mode
  }
  TComMvField  cMvFieldNeighbours[2 * MRG_MAX_NUM_CANDS]; // double length for mv of both lists
  UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
  Int numValidMergeCand = 0;
  const Bool bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);

  for( UInt ui = 0; ui < rpcTempCU->getSlice()->getMaxNumMergeCand(); ++ui )
  {
    uhInterDirNeighbours[ui] = 0;
  }
  UChar uhDepth = rpcTempCU->getDepth( 0 );
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth ); // interprets depth relative to CTU level
#if MCTS_ENC_CHECK
  UInt numSpatialMergeCandidates = 0;
  rpcTempCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand, numSpatialMergeCandidates );
#else
  rpcTempCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours,uhInterDirNeighbours, numValidMergeCand );
#endif
  rpcTempCU->roundMergeCandidates(cMvFieldNeighbours, numValidMergeCand);
  rpcTempCU->xRestrictBipredMergeCand(0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand );

  Int mergeCandBuffer[MRG_MAX_NUM_CANDS];
  for( UInt ui = 0; ui < numValidMergeCand; ++ui )
  {
    mergeCandBuffer[ui] = 0;
  }
#if MCTS_ENC_CHECK
  if (m_pcEncCfg->getTMCTSSEITileConstraint() && rpcTempCU->isLastColumnCTUInTile())
  {
    numValidMergeCand = numSpatialMergeCandidates;
  }
#endif

  Bool bestIsSkip = false;

  UInt iteration;
  UInt iterationBegin = checkSkipOnly ? 1 : 0;
  if ( rpcTempCU->isLosslessCoded(0))
  {
    iteration = 1;
    iterationBegin = 0;
  }
  else
  {
    iteration = 2;
  }
  DEBUG_STRING_NEW(bestStr)

  const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
  for( UInt uiNoResidual = iterationBegin; uiNoResidual < iteration; ++uiNoResidual )
  {
    for( UInt uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand )
    {
      if(!(uiNoResidual==1 && mergeCandBuffer[uiMergeCand]==1))
      {
        if( !(bestIsSkip && uiNoResidual == 0) )
        {
          if ( (uhInterDirNeighbours[uiMergeCand] == 1 || uhInterDirNeighbours[uiMergeCand] == 3) && rpcTempCU->getSlice()->getRefPic( REF_PIC_LIST_0, cMvFieldNeighbours[uiMergeCand<<1].getRefIdx() )->getPOC() == rpcTempCU->getSlice()->getPOC() )
          {
            continue;
          }

          DEBUG_STRING_NEW(tmpStr)
          // set MC parameters
          rpcTempCU->setPredModeSubParts( MODE_INTER, 0, uhDepth ); // interprets depth relative to CTU level
          rpcTempCU->setCUTransquantBypassSubParts( bTransquantBypassFlag, 0, uhDepth );
          rpcTempCU->setChromaQpAdjSubParts( bTransquantBypassFlag ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth );
          rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth ); // interprets depth relative to CTU level
          rpcTempCU->setMergeFlagSubParts( true, 0, 0, uhDepth ); // interprets depth relative to CTU level
          rpcTempCU->setMergeIndexSubParts( uiMergeCand, 0, 0, uhDepth ); // interprets depth relative to CTU level
          rpcTempCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeCand], 0, 0, uhDepth ); // interprets depth relative to CTU level
          rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->setAllMvField( cMvFieldNeighbours[0 + 2*uiMergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level
          rpcTempCU->getCUMvField( REF_PIC_LIST_1 )->setAllMvField( cMvFieldNeighbours[1 + 2*uiMergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level

#if MCTS_ENC_CHECK
          if ( m_pcEncCfg->getTMCTSSEITileConstraint () && (!(m_pcPredSearch->checkTMctsMvp(rpcTempCU))))
          {
            continue;
          }

#endif
          // do MC
          m_pcPredSearch->motionCompensation ( rpcTempCU, m_ppcPredYuvTemp[uhDepth] );
          Bool bColourTrans = (m_pcEncCfg->getRGBFormatFlag() && rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans())? true : false;

          if ( bTransquantBypassFlag && (sps.getBitDepth( CHANNEL_TYPE_LUMA ) != sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
          {
            bColourTrans = false;
          }
          TComYuv* pcTmpPredYuv = m_ppcPredYuvTemp[uhDepth];
          if( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() && (!bTransquantBypassFlag || sps.getBitDepth( CHANNEL_TYPE_LUMA ) == sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
          {
            if(!getEnableInterTUACT())
            {
              m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                         m_ppcOrigYuv    [uhDepth],
                                                         pcTmpPredYuv             ,
                                                         m_ppcResiYuvTemp[uhDepth],
                                                         m_ppcResiYuvBest[uhDepth],
                                                         m_ppcRecoYuvTemp[uhDepth],
                                                         (uiNoResidual != 0),
                                                         m_ppcNoCorrYuv  [uhDepth],
                                                         (bColourTrans? ACT_TRAN_CLR: ACT_ORG_CLR)
                                                         DEBUG_STRING_PASS_INTO(tmpStr)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[uhDepth]
#endif
                                                        );
            }
            else
            {
              m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                         m_ppcOrigYuv    [uhDepth],
                                                         pcTmpPredYuv             ,
                                                         m_ppcResiYuvTemp[uhDepth],
                                                         m_ppcResiYuvBest[uhDepth],
                                                         m_ppcRecoYuvTemp[uhDepth],
                                                         (uiNoResidual != 0),
                                                         m_ppcNoCorrYuv  [uhDepth],
                                                         ACT_TWO_CLR
                                                         DEBUG_STRING_PASS_INTO(tmpStr)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[uhDepth]
#endif
                                                        );
            }
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                       m_ppcOrigYuv    [uhDepth],
                                                       pcTmpPredYuv             ,
                                                       m_ppcResiYuvTemp[uhDepth],
                                                       m_ppcResiYuvBest[uhDepth],
                                                       m_ppcRecoYuvTemp[uhDepth],
                                                       (uiNoResidual != 0),
                                                       m_ppcNoCorrYuv  [uhDepth],
                                                       ACT_ORG_CLR
#if PCC_RDO_EXT
                                                       DEBUG_STRING_PASS_INTO(tmpStr),
                                                       m_ppcOccupancyYuv[uhDepth]);
#else
                                                       DEBUG_STRING_PASS_INTO(tmpStr) );
#endif
          }
          rpcTempCU->setSkipFlagSubParts( rpcTempCU->getQtRootCbf(0) == 0, 0, uhDepth );

#if DEBUG_STRING
          DebugInterPredResiReco(tmpStr, *(m_ppcPredYuvTemp[uhDepth]), *(m_ppcResiYuvBest[uhDepth]), *(m_ppcRecoYuvTemp[uhDepth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif

          if ((uiNoResidual == 0) && (rpcTempCU->getQtRootCbf(0) == 0))
          {
            // If no residual when allowing for one, then set mark to not try case where residual is forced to 0
            mergeCandBuffer[uiMergeCand] = 1;
          }

          Double rdCost1 = rpcTempCU->getTotalCost();
          if ( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() )
          {
            if ( rdCost1 < m_ppcBestCU[uhDepth]->getTmpInterRDCost() )
            {
              m_ppcBestCU[uhDepth]->setTmpInterRDCost(rdCost1);
              m_ppcBestCU[uhDepth]->setInterCSCEnabled(true);
              m_ppcTempCU[uhDepth]->setTmpInterRDCost(rdCost1);
              m_ppcTempCU[uhDepth]->setInterCSCEnabled(true);
            }
          }

          Int orgQP = rpcTempCU->getQP( 0 );
          xCheckDQP( rpcTempCU );
          TComDataCU *rpcTempCUPre = rpcTempCU;

          xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth DEBUG_STRING_PASS_INTO(bestStr) DEBUG_STRING_PASS_INTO(tmpStr));

          Bool bParentUseCSC = ( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans()
                                 && (uhDepth != 0) && (m_ppcBestCU[uhDepth -1]->getInterCSCEnabled()) );
          bParentUseCSC = bParentUseCSC || rpcBestCU->getQtRootCbf(0) == 0 || rpcBestCU->getMergeFlag( 0 );
          if(rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() && !uiNoResidual && rpcTempCUPre->getQtRootCbf(0) && !bParentUseCSC && (!bTransquantBypassFlag || sps.getBitDepth( CHANNEL_TYPE_LUMA ) == sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
          {
            if( rpcTempCU != rpcTempCUPre )
            {
              rpcTempCU->initEstData( uhDepth, orgQP, bTransquantBypassFlag );
              rpcTempCU->copyPartFrom( rpcBestCU, 0, uhDepth );
            }
            if ( !getEnableInterTUACT() )
            {
              bColourTrans = !bColourTrans;
            }
            else
            {
              bColourTrans = m_pcEncCfg->getRGBFormatFlag();
            }
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                       m_ppcOrigYuv    [uhDepth],
                                                       pcTmpPredYuv             ,
                                                       m_ppcResiYuvTemp[uhDepth],
                                                       m_ppcResiYuvBest[uhDepth],
                                                       m_ppcRecoYuvTemp[uhDepth],
                                                       (uiNoResidual != 0),
                                                       m_ppcNoCorrYuv  [uhDepth],
                                                       (bColourTrans? ACT_TRAN_CLR: ACT_ORG_CLR)
                                                       DEBUG_STRING_PASS_INTO(tmpStr)
#if PCC_RDO_EXT
                                                     , m_ppcOccupancyYuv[uhDepth]
#endif
                                                      );
            rpcTempCU->setSkipFlagSubParts( rpcTempCU->getQtRootCbf(0) == 0, 0, uhDepth );
            Double rdCost = rpcTempCU->getTotalCost();
            if(rdCost < m_ppcBestCU[uhDepth]->getTmpInterRDCost() )
            {
              m_ppcBestCU[uhDepth]->setTmpInterRDCost(rdCost);
              m_ppcBestCU[uhDepth]->setInterCSCEnabled(false);
              m_ppcTempCU[uhDepth]->setTmpInterRDCost(rdCost);
              m_ppcTempCU[uhDepth]->setInterCSCEnabled(false);
            }

            xCheckDQP( rpcTempCU );
            if(rpcTempCU->getQtRootCbf(0))
            {
              xCheckBestMode( rpcBestCU, rpcTempCU, uhDepth DEBUG_STRING_PASS_INTO(bestStr) DEBUG_STRING_PASS_INTO(tmpStr) );
            }
          }
          rpcTempCU->initEstData( uhDepth, orgQP, bTransquantBypassFlag );

          if( m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip )
          {
            bestIsSkip = rpcBestCU->getQtRootCbf(0) == 0;
          }
        }
      }
    }

    if(uiNoResidual == 0 && m_pcEncCfg->getUseEarlySkipDetection())
    {
      if(rpcBestCU->getQtRootCbf( 0 ) == 0)
      {
        if( rpcBestCU->getMergeFlag( 0 ))
        {
          *earlyDetectionSkipMode = true;
        }
        else if(m_pcEncCfg->getMotionEstimationSearchMethod() != MESEARCH_SELECTIVE)
        {
          Int absoulte_MV=0;
          for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            if ( rpcBestCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
              TComCUMvField* pcCUMvField = rpcBestCU->getCUMvField(RefPicList( uiRefListIdx ));
              Int iHor = pcCUMvField->getMvd( 0 ).getAbsHor();
              Int iVer = pcCUMvField->getMvd( 0 ).getAbsVer();
              absoulte_MV+=iHor+iVer;
            }
          }

          if(absoulte_MV == 0)
          {
            *earlyDetectionSkipMode = true;
          }
        }
      }
    }
  }
  DEBUG_STRING_APPEND(sDebug, bestStr)
}


#if AMP_MRG
Void TEncCu::xCheckRDCostInter( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize DEBUG_STRING_FN_DECLARE(sDebug), Bool bUseMRG, TComMv *iMVCandList)
#else
Void TEncCu::xCheckRDCostInter( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize, TComMv *iMVCandList )
#endif
{
  DEBUG_STRING_NEW(sTest)

  if(getFastDeltaQp())
  {
    const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
    const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight()>>(sps.getLog2DiffMaxMinCodingBlockSize()), sps.getMaxCUHeight(), 32u);
    if(ePartSize != SIZE_2Nx2N || rpcTempCU->getWidth( 0 ) > fastDeltaQPCuMaxSize)
    {
      return; // only check necessary 2Nx2N Inter in fast deltaqp mode
    }
  }

  // prior to this, rpcTempCU will have just been reset using rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
  UChar uhDepth = rpcTempCU->getDepth( 0 );

  rpcTempCU->setPartSizeSubParts  ( ePartSize,  0, uhDepth );
  rpcTempCU->setPredModeSubParts  ( MODE_INTER, 0, uhDepth );
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth );

#if MCTS_ENC_CHECK
  rpcTempCU->setTMctsMvpIsValid(true);
#endif

#if AMP_MRG
  rpcTempCU->setMergeAMP (true);
  Bool valid = m_pcPredSearch->predInterSearch ( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth] DEBUG_STRING_PASS_INTO(sTest), false, bUseMRG, iMVCandList );
#else
  Bool valid = m_pcPredSearch->predInterSearch ( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, iMVCandList );
#endif

  if( !valid )
  {
    return;
  }

#if AMP_MRG
  if ( !rpcTempCU->getMergeAMP() )
  {
    return;
  }
#endif

#if MCTS_ENC_CHECK
  if (m_pcEncCfg->getTMCTSSEITileConstraint() && (!rpcTempCU->getTMctsMvpIsValid()))
  {
    return;
  }
#endif

  TComDataCU *rpcTempCUPre = NULL;
  Bool   bEnableTrans      = rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans();

  UChar  colourTransform = 0;
  Bool   bRGB              = m_pcEncCfg->getRGBFormatFlag();
  SChar  orgQP             = rpcTempCU->getQP( 0 );
  Bool   bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);
  const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
  if ( bTransquantBypassFlag && (sps.getBitDepth( CHANNEL_TYPE_LUMA ) != sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
  {
    bEnableTrans = false;
  }
  TComYuv* pcTmpPredYuv = m_ppcPredYuvTemp[uhDepth];
  for(UInt i = 0;  i < 2 ; i++)
  {
    colourTransform = (bRGB && bEnableTrans)? (1-i): i;
    if ( i == 0 )
    {
      if ( bEnableTrans )
      {
        if ( !getEnableInterTUACT() )
        {
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], pcTmpPredYuv, m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, m_ppcNoCorrYuv[uhDepth], (colourTransform? ACT_TRAN_CLR: ACT_ORG_CLR) DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[uhDepth]
#endif
                                                    );
        }
        else
        {
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], pcTmpPredYuv, m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, m_ppcNoCorrYuv[uhDepth], ACT_TWO_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[uhDepth]
#endif
                                                    );
        }
      }
      else
      {
#if PCC_RDO_EXT
        m_pcPredSearch->encodeResAndCalcRdInterCU(rpcTempCU, m_ppcOrigYuv[uhDepth], pcTmpPredYuv, m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, m_ppcNoCorrYuv[uhDepth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest), m_ppcOccupancyYuv[uhDepth]);
#else
        m_pcPredSearch->encodeResAndCalcRdInterCU(rpcTempCU, m_ppcOrigYuv[uhDepth], pcTmpPredYuv, m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, m_ppcNoCorrYuv[uhDepth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest));
#endif
      }
    }
    else
    {
      if ( m_pcEncCfg->getRGBFormatFlag() )
      {
        if ( !getEnableInterTUACT() )
        {
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], pcTmpPredYuv, m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, m_ppcNoCorrYuv[uhDepth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[uhDepth]
#endif
                                                    );
        }
        else
        {
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], pcTmpPredYuv, m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, m_ppcNoCorrYuv[uhDepth], ACT_TRAN_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[uhDepth]
#endif
                                                    );
        }
      }
      else
      {
        if ( !getEnableInterTUACT() )
        {
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], pcTmpPredYuv, m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, m_ppcNoCorrYuv[uhDepth], ACT_TRAN_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[uhDepth]
#endif
                                                    );
        }
        else
        {
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], pcTmpPredYuv, m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, m_ppcNoCorrYuv[uhDepth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[uhDepth]
#endif
                                                    );
        }
      }
    }
    rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

#if DEBUG_STRING
    DebugInterPredResiReco(sTest, *(m_ppcPredYuvTemp[uhDepth]), *(m_ppcResiYuvBest[uhDepth]), *(m_ppcRecoYuvTemp[uhDepth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif

    xCheckDQP( rpcTempCU );

    if ( bEnableTrans )
    {
      if ( !i )
      {
        Double rdCost1 = rpcTempCU->getTotalCost();
        if ( rdCost1 < m_ppcBestCU[uhDepth]->getTmpInterRDCost() )
        {
          m_ppcBestCU[uhDepth]->setTmpInterRDCost(rdCost1);
          m_ppcBestCU[uhDepth]->setInterCSCEnabled(true);
          m_ppcTempCU[uhDepth]->setTmpInterRDCost(rdCost1);
          m_ppcTempCU[uhDepth]->setInterCSCEnabled(true);
        }
      }
      else
      {
        Double rdCost = rpcTempCU->getTotalCost();
        if ( rdCost < m_ppcBestCU[uhDepth]->getTmpInterRDCost() )
        {
          m_ppcBestCU[uhDepth]->setTmpInterRDCost(rdCost);
          m_ppcBestCU[uhDepth]->setInterCSCEnabled(false);
          m_ppcTempCU[uhDepth]->setTmpInterRDCost(rdCost);
          m_ppcTempCU[uhDepth]->setInterCSCEnabled(false);
        }
      }
    }
    rpcTempCUPre = rpcTempCU;

    xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));

    if(bEnableTrans)
    {
      if( !rpcBestCU->isInter(0) || rpcBestCU->isSkipped(0) )
      {
        return;
      }
      Bool bParentUseCSC = ( (uhDepth != 0) && (m_ppcBestCU[uhDepth -1]->getInterCSCEnabled()) );
      if(!i && (!rpcTempCUPre->getQtRootCbf(0) || bParentUseCSC))
      {
        return;
      }
      else if(!i && rpcTempCUPre->getQtRootCbf(0))
      {
        if( rpcTempCU != rpcTempCUPre )
        {
          rpcTempCU->initEstData( uhDepth, orgQP, bTransquantBypassFlag);
          rpcTempCU->copyPartFrom( rpcBestCU, 0, uhDepth );
        }
      }
    }
    else
    {
      return;
    }
  }
}

Void TEncCu::xCheckRDCostIntra( TComDataCU *&rpcBestCU,
                                TComDataCU *&rpcTempCU,
                                Double      &cost,
                                PartSize     eSize
                                DEBUG_STRING_FN_DECLARE(sDebug),
                                Bool         bRGBIntraModeReuse
                               )
{
  DEBUG_STRING_NEW(sTest)

  if(getFastDeltaQp())
  {
    const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
    const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight()>>(sps.getLog2DiffMaxMinCodingBlockSize()), sps.getMaxCUHeight(), 32u);
    if(rpcTempCU->getWidth( 0 ) > fastDeltaQPCuMaxSize)
    {
      return; // only check necessary 2Nx2N Intra in fast deltaqp mode
    }
  }

  UInt uiDepth = rpcTempCU->getDepth( 0 );

  rpcTempCU->setSkipFlagSubParts( false, 0, uiDepth );

  rpcTempCU->setPartSizeSubParts( eSize, 0, uiDepth );
  rpcTempCU->setPredModeSubParts( MODE_INTRA, 0, uiDepth );
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uiDepth );
  rpcTempCU->setColourTransformSubParts(false, 0 ,uiDepth);

  Pel resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];

  if( bRGBIntraModeReuse )
  {
    m_pcPredSearch->estIntraPredLumaQTWithModeReuse( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], resiLuma
#if PCC_RDO_EXT
                                                    , m_ppcOccupancyYuv[uiDepth]
#endif
                                                    );
  }
  else
  {
#if PCC_RDO_EXT
    m_pcPredSearch->estIntraPredLumaQT(rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcOccupancyYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], resiLuma DEBUG_STRING_PASS_INTO(sTest));
#else
    m_pcPredSearch->estIntraPredLumaQT(rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], resiLuma DEBUG_STRING_PASS_INTO(sTest));
#endif
  }

  m_ppcRecoYuvTemp[uiDepth]->copyToPicComponent(COMPONENT_Y, rpcTempCU->getPic()->getPicYuvRec(), rpcTempCU->getCtuRsAddr(), rpcTempCU->getZorderIdxInCtu() );

  if (rpcBestCU->getPic()->getChromaFormat()!=CHROMA_400)
  {
    if( bRGBIntraModeReuse )
    {
      m_pcPredSearch->estIntraPredChromaQTWithModeReuse( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], resiLuma
#if PCC_RDO_EXT
                                                        , m_ppcOccupancyYuv[uiDepth]
#endif
                                                        );
    }
    else
    {
#if PCC_RDO_EXT
      m_pcPredSearch->estIntraPredChromaQT(rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcOccupancyYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], resiLuma DEBUG_STRING_PASS_INTO(sTest));
#else
      m_pcPredSearch->estIntraPredChromaQT(rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], resiLuma DEBUG_STRING_PASS_INTO(sTest));
#endif
    }
  }

  m_pcEntropyCoder->resetBits();

  if ( rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnabledFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( rpcTempCU, 0,          true );
  }

  m_pcEntropyCoder->encodeSkipFlag ( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePredMode( rpcTempCU, 0,          true );
  Bool bCodeDQP = getdQPFlag();
  Bool codeChromaQpAdjFlag = getCodeChromaQpAdjFlag();

  m_pcEntropyCoder->encodePaletteModeInfo( rpcTempCU, 0, true, &bCodeDQP, &codeChromaQpAdjFlag );

  m_pcEntropyCoder->encodePartSize( rpcTempCU, 0, uiDepth, true );
  m_pcEntropyCoder->encodePredInfo( rpcTempCU, 0 );
  m_pcEntropyCoder->encodeIPCMInfo(rpcTempCU, 0, true );

  // Encode Coefficients
  m_pcEntropyCoder->encodeCoeff( rpcTempCU, 0, uiDepth, bCodeDQP, codeChromaQpAdjFlag );
  setCodeChromaQpAdjFlag( codeChromaQpAdjFlag );
  setdQPFlag( bCodeDQP );

  m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);

  rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
  rpcTempCU->getTotalBins() = ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
  rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

  xCheckDQP( rpcTempCU );

  cost = rpcTempCU->getTotalCost();

  xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));
}


/** Check R-D costs for a CU with PCM mode.
 * \param rpcBestCU pointer to best mode CU data structure
 * \param rpcTempCU pointer to testing mode CU data structure
 * \returns Void
 *
 * \note Current PCM implementation encodes sample values in a lossless way. The distortion of PCM mode CUs are zero. PCM mode is selected if the best mode yields bits greater than that of PCM mode.
 */
Void TEncCu::xCheckIntraPCM( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU )
{
  if(getFastDeltaQp())
  {
    const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
    const UInt fastDeltaQPCuMaxPCMSize = Clip3((UInt)1<<sps.getPCMLog2MinSize(), (UInt)1<<sps.getPCMLog2MaxSize(), 32u);
    if (rpcTempCU->getWidth( 0 ) > fastDeltaQPCuMaxPCMSize)
    {
      return;   // only check necessary PCM in fast deltaqp mode
    }
  }
  
  UInt uiDepth = rpcTempCU->getDepth( 0 );

  rpcTempCU->setSkipFlagSubParts( false, 0, uiDepth );

  rpcTempCU->setIPCMFlag(0, true);
  rpcTempCU->setIPCMFlagSubParts (true, 0, rpcTempCU->getDepth(0));
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uiDepth );
  rpcTempCU->setPredModeSubParts( MODE_INTRA, 0, uiDepth );
  rpcTempCU->setTrIdxSubParts ( 0, 0, uiDepth );
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uiDepth );

  m_pcPredSearch->IPCMSearch( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth]);

  m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);

  m_pcEntropyCoder->resetBits();

  if ( rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnabledFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( rpcTempCU, 0,          true );
  }

  m_pcEntropyCoder->encodeSkipFlag ( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePredMode ( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePaletteModeInfo( rpcTempCU, 0,    true );
  m_pcEntropyCoder->encodePartSize ( rpcTempCU, 0, uiDepth, true );
  m_pcEntropyCoder->encodeIPCMInfo ( rpcTempCU, 0, true );

  m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);

  rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
  rpcTempCU->getTotalBins() = ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
  rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

  xCheckDQP( rpcTempCU );
  DEBUG_STRING_NEW(a)
  DEBUG_STRING_NEW(b)
  xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(a) DEBUG_STRING_PASS_INTO(b));
}

/** check whether current try is the best with identifying the depth of current try
 * \param rpcBestCU
 * \param rpcTempCU
 * \param uiDepth
 */
Void TEncCu::xCheckBestMode( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth DEBUG_STRING_FN_DECLARE(sParent) DEBUG_STRING_FN_DECLARE(sTest) DEBUG_STRING_PASS_INTO(Bool bAddSizeInfo) )
{
  if( rpcTempCU->getTotalCost() < rpcBestCU->getTotalCost() )
  {
    TComYuv* pcYuv;
    // Change Information data
    TComDataCU* pcCU = rpcBestCU;
    rpcBestCU = rpcTempCU;
    rpcTempCU = pcCU;

    // Change Prediction data
    pcYuv = m_ppcPredYuvBest[uiDepth];
    m_ppcPredYuvBest[uiDepth] = m_ppcPredYuvTemp[uiDepth];
    m_ppcPredYuvTemp[uiDepth] = pcYuv;

    // Change Reconstruction data
    pcYuv = m_ppcRecoYuvBest[uiDepth];
    m_ppcRecoYuvBest[uiDepth] = m_ppcRecoYuvTemp[uiDepth];
    m_ppcRecoYuvTemp[uiDepth] = pcYuv;

    pcYuv = NULL;
    pcCU  = NULL;

    // store temp best CI for next CU coding
    m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]->store(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);


#if DEBUG_STRING
    DEBUG_STRING_SWAP(sParent, sTest)
    const PredMode predMode=rpcBestCU->getPredictionMode(0);
    if ((DebugOptionList::DebugString_Structure.getInt()&DebugStringGetPredModeMask(predMode)) && bAddSizeInfo)
    {
      std::stringstream ss(stringstream::out);
      ss <<"###: " << (predMode==MODE_INTRA?"Intra   ":"Inter   ") << partSizeToString[rpcBestCU->getPartitionSize(0)] << " CU at " << rpcBestCU->getCUPelX() << ", " << rpcBestCU->getCUPelY() << " width=" << UInt(rpcBestCU->getWidth(0)) << std::endl;
      sParent+=ss.str();
    }
#endif
  }
}

Void TEncCu::xCheckDQP( TComDataCU* pcCU )
{
  UInt uiDepth = pcCU->getDepth( 0 );

  const TComPPS &pps = *(pcCU->getSlice()->getPPS());
  if ( pps.getUseDQP() && uiDepth <= pps.getMaxCuDQPDepth() )
  {
    if ( pcCU->getQtRootCbf( 0) || ( pcCU->getPaletteModeFlag(0) && pcCU->getPaletteEscape(COMPONENT_Y, 0) ) )
    {
      m_pcEntropyCoder->resetBits();
      m_pcEntropyCoder->encodeQP( pcCU, 0, false );
      pcCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
      pcCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
      pcCU->getTotalCost() = m_pcRdCost->calcRdCost( pcCU->getTotalBits(), pcCU->getTotalDistortion() );
    }
    else
    {
      pcCU->setQPSubParts( pcCU->getRefQP( 0 ), 0, uiDepth ); // set QP to default QP
    }
  }
}

Void TEncCu::xCopyAMVPInfo (AMVPInfo* pSrc, AMVPInfo* pDst)
{
  pDst->iN = pSrc->iN;
  for (Int i = 0; i < pSrc->iN; i++)
  {
    pDst->m_acMvCand[i] = pSrc->m_acMvCand[i];
  }
}
Void TEncCu::xCopyYuv2Pic(TComPic* rpcPic, UInt uiCUAddr, UInt uiAbsPartIdx, UInt uiDepth, UInt uiSrcDepth, TComDataCU* pcCU )
{
  UInt uiAbsPartIdxInRaster = g_auiZscanToRaster[uiAbsPartIdx];
  UInt uiSrcBlkWidth = rpcPic->getNumPartInCtuWidth() >> (uiSrcDepth);
  UInt uiBlkWidth    = rpcPic->getNumPartInCtuWidth() >> (uiDepth);
  UInt uiPartIdxX = ( ( uiAbsPartIdxInRaster % rpcPic->getNumPartInCtuWidth() ) % uiSrcBlkWidth) / uiBlkWidth;
  UInt uiPartIdxY = ( ( uiAbsPartIdxInRaster / rpcPic->getNumPartInCtuWidth() ) % uiSrcBlkWidth) / uiBlkWidth;
  UInt uiPartIdx = uiPartIdxY * ( uiSrcBlkWidth / uiBlkWidth ) + uiPartIdxX;
  m_ppcRecoYuvBest[uiSrcDepth]->copyToPicYuv( rpcPic->getPicYuvRec (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);
  if( pcCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() && m_pcEncCfg->getRGBFormatFlag() )
  {
    TComRectangle cuCompRect;
    cuCompRect.x0     = 0;
    cuCompRect.y0     = 0;
    cuCompRect.width  = (pcCU->getSlice()->getSPS()->getMaxCUWidth()>>uiSrcDepth);
    cuCompRect.height = (pcCU->getSlice()->getSPS()->getMaxCUHeight()>>uiSrcDepth);

    for ( UInt ch = 0; ch < MAX_NUM_COMPONENT; ch++ )
    {
      m_ppcRecoYuvBest[uiSrcDepth]->copyPartToPartComponentMxN( ComponentID( ch ), m_ppcRecoYuvTemp[uiSrcDepth], cuCompRect );
    }
    m_ppcRecoYuvTemp[uiSrcDepth]->DefaultConvertPix( cuCompRect.x0, cuCompRect.y0, cuCompRect.width, pcCU->getSlice()->getSPS()->getBitDepths() );
    m_ppcRecoYuvTemp[uiSrcDepth]->copyToPicYuv( rpcPic->getPicYuvCSC(), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);
  }

  m_ppcPredYuvBest[uiSrcDepth]->copyToPicYuv( rpcPic->getPicYuvPred (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);
}

Void TEncCu::xCopyYuv2Tmp( UInt uiPartUnitIdx, UInt uiNextDepth )
{
  UInt uiCurrDepth = uiNextDepth - 1;
  m_ppcRecoYuvBest[uiNextDepth]->copyToPartYuv( m_ppcRecoYuvTemp[uiCurrDepth], uiPartUnitIdx );
  m_ppcPredYuvBest[uiNextDepth]->copyToPartYuv( m_ppcPredYuvBest[uiCurrDepth], uiPartUnitIdx);
}

/** Function for filling the PCM buffer of a CU using its original sample array
 * \param pCU pointer to current CU
 * \param pOrgYuv pointer to original sample array
 */
Void TEncCu::xFillPCMBuffer     ( TComDataCU* pCU, TComYuv* pOrgYuv )
{
  const ChromaFormat format = pCU->getPic()->getChromaFormat();
  const UInt numberValidComponents = getNumberValidComponents(format);
  for (UInt componentIndex = 0; componentIndex < numberValidComponents; componentIndex++)
  {
    const ComponentID component = ComponentID(componentIndex);

    const UInt width  = pCU->getWidth(0)  >> getComponentScaleX(component, format);
    const UInt height = pCU->getHeight(0) >> getComponentScaleY(component, format);

    Pel *source      = pOrgYuv->getAddr(component, 0, width);
    Pel *destination = pCU->getPCMSample(component);

    const UInt sourceStride = pOrgYuv->getStride(component);

    for (Int line = 0; line < height; line++)
    {
      for (Int column = 0; column < width; column++)
      {
        destination[column] = source[column];
      }

      source      += sourceStride;
      destination += width;
    }
  }
}

#if ADAPTIVE_QP_SELECTION
/** Collect ARL statistics from one block
  */
Int TEncCu::xTuCollectARLStats(TCoeff* rpcCoeff, TCoeff* rpcArlCoeff, Int NumCoeffInCU, Double* cSum, UInt* numSamples )
{
  for( Int n = 0; n < NumCoeffInCU; n++ )
  {
    TCoeff u = abs( rpcCoeff[ n ] );
    TCoeff absc = rpcArlCoeff[ n ];

    if( u != 0 )
    {
      if( u < LEVEL_RANGE )
      {
        cSum[ u ] += ( Double )absc;
        numSamples[ u ]++;
      }
      else
      {
        cSum[ LEVEL_RANGE ] += ( Double )absc - ( Double )( u << ARL_C_PRECISION );
        numSamples[ LEVEL_RANGE ]++;
      }
    }
  }

  return 0;
}

//! Collect ARL statistics from one CTU
Void TEncCu::xCtuCollectARLStats(TComDataCU* pCtu )
{
  Double cSum[ LEVEL_RANGE + 1 ];     //: the sum of DCT coefficients corresponding to data type and quantization output
  UInt numSamples[ LEVEL_RANGE + 1 ]; //: the number of coefficients corresponding to data type and quantization output

  TCoeff* pCoeffY = pCtu->getCoeff(COMPONENT_Y);
  TCoeff* pArlCoeffY = pCtu->getArlCoeff(COMPONENT_Y);
  const TComSPS &sps = *(pCtu->getSlice()->getSPS());

  const UInt uiMinCUWidth = sps.getMaxCUWidth() >> sps.getMaxTotalCUDepth(); // NOTE: ed - this is not the minimum CU width. It is the square-root of the number of coefficients per part.
  const UInt uiMinNumCoeffInCU = 1 << uiMinCUWidth;                          // NOTE: ed - what is this?

  memset( cSum, 0, sizeof( Double )*(LEVEL_RANGE+1) );
  memset( numSamples, 0, sizeof( UInt )*(LEVEL_RANGE+1) );

  // Collect stats to cSum[][] and numSamples[][]
  for(Int i = 0; i < pCtu->getTotalNumPart(); i ++ )
  {
    UInt uiTrIdx = pCtu->getTransformIdx(i);

    if(pCtu->isInter(i) && pCtu->getCbf( i, COMPONENT_Y, uiTrIdx ) )
    {
      xTuCollectARLStats(pCoeffY, pArlCoeffY, uiMinNumCoeffInCU, cSum, numSamples);
    }//Note that only InterY is processed. QP rounding is based on InterY data only.

    pCoeffY  += uiMinNumCoeffInCU;
    pArlCoeffY  += uiMinNumCoeffInCU;
  }

  for(Int u=1; u<LEVEL_RANGE;u++)
  {
    m_pcTrQuant->getSliceSumC()[u] += cSum[ u ] ;
    m_pcTrQuant->getSliceNSamples()[u] += numSamples[ u ] ;
  }
  m_pcTrQuant->getSliceSumC()[LEVEL_RANGE] += cSum[ LEVEL_RANGE ] ;
  m_pcTrQuant->getSliceNSamples()[LEVEL_RANGE] += numSamples[ LEVEL_RANGE ] ;
}
#endif

// SCM new added functions
// ====================================================================================================================
// Static function
// ====================================================================================================================
/** Calculates the minimum of the H and V luma activity of a CU in the source image.
*\param   pcCU
*\param   uiAbsPartIdx
*\param   ppcOrigYuv
*\returns Intermediate_Int
*
*/
static Intermediate_Int
CalculateMinimumHVLumaActivity(TComDataCU *pcCU, const UInt uiAbsPartIdx, const TComYuv * const * const ppcOrigYuv)
{
  const TComYuv *pOrgYuv = ppcOrigYuv[pcCU->getDepth(uiAbsPartIdx)];
  const Int      stride  = pOrgYuv->getStride(COMPONENT_Y);
  const Int      width   = pcCU->getWidth(uiAbsPartIdx);
  const Int      height  = pcCU->getHeight(uiAbsPartIdx);

  // Get activity
  Intermediate_Int hAct = 0;
  const Pel *pY = pOrgYuv->getAddr(COMPONENT_Y, uiAbsPartIdx);
  for (Int y = 0; y < height; y++)
  {
    for (Int x = 1; x < width; x++)
    {
      hAct += abs( pY[x] - pY[x - 1]);
    }
    pY += stride;
  }

  Intermediate_Int vAct = 0;
  pY = pOrgYuv->getAddr(COMPONENT_Y, 0) + stride;
  for (Int y = 1; y < height; y++)
  {
    for (Int x = 0; x < width; x++)
    {
      vAct += abs(pY[x] - pY[x - stride]);
    }
    pY += stride;
  }

  return std::min(hAct, vAct);
}

Void TEncCu::setEnableIntraTUACT(UInt depth, TComSlice* pcSlice)
{
  if( m_pcEncCfg->getRGBFormatFlag() )
  {
    m_bEnableIntraTUACTRD = true;
    if( (m_pcEncCfg->getGOPSize() == 8 && pcSlice->getDepth() == 3) || (m_pcEncCfg->getGOPSize() == 4 && pcSlice->getDepth() != 2) || (depth < 3) )
    {
      m_bEnableIntraTUACTRD = false;
    }
  }
  else
  {
    m_bEnableIntraTUACTRD = false;
  }
}

Void TEncCu::setEnableIBCTUACT(UInt depth, TComSlice* pcSlice)
{
  if( m_pcEncCfg->getRGBFormatFlag() )
  {
    m_bEnableIBCTUACTRD = true;
    if( (m_pcEncCfg->getGOPSize() == 8 && pcSlice->getDepth() == 3) || (m_pcEncCfg->getGOPSize() == 4 && pcSlice->getDepth() != 2) || (depth < 3) )
    {
      m_bEnableIBCTUACTRD = false;
    }
  }
  else
  {
    m_bEnableIBCTUACTRD = false;
  }
}

Void TEncCu::setEnableInterTUACT(UInt depth, TComSlice* pcSlice)
{
  if( m_pcEncCfg->getRGBFormatFlag() )
  {
    m_bEnableInterTUACTRD = true;
    if( (pcSlice->getDepth() == 0) || (depth < 2) )
    {
      m_bEnableInterTUACTRD = false;
    }
  }
  else
  {
    m_bEnableInterTUACTRD = false;
  }
}

Void TEncCu::xCheckRDCostIntraBCMerge2Nx2N( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU )
{
  TComMvField  cMvFieldNeighbours[2 * MRG_MAX_NUM_CANDS]; // double length for mv of both lists
  UChar interDirNeighbours[MRG_MAX_NUM_CANDS];
  Int numValidMergeCand = 0;
  const Bool bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);

  for( UInt ui = 0; ui < rpcTempCU->getSlice()->getMaxNumMergeCand(); ++ui )
  {
    interDirNeighbours[ui] = 0;
  }

  Int xPos = rpcTempCU->getCUPelX();
  Int yPos = rpcTempCU->getCUPelY();
  Int width = rpcTempCU->getWidth( 0 );
  Int height = rpcTempCU->getHeight( 0 );
  UChar depth = rpcTempCU->getDepth( 0 );
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, depth );
#if MCTS_ENC_CHECK
  UInt numSpatialMergeCandidates = 0;
  rpcTempCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours,interDirNeighbours, numValidMergeCand, numSpatialMergeCandidates );
#else
  rpcTempCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours,interDirNeighbours, numValidMergeCand );
#endif
  rpcTempCU->roundMergeCandidates(cMvFieldNeighbours, numValidMergeCand);
  rpcTempCU->xRestrictBipredMergeCand( 0, cMvFieldNeighbours, interDirNeighbours, numValidMergeCand );

  Int mergeCandBuffer[MRG_MAX_NUM_CANDS];
  for( UInt ui = 0; ui < numValidMergeCand; ++ui )
  {
    mergeCandBuffer[ui] = 0;
  }

#if MCTS_ENC_CHECK
  if (m_pcEncCfg->getTMCTSSEITileConstraint() && rpcTempCU->isLastColumnCTUInTile())
  {
    numValidMergeCand = numSpatialMergeCandidates;
  }
#endif

  Bool bestIsSkip = false;

  UInt iteration;
  if ( rpcTempCU->isLosslessCoded(0))
  {
    iteration = 1;
  }
  else
  {
    iteration = 2;
  }

  const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
  for( UInt noResidual = 0; noResidual < iteration; ++noResidual )
  {
    for( UInt mergeCand = 0; mergeCand < numValidMergeCand; ++mergeCand )
    {
      if(!(noResidual==1 && mergeCandBuffer[mergeCand]==1))
      {
        if( !(bestIsSkip && noResidual == 0) )
        {
          if ( interDirNeighbours[mergeCand] != 1 )
          {
            continue;
          }

          if ( rpcTempCU->getSlice()->getRefPic( REF_PIC_LIST_0, cMvFieldNeighbours[mergeCand<<1].getRefIdx() )->getPOC() != rpcTempCU->getSlice()->getPOC() )
          {
            continue;
          }

          if ( !m_pcPredSearch->isBlockVectorValid( xPos, yPos, width, height, rpcTempCU,
            0, 0, (cMvFieldNeighbours[mergeCand<<1].getHor() >> 2), (cMvFieldNeighbours[mergeCand<<1].getVer()>>2), sps.getMaxCUWidth() ) )
          {
            continue;
          }

          // set MC parameters
          rpcTempCU->setPredModeSubParts( MODE_INTER, 0, depth ); // interprets depth relative to LCU level
          rpcTempCU->setCUTransquantBypassSubParts( bTransquantBypassFlag, 0, depth );
          rpcTempCU->setChromaQpAdjSubParts( bTransquantBypassFlag ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, depth );
          rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, depth ); // interprets depth relative to LCU level
          rpcTempCU->setMergeFlagSubParts( true, 0, 0, depth ); // interprets depth relative to LCU level
          rpcTempCU->setMergeIndexSubParts( mergeCand, 0, 0, depth ); // interprets depth relative to LCU level
          rpcTempCU->setInterDirSubParts( interDirNeighbours[mergeCand], 0, 0, depth ); // interprets depth relative to LCU level
          rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->setAllMvField( cMvFieldNeighbours[0 + 2*mergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level
          rpcTempCU->getCUMvField( REF_PIC_LIST_1 )->setAllMvField( cMvFieldNeighbours[1 + 2*mergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level

#if MCTS_ENC_CHECK
          if ( m_pcEncCfg->getTMCTSSEITileConstraint () && (!(m_pcPredSearch->checkTMctsMvp(rpcTempCU))))
          {
            continue;
          }
#endif

          // do MC
          m_pcPredSearch->motionCompensation ( rpcTempCU, m_ppcPredYuvTemp[depth] );
          Bool bColourTrans = (m_pcEncCfg->getRGBFormatFlag() && rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans())? true : false;

          if( bTransquantBypassFlag && (sps.getBitDepth( CHANNEL_TYPE_LUMA ) != sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
          {
            bColourTrans = false;
          }

          TComYuv* pcTmpPredYuv = m_ppcPredYuvTemp[depth];
          if( rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() && (!bTransquantBypassFlag || sps.getBitDepth( CHANNEL_TYPE_LUMA ) == sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
          {
            if(!getEnableIBCTUACT())
            {
              m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                         m_ppcOrigYuv    [depth],
                                                         pcTmpPredYuv           ,
                                                         m_ppcResiYuvTemp[depth],
                                                         m_ppcResiYuvBest[depth],
                                                         m_ppcRecoYuvTemp[depth],
                                                         (noResidual != 0),
                                                         m_ppcNoCorrYuv  [depth],
                                                         (bColourTrans? ACT_TRAN_CLR: ACT_ORG_CLR)
                                                         DEBUG_STRING_PASS_INTO(tmpStr)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                        );
            }
            else
            {
              m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                         m_ppcOrigYuv    [depth],
                                                         pcTmpPredYuv           ,
                                                         m_ppcResiYuvTemp[depth],
                                                         m_ppcResiYuvBest[depth],
                                                         m_ppcRecoYuvTemp[depth],
                                                         (noResidual != 0),
                                                         m_ppcNoCorrYuv  [depth],
                                                         ACT_TWO_CLR
                                                         DEBUG_STRING_PASS_INTO(tmpStr)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                        );
            }
          }
          else
          {
#if PCC_RDO_EXT
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                       m_ppcOrigYuv    [depth],
                                                       pcTmpPredYuv           ,
                                                       m_ppcResiYuvTemp[depth],
                                                       m_ppcResiYuvBest[depth],
                                                       m_ppcRecoYuvTemp[depth],
                                                       (noResidual != 0),
                                                       m_ppcNoCorrYuv  [depth],
                                                       ACT_ORG_CLR
                                                       DEBUG_STRING_PASS_INTO(tmpStr),
                                                       m_ppcOccupancyYuv[depth]);
#else
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                       m_ppcOrigYuv    [depth],
                                                       pcTmpPredYuv           ,
                                                       m_ppcResiYuvTemp[depth],
                                                       m_ppcResiYuvBest[depth],
                                                       m_ppcRecoYuvTemp[depth],
                                                       (noResidual != 0),
                                                       m_ppcNoCorrYuv  [depth],
                                                       ACT_ORG_CLR
                                                       DEBUG_STRING_PASS_INTO(tmpStr) );
#endif
          }

          if ((noResidual == 0) && (rpcTempCU->getQtRootCbf(0) == 0))
          {
            // If no residual when allowing for one, then set mark to not try case where residual is forced to 0
            mergeCandBuffer[mergeCand] = 1;
          }

          rpcTempCU->setSkipFlagSubParts( rpcTempCU->getQtRootCbf(0) == 0, 0, depth );
          Int orgQP = rpcTempCU->getQP( 0 );
          xCheckDQP( rpcTempCU );
          TComDataCU *rpcTempCUPre = rpcTempCU;
          xCheckBestMode(rpcBestCU, rpcTempCU, depth DEBUG_STRING_PASS_INTO(bestStr) DEBUG_STRING_PASS_INTO(tmpStr));

          if(rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() && !noResidual && rpcTempCUPre->getQtRootCbf(0) )
          {
            if( rpcTempCU != rpcTempCUPre )
            {
              rpcTempCU->initEstData( depth, orgQP, bTransquantBypassFlag );
              rpcTempCU->copyPartFrom( rpcBestCU, 0, depth );
            }
            if(!getEnableIBCTUACT())
            {
              bColourTrans = !bColourTrans;
            }
            else
            {
              bColourTrans = m_pcEncCfg->getRGBFormatFlag();
            }
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                       m_ppcOrigYuv    [depth],
                                                       pcTmpPredYuv           ,
                                                       m_ppcResiYuvTemp[depth],
                                                       m_ppcResiYuvBest[depth],
                                                       m_ppcRecoYuvTemp[depth],
                                                       (noResidual != 0),
                                                       m_ppcNoCorrYuv  [depth],
                                                       (bColourTrans? ACT_TRAN_CLR: ACT_ORG_CLR)
                                                       DEBUG_STRING_PASS_INTO(tmpStr)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
            rpcTempCU->setSkipFlagSubParts( rpcTempCU->getQtRootCbf(0) == 0, 0, depth );
            //Double rdCost = rpcTempCU->getTotalCost();

            xCheckDQP( rpcTempCU );
            if(rpcTempCU->getQtRootCbf(0))
            {
              xCheckBestMode( rpcBestCU, rpcTempCU, depth );
            }
          }
          rpcTempCU->initEstData( depth, orgQP, bTransquantBypassFlag );

          if( m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip )
          {
            bestIsSkip = rpcBestCU->getQtRootCbf(0) == 0;
          }
        }
      }
    }
  }
}

Void TEncCu::xCheckRDCostIntraCSC( TComDataCU     *&rpcBestCU,
                                   TComDataCU     *&rpcTempCU,
                                    Double         &cost,
                                   PartSize       eSize,
                                   ACTRDTestTypes eACTRDTestType,
                                   Bool           bReuseIntraMode
                                   DEBUG_STRING_FN_DECLARE(sDebug)
)
{
  DEBUG_STRING_NEW(sTest)

  UInt depth = rpcTempCU->getDepth( 0 );
  rpcTempCU->setSkipFlagSubParts( false, 0, depth );
  rpcTempCU->setPartSizeSubParts( eSize, 0, depth );
  rpcTempCU->setPredModeSubParts( MODE_INTRA, 0, depth );
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, depth );

  m_pcPredSearch->estIntraPredQTCT( rpcTempCU, m_ppcOrigYuv[depth], m_ppcPredYuvTemp[depth], m_ppcResiYuvTemp[depth], m_ppcRecoYuvTemp[depth], eACTRDTestType, bReuseIntraMode DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                              , m_ppcOccupancyYuv[depth]
#endif
                                   );

  m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[depth][CI_CURR_BEST]);
  m_pcEntropyCoder->resetBits();

  if ( rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnabledFlag() )
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( rpcTempCU, 0, true );
  }

  m_pcEntropyCoder->encodeSkipFlag ( rpcTempCU, 0, true );

  m_pcEntropyCoder->encodePredMode( rpcTempCU, 0, true );

  Bool bCodeDQP = getdQPFlag();
  Bool codeChromaQpAdjFlag = getCodeChromaQpAdjFlag();

  m_pcEntropyCoder->encodePaletteModeInfo( rpcTempCU, 0, true, &bCodeDQP, &codeChromaQpAdjFlag );

  m_pcEntropyCoder->encodePartSize( rpcTempCU, 0, depth, true );
  m_pcEntropyCoder->encodePredInfo( rpcTempCU, 0 );
  m_pcEntropyCoder->encodeIPCMInfo(rpcTempCU, 0, true );

  // Encode Coefficients
  m_pcEntropyCoder->encodeCoeff( rpcTempCU, 0, depth, bCodeDQP, codeChromaQpAdjFlag );
  setCodeChromaQpAdjFlag( codeChromaQpAdjFlag );
  setdQPFlag( bCodeDQP );

  m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[depth][CI_TEMP_BEST]);

  rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
  rpcTempCU->getTotalBins() = ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
  rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

  xCheckDQP( rpcTempCU );

  cost = rpcTempCU->getTotalCost();

  xCheckBestMode(rpcBestCU, rpcTempCU, depth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));
}

Void TEncCu::xCheckRDCostIntraBC( TComDataCU *&rpcBestCU,
                                  TComDataCU *&rpcTempCU,
                                  Bool         bUse1DSearchFor8x8,
                                  PartSize     eSize,
                                  Double      &rdCost,
                                  Bool         testPredOnly,
                                  TComMv      *iMVCandList
                                  DEBUG_STRING_FN_DECLARE(sDebug))
{
  DEBUG_STRING_NEW(sTest)
  UInt depth = rpcTempCU->getDepth( 0 );

  rpcTempCU->setDepthSubParts( depth, 0 );
  rpcTempCU->setSkipFlagSubParts( false, 0, depth );
  rpcTempCU->setPartSizeSubParts( eSize, 0, depth );
  rpcTempCU->setPredModeSubParts( MODE_INTER, 0, depth );

  rpcTempCU->setIntraDirSubParts( CHANNEL_TYPE_LUMA, DC_IDX, 0, depth );
  rpcTempCU->setIntraDirSubParts( CHANNEL_TYPE_CHROMA, DC_IDX, 0, depth );
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, depth );

#if MCTS_ENC_CHECK
  rpcTempCU->setTMctsMvpIsValid(true);
#endif

  // intra BV search
  const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
  if( rpcTempCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans () && m_pcEncCfg->getRGBFormatFlag() && (!rpcTempCU->getCUTransquantBypass(0) || (sps.getBitDepth( CHANNEL_TYPE_LUMA ) == sps.getBitDepth( CHANNEL_TYPE_CHROMA ))) )
  {
    m_ppcOrigYuv[depth]->copyFromPicYuv( rpcTempCU->getPic()->getPicYuvResi(), rpcTempCU->getCtuRsAddr(), rpcTempCU->getZorderIdxInCtu() );
    rpcTempCU->getPic()->exchangePicYuvRec();  //YCgCo reference picture
  }

  Bool bValid = m_pcPredSearch->predIntraBCSearch ( rpcTempCU,
                                                    m_ppcOrigYuv[depth],
                                                    m_ppcPredYuvTemp[depth],
                                                    m_ppcResiYuvTemp[depth],
                                                    m_ppcRecoYuvTemp[depth]
                                                    DEBUG_STRING_PASS_INTO(sTest),
                                                    bUse1DSearchFor8x8,
                                                    false,
                                                    testPredOnly
                                                  );

  if ( bValid && (rpcTempCU->getWidth( 0 ) <= 16) && (eSize == SIZE_2NxN || eSize == SIZE_Nx2N) )
  {
    Int dummyWidth, dummyHeight;
    UInt partAddr = 0;
    rpcTempCU->getPartIndexAndSize( 1, partAddr, dummyWidth, dummyHeight );
    iMVCandList[0] = rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->getMv( 0 );
    iMVCandList[1] = rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->getMv( partAddr );
  }
  if( rpcTempCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans() && m_pcEncCfg->getRGBFormatFlag() && (!rpcTempCU->getCUTransquantBypass(0) || (sps.getBitDepth( CHANNEL_TYPE_LUMA ) == sps.getBitDepth( CHANNEL_TYPE_CHROMA ))) )
  {
    m_ppcOrigYuv[depth]->copyFromPicYuv( rpcTempCU->getPic()->getPicYuvOrg(), rpcTempCU->getCtuRsAddr(), rpcTempCU->getZorderIdxInCtu() );
    rpcTempCU->getPic()->exchangePicYuvRec();
  }
  
  if (bValid)
  {
#if MCTS_ENC_CHECK
    if (m_pcEncCfg->getTMCTSSEITileConstraint() && (!rpcTempCU->getTMctsMvpIsValid()))
    {
      return;
    }
#endif
    TComDataCU *rpcTempCUPre = NULL;
    Bool   bEnableTrans      = rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans();
    UChar  colourTransform = 0;
    Bool   bRGB              = m_pcEncCfg->getRGBFormatFlag();
    Double dCostFst          = MAX_DOUBLE;
    SChar   orgQP            = rpcTempCU->getQP( 0 );
    Bool   bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);
    if ( bTransquantBypassFlag && (sps.getBitDepth( CHANNEL_TYPE_LUMA ) != sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
    {
      bEnableTrans = false;
    }
    TComYuv* pcTmpPredYuv = m_ppcPredYuvTemp[depth];
    for(UInt i = 0;  i < (testPredOnly ? 1 : 2) ; i++)
    {
      colourTransform = ( bRGB && bEnableTrans )? (1-i): i;

      if( colourTransform && m_pcEncCfg->getRGBFormatFlag())
      {
        m_pcPredSearch->motionCompensation ( rpcTempCU, m_ppcPredYuvTemp[depth] );
      }

      if ( i == 0 )
      {
        if ( bEnableTrans )
        {
          if ( !getEnableIBCTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], (colourTransform? ACT_TRAN_CLR: ACT_ORG_CLR) DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                      , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TWO_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                      , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
        else
        {
#if PCC_RDO_EXT
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest), m_ppcOccupancyYuv[depth] );
#else
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest) );
#endif
        }
      }
      else
      {
        if ( m_pcEncCfg->getRGBFormatFlag() )
        {
          if ( !getEnableIBCTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TRAN_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
        else
        {
          if ( !getEnableIBCTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TRAN_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
      }
      rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );
      rdCost = rpcTempCU->getTotalCost();

#if DEBUG_STRING
      DebugInterPredResiReco(sTest, *(m_ppcPredYuvTemp[depth]), *(m_ppcResiYuvBest[depth]), *(m_ppcRecoYuvTemp[depth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif

      if ( bEnableTrans )
      {
        if ( !i )
        {
          if ( rdCost < m_ppcBestCU[depth]->getTmpIntraBCRDCost() )
          {
            m_ppcBestCU[depth]->setTmpIntraBCRDCost(rdCost);
            m_ppcBestCU[depth]->setIntraBCCSCEnabled(true);
            m_ppcTempCU[depth]->setTmpIntraBCRDCost(rdCost);
            m_ppcTempCU[depth]->setIntraBCCSCEnabled(true);
          }
        }
        else
        {
          if ( rdCost < m_ppcBestCU[depth]->getTmpIntraBCRDCost() )
          {
            m_ppcBestCU[depth]->setTmpIntraBCRDCost(rdCost);
            m_ppcBestCU[depth]->setIntraBCCSCEnabled(false);
            m_ppcTempCU[depth]->setTmpIntraBCRDCost(rdCost);
            m_ppcTempCU[depth]->setIntraBCCSCEnabled(false);
          }
        }
      }
      xCheckDQP( rpcTempCU );
      rpcTempCUPre = rpcTempCU;
      xCheckBestMode(rpcBestCU, rpcTempCU, depth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));
      if( bEnableTrans )
      {
        if( testPredOnly || !rpcBestCU->isIntraBC(0) )
        {
          return;
        }
        Bool bParentUseCSC = ( ((depth > 2)) && (m_ppcBestCU[depth -1]->getIntraBCCSCEnabled()) );
        if(!i && (!rpcTempCUPre->getQtRootCbf(0) || bParentUseCSC))
        {
          return;
        }
        else if(!i)
        {
          dCostFst = rdCost;
          if( rpcTempCU != rpcTempCUPre )
          {
            rpcTempCU->initEstData( depth, orgQP, bTransquantBypassFlag );
            rpcTempCU->copyPartFrom( rpcBestCU, 0, depth );
          }
        }
        else
        {
          if ( rpcTempCU == rpcTempCUPre )
          {
            rdCost = dCostFst;
          }
        }
      }
      else
      {
        return;
      }
    }
  }
  else
  {
    rdCost = MAX_DOUBLE;
  }
}

Void TEncCu::xCheckRDCostIntraBCMixed( TComDataCU *&rpcBestCU,
                                       TComDataCU *&rpcTempCU,
                                       PartSize     eSize,
                                       Double      &rdCost
                                       DEBUG_STRING_FN_DECLARE( sDebug ),
                                       TComMv* iMVCandList
                                     )
{
  DEBUG_STRING_NEW( sTest )
  UInt depth = rpcTempCU->getDepth( 0 );
  rpcTempCU->setPredModeSubParts( MODE_INTER, 0, depth );
  rpcTempCU->setDepthSubParts( depth, 0 );
  rpcTempCU->setSkipFlagSubParts( false, 0, depth );
  rpcTempCU->setPartSizeSubParts( eSize, 0, depth );
  rpcTempCU->setIntraDirSubParts( CHANNEL_TYPE_LUMA, DC_IDX, 0, depth );
  rpcTempCU->setIntraDirSubParts( CHANNEL_TYPE_CHROMA, DC_IDX, 0, depth );
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass( 0 ) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, depth );

#if MCTS_ENC_CHECK
  rpcTempCU->setTMctsMvpIsValid(true);
#endif

  Bool bValid = m_pcPredSearch->predMixedIntraBCInterSearch( rpcTempCU,
                                                             m_ppcOrigYuv[depth],
                                                             m_ppcPredYuvTemp[depth],
                                                             m_ppcResiYuvTemp[depth],
                                                             m_ppcRecoYuvTemp[depth]
                                                             DEBUG_STRING_PASS_INTO( sTest ),
                                                             iMVCandList,
                                                             false
                                                           );

  if ( bValid )
  {
#if MCTS_ENC_CHECK
    if (m_pcEncCfg->getTMCTSSEITileConstraint() && (!rpcTempCU->getTMctsMvpIsValid()))
    {
      return;
    }
#endif
    TComDataCU *rpcTempCUPre = NULL;
    Bool   bEnableTrans      = rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans();
    UChar  colourTransform = 0;
    Bool   bRGB              = m_pcEncCfg->getRGBFormatFlag();
    Double dCostFst          = MAX_DOUBLE;
    SChar   orgQP            = rpcTempCU->getQP( 0 );
    Bool   bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);
    const TComSPS &sps       = *(rpcTempCU->getSlice()->getSPS());
    if ( bTransquantBypassFlag && (sps.getBitDepth( CHANNEL_TYPE_LUMA ) != sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
    {
      bEnableTrans = false;
    }
    TComYuv* pcTmpPredYuv = m_ppcPredYuvTemp[depth];

    for(UInt i = 0;  i < 2 ; i++)
    {
      colourTransform = (bRGB && bEnableTrans)? (1-i): i;
      if ( i == 0 )
      {
        if ( bEnableTrans )
        {
          if ( !getEnableInterTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], (colourTransform? ACT_TRAN_CLR: ACT_ORG_CLR) DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TWO_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
        else
        {
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
              , m_ppcOccupancyYuv[depth]
#endif
              );
        }
      }
      else
      {
        if ( m_pcEncCfg->getRGBFormatFlag() )
        {
          if ( !getEnableInterTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TRAN_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
        else
        {
          if ( !getEnableInterTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TRAN_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
      }
      rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );
      rdCost = rpcTempCU->getTotalCost();

#if DEBUG_STRING
      DebugInterPredResiReco(sTest, *(m_ppcPredYuvTemp[depth]), *(m_ppcResiYuvBest[depth]), *(m_ppcRecoYuvTemp[depth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif
      if ( bEnableTrans )
      {
        if ( !i )
        {
          if ( rdCost < m_ppcBestCU[depth]->getTmpIntraBCRDCost() )
          {
            m_ppcBestCU[depth]->setTmpIntraBCRDCost(rdCost);
            m_ppcBestCU[depth]->setIntraBCCSCEnabled(true);
            m_ppcTempCU[depth]->setTmpIntraBCRDCost(rdCost);
            m_ppcTempCU[depth]->setIntraBCCSCEnabled(true);
          }
        }
        else
        {
          if ( rdCost < m_ppcBestCU[depth]->getTmpIntraBCRDCost() )
          {
            m_ppcBestCU[depth]->setTmpIntraBCRDCost(rdCost);
            m_ppcBestCU[depth]->setIntraBCCSCEnabled(false);
            m_ppcTempCU[depth]->setTmpIntraBCRDCost(rdCost);
            m_ppcTempCU[depth]->setIntraBCCSCEnabled(false);
          }
        }
      }
      xCheckDQP( rpcTempCU );
      rpcTempCUPre = rpcTempCU;
      xCheckBestMode(rpcBestCU, rpcTempCU, depth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));
      if( bEnableTrans )
      {
        Bool bParentUseCSC = ( ((depth > 2)) && (m_ppcBestCU[depth -1]->getIntraBCCSCEnabled()) );
        if(!i && (!rpcTempCUPre->getQtRootCbf(0) || bParentUseCSC))
        {
          return;
        }
        else if(!i)
        {
          dCostFst = rdCost;
          if( rpcTempCU != rpcTempCUPre )
          {
            rpcTempCU->initEstData( depth, orgQP, bTransquantBypassFlag );
            rpcTempCU->copyPartFrom( rpcBestCU, 0, depth );
          }
        }
        else
        {
          if ( rpcTempCU == rpcTempCUPre )
          {
            rdCost = dCostFst;
          }
        }
      }
      else
      {
        return;
      }
    }
  }
  else
  {
    rdCost = MAX_DOUBLE;
  }
}

Void TEncCu::xCheckRDCostHashInter( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, Bool& isPerfectMatch DEBUG_STRING_FN_DECLARE(sDebug) )
{
  DEBUG_STRING_NEW(sTest)

  isPerfectMatch = false;
  UChar depth = rpcTempCU->getDepth( 0 );
  rpcTempCU->setDepthSubParts( depth, 0 );
  rpcTempCU->setSkipFlagSubParts( false, 0, depth );
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, depth );
  rpcTempCU->setPredModeSubParts( MODE_INTER, 0, depth );
  rpcTempCU->setMergeFlagSubParts( false, 0, 0, depth );

  if ( m_pcPredSearch->predInterHashSearch( rpcTempCU, m_ppcPredYuvTemp[depth], isPerfectMatch ) )
  {
    TComDataCU *rpcTempCUPre = NULL;
    Bool   bEnableTrans      = rpcBestCU->getSlice()->getPPS()->getPpsScreenExtension().getUseColourTrans();
    UChar  colourTransform = 0;
    Bool   bRGB              = m_pcEncCfg->getRGBFormatFlag();
    SChar   orgQP            = rpcTempCU->getQP( 0 );
    Bool   bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);
    const TComSPS &sps=*(rpcTempCU->getSlice()->getSPS());
    if ( bTransquantBypassFlag && (sps.getBitDepth( CHANNEL_TYPE_LUMA ) != sps.getBitDepth( CHANNEL_TYPE_CHROMA )) )
    {
      bEnableTrans = false;
    }

    TComYuv* pcTmpPredYuv = m_ppcPredYuvTemp[depth];
    for(UInt i = 0;  i < 2 ; i++)
    {
      colourTransform = ( bRGB && bEnableTrans )? (1-i): i;
      if ( i == 0 )
      {
        if ( bEnableTrans )
        {
          if ( !getEnableInterTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], (colourTransform? ACT_TRAN_CLR: ACT_ORG_CLR) DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TWO_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
        else
        {
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                    );
        }
      }
      else
      {
        if ( m_pcEncCfg->getRGBFormatFlag() )
        {
          if ( !getEnableInterTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TRAN_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
        else
        {
          if ( !getEnableInterTUACT() )
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_TRAN_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
          else
          {
            m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[depth], pcTmpPredYuv, m_ppcResiYuvTemp[depth], m_ppcResiYuvBest[depth], m_ppcRecoYuvTemp[depth], false, m_ppcNoCorrYuv[depth], ACT_ORG_CLR DEBUG_STRING_PASS_INTO(sTest)
#if PCC_RDO_EXT
                                                       , m_ppcOccupancyYuv[depth]
#endif
                                                      );
          }
        }
      }
      rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

      if(!i)
      {
        Double rdCost1 = rpcTempCU->getTotalCost();
        if( rdCost1 < m_ppcBestCU[depth]->getTmpInterRDCost() )
        {
          m_ppcBestCU[depth]->setTmpInterRDCost(rdCost1);
          m_ppcBestCU[depth]->setInterCSCEnabled(true);
          m_ppcTempCU[depth]->setTmpInterRDCost(rdCost1);
          m_ppcTempCU[depth]->setInterCSCEnabled(true);
        }
      }
      else
      {
        Double rdCost = rpcTempCU->getTotalCost();
        if( rdCost < m_ppcBestCU[depth]->getTmpInterRDCost() )
        {
          m_ppcBestCU[depth]->setTmpInterRDCost(rdCost);
          m_ppcBestCU[depth]->setInterCSCEnabled(false);
          m_ppcTempCU[depth]->setTmpInterRDCost(rdCost);
          m_ppcTempCU[depth]->setInterCSCEnabled(false);
        }
      }

      xCheckDQP( rpcTempCU );
      rpcTempCUPre = rpcTempCU;

      xCheckBestMode( rpcBestCU, rpcTempCU, depth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest) );
      if(bEnableTrans)
      {
        if( !rpcBestCU->isInter(0) || rpcBestCU->isSkipped(0)  )
        {
          return;
        }
        Bool bParentUseCSC = ( (depth != 0) && (m_ppcBestCU[depth -1]->getInterCSCEnabled()) );
        if(!i && (!rpcTempCUPre->getQtRootCbf(0) || bParentUseCSC))
        {
          return;
        }
        else if(!i && rpcTempCUPre->getQtRootCbf(0))
        {
          if( rpcTempCU != rpcTempCUPre )
          {
            rpcTempCU->initEstData( depth, orgQP, bTransquantBypassFlag);
            rpcTempCU->copyPartFrom( rpcBestCU, 0, depth );
          }
        }
      }
      else
      {
        return;
      }
    }
  }
}

UInt TEncCu::xCheckPaletteMode(TComDataCU *&rpcBestCU, TComDataCU *&rpcTempCU, Bool forcePalettePrediction, UInt iterNumber, UInt *paletteSize)
{
  UInt depth = rpcTempCU->getDepth( 0 );

  rpcTempCU->setSkipFlagSubParts( false, 0, depth );
  rpcTempCU->setIPCMFlag(0, false);
  rpcTempCU->setIPCMFlagSubParts (false, 0, rpcTempCU->getDepth(0));
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, depth );
  rpcTempCU->setPredModeSubParts( MODE_INTRA, 0, depth );
  rpcTempCU->setTrIdxSubParts ( 0, 0, depth );
  rpcTempCU->setPaletteModeFlagSubParts(true, 0, rpcTempCU->getDepth(0));
  rpcTempCU->setChromaQpAdjSubParts( rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, depth );

#if PCC_RDO_EXT
  UInt testedModes=m_pcPredSearch->paletteSearch(rpcTempCU, m_ppcOrigYuv[depth], m_ppcPredYuvTemp[depth], m_ppcResiYuvTemp[depth],
                                                 m_ppcRecoYuvTemp[depth], forcePalettePrediction, iterNumber, paletteSize, m_ppcOccupancyYuv[depth]);
#else
  UInt testedModes=m_pcPredSearch->paletteSearch(rpcTempCU, m_ppcOrigYuv[depth], m_ppcPredYuvTemp[depth], m_ppcResiYuvTemp[depth],
                                                 m_ppcRecoYuvTemp[depth], forcePalettePrediction, iterNumber, paletteSize);
#endif

  xCheckDQP( rpcTempCU );
  DEBUG_STRING_NEW(a)
  DEBUG_STRING_NEW(b)
  xCheckBestMode(rpcBestCU, rpcTempCU, depth DEBUG_STRING_PASS_INTO(a) DEBUG_STRING_PASS_INTO(b));

  return (testedModes);
}

//! \}


} // namespace pcc_hm
