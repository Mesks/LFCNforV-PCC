#ifndef __occupancyDciDudgeCPP__
#define __occupancyDciDudgeCPP__

extern unsigned char** occupancyData;
extern int             occupancyHeight;
extern int             occupancyWidth;
extern int             occupancyFrameNumber;
extern bool            QPr5;
extern int             OorGorA;

void occupancyDciInit( int width, int height, int frameNum, string name ) {
  QPr5            = false;
  occupancyWidth  = width / 4;
  occupancyHeight = height / 4;
  occupancyFrameNumber = frameNum / 2;
  if ( strstr( name.c_str(), "_GOF0_attribute" ) != NULL )
    OorGorA = 1;
  else if ( strstr( name.c_str(), "_GOF0_geometry" ) != NULL )
    OorGorA = 0;
  else if ( strstr( name.c_str(), "_GOF0_occupancy" ) != NULL )
    OorGorA = -1;

  int npos_GOF0_, npos_x;
  npos_GOF0_           = name.find( "_GOF0_" );
  char*  temp          = new char;
  string occupancyName = name.substr( 0, npos_GOF0_ ) + "_GOF0_occupancy_" + itoa( occupancyWidth, temp, 10 ) + "x" +
                         itoa( occupancyHeight, temp, 10 ) + "_8bit_p420.yuv";
  FILE* oriOccupancy = fopen( occupancyName.data(), "rb+" );
  if ( oriOccupancy == nullptr ) {
    QPr5 = true;
    occupancyWidth *= 2;
    occupancyHeight *= 2;
    occupancyName = name.substr( 0, npos_GOF0_ ) + "_GOF0_occupancy_" + itoa( occupancyWidth, temp, 10 ) + "x" +
                    itoa( occupancyHeight, temp, 10 ) + "_8bit_p420.yuv";
    oriOccupancy = fopen( occupancyName.data(), "rb+" );
  }
  occupancyData = new unsigned char*[occupancyFrameNumber];
  for ( int i = 0; i < occupancyFrameNumber; i++ )
    occupancyData[i] = (unsigned char*)malloc( sizeof( unsigned char ) * occupancyHeight * occupancyWidth * 3 / 2 );

  if ( oriOccupancy != nullptr )
    for ( int i = 0; i < occupancyFrameNumber; i++ )
      fread( occupancyData[i], sizeof( unsigned char ), occupancyWidth * occupancyHeight * 3 / 2, oriOccupancy );
  if ( oriOccupancy != nullptr ) fclose( oriOccupancy );
}

void occupancyDciDestroy() {
  for ( int i = 0; i < occupancyFrameNumber; i++ ) { delete[] occupancyData[i]; }
  delete[] occupancyData;
}

void checkData( int width, int height, int frameNum ) {
  for ( int i = 0; i < frameNum; i++ ) {
    for ( int j = 0; j < height; j++ ) {
      cout << "¡¾" << j << "¡¿: ";
      for ( int k = 0; k < width; k++ ) { cout << occupancyData[i][j * occupancyWidth + k] - 0 << " "; }
      cout << endl;
    }
    cout << "-========================-" << endl;
  }
}

#endif  // __occupancyDciDudge__