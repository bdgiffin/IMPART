#ifndef COLOR_MAP_H
#define COLOR_MAP_H

#include<vector>
#include<utility> // std::pair

class ColorMap {
public:

  // table of interpolated coordinate vs. color values
  std::vector<std::pair<float,std::array<float,3> > > colorTable;

  // find the color range in the table, and interpolate the RGB value
  void getColor(float x, float (&rgb)[3]) {
    if (x <= colorTable.front().first) {
      auto& firstColor = colorTable.front().second;
      rgb[0] = firstColor[0];
      rgb[1] = firstColor[1];
      rgb[2] = firstColor[2];
    } else if (x >= colorTable.back().first) {
      auto& lastColor = colorTable.back().second;
      rgb[0] = lastColor[0];
      rgb[1] = lastColor[1];
      rgb[2] = lastColor[2];
    } else {
      int i = 0;
      while(x > colorTable[i+1].first) i++;
      float x1     = colorTable[i].first;
      auto& color1 = colorTable[i].second;
      float x2     = colorTable[i+1].first;
      auto& color2 = colorTable[i+1].second;
      if (x2 > x1) {
	float t = (x - x1) / (x2 - x1);
	rgb[0] = color1[0] + t*(color2[0] - color1[0]);
	rgb[1] = color1[1] + t*(color2[1] - color1[1]);
	rgb[2] = color1[2] + t*(color2[2] - color1[2]);
      } else {
	rgb[0] = 0.5*(color1[0]+color2[0]);
	rgb[1] = 0.5*(color1[1]+color2[1]);
	rgb[2] = 0.5*(color1[2]+color2[2]);
      }
    }
  } // getColor()
  
}; // ColorMap

#endif // COLOR_MAP_H
