#ifndef RADIATION_H
#define RADIATION_H

class Radiation {
public:

  void initialize() {
    heat = 0.0;
  } // initialize()
  
  void timeIntegrate(float dt) {
    
  } // timeIntegrate()

  float position[2]; // position of the radiation source
  float heat;        // rate of heat loss through radiation
}; // Radiation

#endif // RADIATION_H
