
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// NB: shamelessly use the fact that as of 6.18 (?) inheritance from 
// TObject is not required for scripting if a class instance does not 
// need to be streamed out; therefore no changes to the class definitions
// in BeastMagneticField.h;
#pragma link C++ class BeastMagneticField+;
#pragma link C++ class BeastMagneticFieldCell+;

#endif
