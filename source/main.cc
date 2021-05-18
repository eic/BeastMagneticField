
#include <stdio.h>
#include <stdlib.h>

#include <BeastMagneticField.h>

// --------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  if (argc == 4) { 
    auto bmf = new BeastMagneticField(argv[1]);

    if (bmf->ValidMapImported()) {
      // Turn linear interpolation on;
      bmf->UseInterpolation();
      double r = atof(argv[2]), z = atof(argv[3]), br, bz;
      bool ret = bmf->GetFieldValue(r, z, br, bz);
      if (!ret)
	printf("{r,z} coordinates out of the map region!\n");
      else 
	//printf("%8.4f %8.4f%s\n", br, bz, bmf->IsInsideTheBore(r, 0.0) ? "" : " # outside of the bore");
	printf("%8.4f %8.4f [T]\n", br, bz);//, bmf->IsInsideTheBore(r, 0.0) ? "" : " # outside of the bore");
    }
    else
      printf("File '%s' either does not exist or is not a valid field map!\n", argv[1]);
  } else {
    printf("\n   usage (cylindrical coordinates {r,z} in [cm]): \n\n     %s ASCII-field-map-file r z\n", argv[0]);
    printf("\n   output (radial and longitudinal field components in [T]): br bz\n\n");
  } //if
} // main()

// --------------------------------------------------------------------------------------
