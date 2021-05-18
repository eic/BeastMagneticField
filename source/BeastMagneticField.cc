
#include <set>

#include <stdio.h>
#include <string.h>
#include <assert.h>

// Well, one can of course determine this from the input data, but what's the point?;
#define _CELL_SIZE_               (2.0)

// See the original solenoid.C script; prefer to indicate a cutoff;
//#define _CRYOSTAT_INNER_RADIUS_ (137.0)

#include <BeastMagneticField.h>

// For the time being keep this code, which imports the original 
// 5.7MB Elmer file with the extra columns and odd formatting and dumps the simplified ASCII file;
//#define _USE_ORIGINAL_ELMER_FILE_

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------

BeastMagneticField::BeastMagneticField(const char *fname, bool enforce_axial_symmetry): 
  mCellSize(_CELL_SIZE_), mScale(1.0), mRdim(0), mZdim(0), mRmin(0.0), mZmin(0.0), 
  mUseInterpolation(false), mFieldMap(0)
{
  FILE *fmap = fopen(fname, "r");
  if (!fmap) {
    printf("Failed to open '%s' for reading!\n", fname);
    return;
  } //if
#ifdef _USE_ORIGINAL_ELMER_FILE_
  FILE *fout = fopen("mfield.4col.dat", "w");
#endif

  {
    // I kind of know that the strings in this file are short, right?;
    char buffer[1024];
    std::vector<BeastMagneticFieldCell> cells;
    std::set<double> rset, zset;

    // Prefer to import everything, and only after this start handling data; 
    while (fgets(buffer, 1024-1, fmap)) {
      double r, z, br, bz;
      // C-style, sorry;
#ifdef _USE_ORIGINAL_ELMER_FILE_
      double dummy;
      sscanf(buffer, "%le %lf %lf %le %le %le %le", &dummy, &dummy, &dummy, &r, &z, &br, &bz);

      // Convert [m] to [cm]; store the values;
      r *= 100.0; z *= 100.0;
      fprintf(fout, "%6.1f %6.1f %7.4f %7.4f\n", r, z, br, bz);
#else
      sscanf(buffer, "%lf %lf %lf %lf", &r, &z, &br, &bz);
#endif      

      // Elmer gives a (very) small radial component on the symmetry axis; fix this;
      if (!r && enforce_axial_symmetry) br = 0.0;
      
      rset.insert(r); zset.insert(z);

      cells.push_back(BeastMagneticFieldCell(r, z, br, bz));
    } // while

    mRdim = (int)rint((*rset.rbegin() - *rset.begin())/mCellSize) + 1;
    mZdim = (int)rint((*zset.rbegin() - *zset.begin())/mCellSize) + 1;
    mRmin = *rset.begin(); mZmin = *zset.begin();

    // Sanity check; 
    if (mRdim*mZdim != cells.size()) {
      // And mFieldMap[] in this case will remain NULL;
      printf("Line count in the input file (%d) does not match the grid size (%d x %d)\n",
	     (int)cells.size(), mRdim, mZdim);
      return;
    } //if

    // Populate the internal array;
    {
      unsigned dim = mRdim*mZdim;
      bool ok[dim]; memset(ok, false, sizeof(ok));
      mFieldMap = new std::pair<double, double>[dim];

      for(auto cell: cells) {
	// Want the closest value here -> use rint();
	int ir = (int)rint((cell.mR - mRmin)/mCellSize);
	int iz = (int)rint((cell.mZ - mZmin)/mCellSize);
	// FIXME: do it better; should not happen anyway;
	assert(ir >=0 && ir < (int)mRdim && iz >= 0 && iz < (int)mZdim);

	// Calculate 1D index;
	unsigned id = ir*mZdim + iz;
	// Check for double-counting; if never happens, all the cells are filled; good;
	if (ok[id]) {
	  printf("Duplicate cells in the map file (%d,%d)!\n", ir, iz);
	  return;
	} //if
	ok[id] = true; mFieldMap[id] = std::make_pair(cell.mBR, cell.mBZ);
      } //for cell
    } 
  } 

  fclose(fmap);
#ifdef _USE_ORIGINAL_ELMER_FILE_
  fclose(fout);
#endif
} // BeastMagneticField::BeastMagneticField()

// --------------------------------------------------------------------------------------

// Ignore the z-coordinate I guess?;
//bool BeastMagneticField::IsInsideTheBore(double x, double y/*, double z*/) const
//{
//return (sqrt(x*x+y*y) < _CRYOSTAT_INNER_RADIUS_);
//} // BeastMagneticField::IsInsideTheBore()

// --------------------------------------------------------------------------------------

bool BeastMagneticField::GetFieldValue(double   x, double   y, double   z, 
				       double &bx, double &by, double &bz) const
{
  double r = sqrt(x*x + y*y), br;

  bool ret = GetFieldValue(r, z, br, bz);
  if (!ret) return false;

  double fi = atan2(y, x);
  bx = br*cos(fi); by = br*sin(fi);

  return true;
} // BeastMagneticField::GetFieldValue()

// --------------------------------------------------------------------------------------

bool BeastMagneticField::GetFieldValue(double r, double z, double &br, double &bz) const
{
  // Reset per default, in one place;
  br = bz = 0.0;

  // Sanity check;
  if (r < 0.0) return false;

  // Well, this is not exactly clean, but kind of makes sense;
  if (mUseInterpolation) {
    // Want a 2x2 cell group here -> use floor();
    int ir = (int)floor((r - mRmin)/mCellSize);
    int iz = (int)floor((z - mZmin)/mCellSize);

    // Range check; NB: watch '-1' in this case;
    if (ir < 0 || ir >= (int)mRdim-1 || iz < 0 || iz >= (int)mZdim-1) return false;

    // Arrange by hand something like a 2D linear interpolation;
    double qdr = (r - mRmin - ir*mCellSize)/mCellSize, qdz = (z - mZmin - iz*mCellSize)/mCellSize; 
    double wt[2][2] = {{(1.0-qdr)*(1.0-qdz), (1.0-qdr)*qdz}, 
		       {     qdr *(1.0-qdz),      qdr *qdz}};
    double wtsum = 0.0;
    for(unsigned ip=0; ip<2; ip++)
      for(unsigned iq=0; iq<2; iq++)
	wtsum += wt[ip][iq];
    for(unsigned ip=0; ip<2; ip++)
      for(unsigned iq=0; iq<2; iq++)
	wt[ip][iq] /= wtsum;
    
    for(unsigned ip=0; ip<2; ip++)
      for(unsigned iq=0; iq<2; iq++) {
	const std::pair<double, double> &cell = mFieldMap[(ir+ip)*mZdim + iz+iq];

	br += wt[ip][iq]*cell.first;
	bz += wt[ip][iq]*cell.second;
      } //for ip..iq
  } else {
    // Want the closest value here -> use rint();
    int ir = (int)  rint((r - mRmin)/mCellSize);
    int iz = (int)  rint((z - mZmin)/mCellSize);

    // Range check;
    if (ir < 0 || ir >= (int)mRdim   || iz < 0 || iz >= (int)mZdim)   return false;
    
    unsigned id = ir*mZdim + iz;
    const std::pair<double, double> &cell = mFieldMap[id];
    br = cell.first; bz = cell.second;
  } //if 

  // Rescale in a single place;
  br *= mScale; bz *= mScale;

  return true;
} // BeastMagneticField::GetFieldValue()

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
