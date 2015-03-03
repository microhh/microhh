// 1.1 All classes, structs, typedefs and enums go PascalCase:
class ThermoMoist;
struct StatsVar;
typedef std::map<std::string, Field3d *> FieldMap;
enum BoundaryType {DirichletType, NeumannType, FluxType};


// 2.1 All functions go camelCase, function names should be informative:
void execStats();
void initPrognosticVariable;

// 2.2 All accessors and mutators contain the name of the variable without changing the case:
bool get_calcMeanProfs();
void set_calcMeanProfs(bool);
void get_gridData();

// 2.3 Functions that distinguish multiple versions for different spatial orders use underscore suffix:
void interpolate_4th();
void calcFlux_2nd();

// 2.3 Functions that refer to the used name of quantity can also apply underscores:
void calcMask_ql();
void calcMask_qlcore();

// 2.4 Inline functions that are used inside of kernels can be kept short for readability:
inline double interp4(const double, const double, const double, const double);
a[ijk] = u[ijk]*interp4(s[ijk-ii2], s[ijk-ii1], s[ijk], s[ijk+ii1]);

// 3.1 Variable names go at the wish of the user:
bool calcProfs = true;
int icells;


// 4.1 GPU versions of existing functions and variables are written with a suffix _g:
void boundaryCyclic_g();

// 4.2 GPU kernels are stored in namespaces with the name of the class and a suffix _g:
namespace ThermoMoist_g
{
  calcBuoyancy();
}
