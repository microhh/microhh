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
void set_calcMeanProfs(true);
void get_gridData();

// 2.3 Functions that distinguish multiple versions for different spatial orders use underscore suffix:
void interpolate_4th();
void calcFlux_2nd();


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
