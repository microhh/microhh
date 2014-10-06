// 1. All classes, structs, typedefs and enums go PascalCase:
class ThermoMoist;
struct StatsVar;
typedef std::map<std::string, Field3d *> FieldMap;
enum BoundaryType {DirichletType, NeumannType, FluxType};


// 2. All functions go camelCase, function names should be informative:
void execStats();
void initPrognosticVariable;

// 2.1 All accessors and mutators contain the name of the variable without changing the case:
void get_calcProfs();
void set_calcProfs(true);
void get_gridData();

// 2.2 Functions that distinguish multiple versions for different spatial orders use underscore suffix:
void interpolate_4th();
void calcFlux_2nd();


// 3. Variable names go the the wish of the user:
bool calcProfs = true;
int icells;

