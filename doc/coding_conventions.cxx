// NAMING CONVENTIONS

// Coding conventions: use underscores to separate words.
// Name of class, struct, enum, namespace or typedef will begin with an uppercase letter.
// Functions and variables have only lowercase letters.
struct Thermo;
class Thermo_moist;
typedef std::map<std::string, Field3d*> Field_map;

void exec_stats();
void init_prognostic_variable();

// Accessors and mutators contain the exact name of the variable:
bool calc_mean_profs = true;
bool get_calc_mean_profs();
void set_calc_mean_profs(bool);

// Functions that refer to the used name of quantity use the exact same name:
void calc_mask_ql();

// Inline functions that are used inside of kernels can be kept short for readability:
inline double i_4(const double, const double, const double, const double);
a[ijk] = u[ijk] * i_4( s[ijk-ii2], s[ijk-ii1], s[ijk], s[ijk+ii1] );

// GPU versions of existing functions and variables are written with a suffix _g:
void boundary_cyclic_g();


// INDENTATION AND WHITE SPACE
// Only use FOUR spaces for indentation. Add whitespace after if, for, while, catch.
if (name == "b")
    return true;
else
    return false;

// OTHER IMPORTANT CONVENTIONS
// Only initialize variables when you need them, thus in the innermost scope. Declare 
// everything const that will be never changed.

// Correct:
for (int i=grid->istart; i<grid->iend; ++i)
{
    const int ijk = i + j*jj + k*kk;
    u[ijk] = a*s[ijk];
}

// Wrong:
int i;
int ijk;

for (i=grid->istart; i<grid->iend; ++i)
{
    ijk = i + j*jj + k*kk;
    u[ijk] = a*s[ijk];
}

