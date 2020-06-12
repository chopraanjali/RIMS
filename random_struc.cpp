/* Program to initialize system structure using random placement approach for molecules

  - the molecules are placed randomly and with a given minimum distance criteria

    *NOTE : larger molecules are generated first with much larger minimum distance (30% of given distance - hardcoded for now)
            then creates smaller molecule with actual minimum distance requirement

  - writes output file in xyz format

    * requires params.txt file in the current folder
    *  NOTE: params.txt file contains all information
    * requires single molecule File(s) in xyz format in the current folder
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <chrono>

#include "Instrumentor.h"

using namespace std;

#define PROFILING 1
#if PROFILING
#define PROFILE_SCOPE(name) InstrumentationTimer timer##__LINE__(name)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__FUNCTION__)
#else
#define PROFILE_SCOPE(name)
#endif


//GLOBAL PARAMETERS
ofstream outFile;                                 	//output File
double const pi = 3.14159265358979323;


//FUNCTIONS
struct params
{
  double                density;
  double                percent_solute;
  double                box_dim;
  double                r_box_dim;
  double                min_dist;
  string                out_file;
  string                solute_filename;
  string                solvent_filename;
  unsigned long int     seed;
  int                   tot_molec;
  int                   dimensions;
  int                   inp_mols;
  int                   solute_mols;
  int                   solvent_mols;
  int                   total_atoms_in_system;
 };

 struct Molecule
{
  int     num_of_atoms;       // total # of atoms in file
  int     num_of_molecules;   // total molecules to create for each type
  int     total_atoms;        // molecules*atoms in each molecule
  string  *elems;             // elems - array for elements from input files
  string  file_name;
  double  **xyz_co;           // array for storing atom coordinates from input files
  double  **com_co;           // array for generated com coordinates
  double  **com_angles;       // array for generated com angles
  double  *barycenter;        // array for storing barycenter of each molecule
} ;


struct coordinates
{
  double x, y, z;
  string elems;
};

// vector of struct to store final coordinates
vector<coordinates> coords;

//adding a vector element to avoid accessing an empty vector error
void add_coords(vector<coordinates> &Coords)
{
  coordinates coord;
  coord.x = 0;
  coord.y = 0;
  coord.z = 0;
  coord.elems = "0";

  Coords.push_back(coord);
}

//------------------BEGINNING OF FUNCTION DECLARATIONS--------------------------------


//-------------------------------------------------------------------------------------
//params.txt file Parser class
//-------------------------------------------------------------------------------------

class Parameters
{
public:
    // clear all values
    void Clear();

    // load a parameter file
    bool Load(const string& File);

    // check if value associated with given key exists
    bool Contains(const string& key) const;

    // get value associated with given key
    bool Get(const string& key, string& value) const;
    bool Get(const string& key, int&    value) const;
    bool Get(const string& key, long unsigned int&   value) const;
    bool Get(const string& key, double& value) const;
    bool Get(const string& key, bool&   value) const;

private:
    // the container
    map<string,string> data;

    // remove leading and trailing tabs and spaces
    static string Trim(const string& str);
};

void Parameters::Clear()
{
    data.clear();
}

bool Parameters::Load(const string& file)
{
    ifstream inFile(file.c_str());

    if (!inFile.good())
    {
        cout << "Cannot read configuration file " << file << endl;
        return false;
    }

    while (inFile.good() && ! inFile.eof())
    {
        string line;
        getline(inFile, line);

        // filter out comments
        if (!line.empty())
        {
            int pos = line.find('#');

            if (pos != string::npos)
            {
                line = line.substr(0, pos);
            }
        }

        // split line into key and value
        if (!line.empty())
        {
            int pos = line.find('=');

            if (pos != string::npos)
            {
                string key     = Trim(line.substr(0, pos));
               //cout << key << endl;
                string value   = Trim(line.substr(pos + 1));
                //cout << value << endl;

                if (!key.empty() && !value.empty())
                {
                    data[key] = value;
                }
            }
        }
    }

    return true;
}

bool Parameters::Contains(const string& key) const
{
    return data.find(key) != data.end();
}

bool Parameters::Get(const string& key, string& value) const
{
    map<string,string>::const_iterator iter = data.find(key);

    if (iter != data.end())
    {
        value = iter->second;
        return true;
    }
    else
    {
        return false;
    }
}

bool Parameters::Get(const string& key, int& value) const
{
    string str;

    if (Get(key, str))
    {
        value = atoi(str.c_str());
        return true;
    }
    else
    {
        return false;
    }
}

bool Parameters::Get(const string& key, long unsigned int& value) const
{
    string str;

    if (Get(key, str))
    {
        value = atol(str.c_str());
        return true;
    }
    else
    {
        return false;
    }
}

bool Parameters::Get(const string& key, double& value) const
{
    string str;

    if (Get(key, str))
    {
        value = atof(str.c_str());
        return true;
    }
    else
    {
        return false;
    }
}

bool Parameters::Get(const string& key, bool& value) const
{
    string str;

    if (Get(key, str))
    {
        value = (str == "true");
        return true;
    }
    else
    {
        return false;
    }
}

string Parameters::Trim(const string& str)
{
    int first = str.find_first_not_of(" \t");

    if (first != string::npos)
    {
        int last = str.find_last_not_of(" \t");

        return str.substr(first, last - first + 1);
    }
    else
    {
        return "";
    }
}
//-----------------------END OF FILE PARSER STUFF---------------------------------------


//--------------------read parameters from params.txt file------------------------------
params readParamsFile()
{

  params p;
  Parameters par;

  par.Load("params.txt");

  if (par.Get("density", p.density) &&
      par.Get("percent_solute", p.percent_solute) &&
      par.Get("box_dim", p.box_dim) &&
      par.Get("min_dist", p.min_dist) &&
      par.Get("out_file", p.out_file) &&
      par.Get("seed", p.seed) &&
      par.Get("dimensions", p.dimensions) &&
      par.Get("inp_mols", p.inp_mols) &&
      par.Get("solute_filename", p.solute_filename) &&
      par.Get("solvent_filename", p.solvent_filename))
  {
    cout << "Successfully read all parameters from params.txt ! " << endl;
  }
  else
  {
    cout << "Missing parameter in configuration file. Please check." << endl;
  }

  p.tot_molec = (int)((p.density * p.box_dim * p.box_dim * p.box_dim) + 0.5); //cast to int for proper rounding and conversion
  p.r_box_dim = 1.0 / p.box_dim;

  p.solute_mols = (int)((p.tot_molec * p.percent_solute / 100) + 0.5);
  if (p.inp_mols == 1)
  {
    p.solute_mols = 0;
  }
  p.solvent_mols = p.tot_molec - p.solute_mols;

  if (p.tot_molec < 1)
  {
    cerr << "Error : Number density too low!" << endl;
    exit(1);
  }


  cout << "Chosen parameters for creating initial structure | "
       << "Number density : " << p.density << " | "
       << "Box side : " << p.box_dim << " | "
       << "Seed : " << p.seed << endl;
  cout << "Creating " << p.solute_mols << " molecules of solute and " << p.solvent_mols << " molecules of solvent from a total of " << p.tot_molec << endl;

  return p;
}

//---------------------ALLOCATE MEMORY--------------------------------------------------
void allocateMemory(struct Molecule *mol, struct params &s)
{

  for (int i = 0; i < s.inp_mols; i++)
  {

    mol[i].elems = new string[mol[i].num_of_atoms];
    mol[i].xyz_co = new(nothrow) double*[mol[i].num_of_atoms];
                    for (int j = 0; j < mol[i].num_of_atoms; j++)
                    {
                    mol[i].xyz_co[j] = new double[s.dimensions];
                    }
    mol[i].com_co = new(nothrow) double*[mol[i].num_of_molecules];
    mol[i].com_angles = new(nothrow) double*[mol[i].num_of_molecules];
                    for (int j = 0; j < mol[i].num_of_molecules; j++)
                    {
                    mol[i].com_co[j] = new double[s.dimensions];
                    mol[i].com_angles[j] = new double[s.dimensions];
                    }
    mol[i].barycenter = new(nothrow) double[s.dimensions];
  }

	//cout << "Memory allocated!" << endl;
}

//-----------------CHECKS AND OPENS FILES IN CURRENT DIRECTORY (IMP FOR MEM ALLOCATION)---------------------------
void checkAndOpenFiles(params &s, struct Molecule *mol)
{
  PROFILE_FUNCTION();
  int mols = 0;
  string line, a;
  s.total_atoms_in_system = 0;

  ifstream file;
  for (int i = 0; i < s.inp_mols ; i++)
  {
        if (i == 0) { a = s.solvent_filename;}
        else { a = s.solute_filename;}
        file.open(a.c_str());

        if (!file) //unable to open file error
        {
          cerr << "Unable to open file " << a << endl;
          cerr << "Check if it exists in the current folder!" << endl;
          exit(1); // call system to stop
        }

        getline(file, line); //gets num of atoms
        mol[i].num_of_atoms = abs(atoi(line.c_str()));

        if (i == 0) { mol[i].num_of_molecules = s.solvent_mols;}
        else { mol[i].num_of_molecules = s.solute_mols;}

        mol[i].total_atoms = mol[i].num_of_atoms * mol[i].num_of_molecules;

        //cerr << "Number of molecules " << mol[i].num_of_molecules << endl;
        s.total_atoms_in_system += mol[i].total_atoms;
        //cout << "Molecule " << i+1 << " has atoms " << mol[i].total_atoms << endl;
        mol[i].file_name = a;
        //(s.file_name).push_back(a);
        file.close();
  }
}


//---------------------READS XYZ DATA FROM INPUT FILES AND STORES IT----------------------------------------
void readMoleculeXyzFile(struct params &s, struct Molecule *mol)
{
   PROFILE_FUNCTION();
    //Reading single molecule file
    istringstream iss;
    string line;

    for (int j = 0; j < s.inp_mols; j++)
    {
      //cerr << "Filename : " << s.file_name[j] << endl;
      //ifstream inFile(s.file_name[j]);

      ifstream inFile(mol[j].file_name);
      getline(inFile, line);                              //gets num_of_atoms
      getline(inFile,line);                               //gets comment line
      for (int i = 0; i < mol[j].num_of_atoms; i++) 	    //assigns coordinates from file to array
      {
        iss.clear();
        getline(inFile,line);

      if(line.length() == 0)        	//empty line error
      {
        cerr << "Error: empty line " << endl;
        exit(3);
      }

      iss.str(line);
      iss >> mol[j].elems[i];
      iss >> mol[j].xyz_co[i][0];
      iss >> mol[j].xyz_co[i][1];
      iss >> mol[j].xyz_co[i][2];

        //cout << mol[j].elems[i] << "|" << mol[j].xyz_co[i][0] << "|" << mol[j].xyz_co[i][1] << "|" << mol[j].xyz_co[i][2] << endl;
      }
      inFile.close();
    }
    //cout << "Done reading files" << endl;
}

//---------------------MOLECULE SORTING FOR CREATING MORE SPREAD OUT STRUCTURE ---------------------------------------
bool compareMolecules(Molecule lhs, Molecule rhs)
{
  return lhs.num_of_atoms > rhs.num_of_atoms;
}

void sortMolecules(struct Molecule *mol, struct params &s)
{
  PROFILE_FUNCTION();
  sort(mol, mol+s.inp_mols, compareMolecules);
  for (int i = 0 ; i < s.inp_mols; i++)
  {
    //cerr << "Molecules sorted!" << endl;
    //cerr << mol[i].num_of_atoms << endl;
  }
}



//---------------------DISTANCE CHECKING BETWEEN EVERY POSSIBLE PAIR OF ATOMS ---------------------------------------
double distChecker(double xyz_temp[][3], vector<coordinates> coords, struct params &s, int j, int k)
{
  PROFILE_FUNCTION();
	double DIST;
	double dist_comp[3];

    //cerr << xyz_temp[j][dim] << endl;
		dist_comp[0] = coords[k].x - xyz_temp[j][0];
		dist_comp[1] = coords[k].y - xyz_temp[j][1];
		dist_comp[2] = coords[k].z - xyz_temp[j][2];

	DIST = sqrt ((dist_comp[0]*dist_comp[0])  + (dist_comp[1]*dist_comp[1]) + (dist_comp[2]*dist_comp[2]));

	DIST -= static_cast<int>(DIST*s.r_box_dim + 0.5)*s.box_dim ;		//effective distance with PBC
	DIST = fabs(DIST);											                        //absolute value of distance
//	cerr << "Distance is : " << DIST << endl;
	return DIST;
}

//--------------------IMPLEMENTING PERIODIC BOUNDARY CONDITIONS ---------------------------------------
double makePeriodic(double a, double length)
{
  PROFILE_FUNCTION();
	if (a < 0.0)
  {
	  a = a + length;
  }
  else if (a > length)
  {
	  a = a - length;
  }
  else
  {
	 a = a;
	}
	return a;
}

//---------------FIND BARYCENTER FOR EACH TYPE OF GIVEN MOLECULE FOR APPROPRIATE ROTATION------------------------------
void findBarycenter(struct Molecule *mol, struct params &s)
{
  PROFILE_FUNCTION();
  double min[s.inp_mols][s.dimensions]{};
  double max[s.inp_mols][s.dimensions]{};

  for (int i = 0; i < s.inp_mols; i++)
  {
    double x[mol[i].num_of_atoms]{};
    double y[mol[i].num_of_atoms]{};
    double z[mol[i].num_of_atoms]{};
    for (int j = 0; j < mol[i].num_of_atoms; j++)
    {
      x[j] = mol[i].xyz_co[j][0];
      y[j] = mol[i].xyz_co[j][1];
      z[j] = mol[i].xyz_co[j][2];


    min[i][0] = *min_element(x, x+mol[i].num_of_atoms);
    min[i][1] = *min_element(y, y+mol[i].num_of_atoms);
    min[i][2] = *min_element(z, z+mol[i].num_of_atoms);

    max[i][0] = *max_element(x, x+mol[i].num_of_atoms);
    max[i][1] = *max_element(y, y+mol[i].num_of_atoms);
    max[i][2] = *max_element(z, z+mol[i].num_of_atoms);

    }

    for (int dim = 0; dim < s.dimensions; dim++)
    {
      //cerr << "min " << min[i][dim] << " & max " << max[i][dim] << endl;
      mol[i].barycenter[dim] = (max[i][dim] + min[i][dim])/2.0;
      //cerr << "barycenter is for mol " << i << " is " << mol[i].barycenter[dim] <<endl;
    }
  }

}

//---------------GENERATE RANDOM COORDINATES AND CREATE INITIAL STRUCTURE------------------------------
void generateSystem(struct Molecule *mol, struct params &s)
{
  PROFILE_FUNCTION();

  //RNG generator : generates uniformly distributed floating point number in given range
  //  http://www.cplusplus.com/reference/random/uniform_real_distribution/

  bool generation_successful;
  //cout << "Initializing random number generator with seed " << s.seed << endl;
  mt19937 generator(s.seed);
  uniform_real_distribution<double> distribution(0.0, 1.0);
  //cout << "PRNG initiated!" << endl;

  findBarycenter(mol,s);
  add_coords(coords);
  int curr_atoms = 1;
  int n = s.inp_mols;
	for (int p = 0; p < s.inp_mols ; p++)
	{
		cout << "Creating " << mol[p].total_atoms << " total atoms of " << mol[p].num_of_molecules << " molecules of " << mol[p].file_name << endl;
		double min_dist = s.min_dist + (n * 0.3) ;        //creates bigger molecule with much larger min dist requirement
    //double min_dist_sqr = (s.min_dist_sqr * s.min_dist_sqr) + (n * n * 0.3 * 0.3) + ( 2.0 * s.min_dist_sqr * n * 0.3 );
		for (int j = 0; j < mol[p].num_of_atoms; j++)     //translates molecule to origin based on user's chosen com_id
		{
			mol[p].xyz_co[j][0] = mol[p].xyz_co[j][0] - mol[p].barycenter[0];
			mol[p].xyz_co[j][1] = mol[p].xyz_co[j][1] - mol[p].barycenter[1];
			mol[p].xyz_co[j][2] = mol[p].xyz_co[j][2] - mol[p].barycenter[2];
		}

		int curr_total_molec = 0;
		for (int i = 0 ; i < mol[p].num_of_molecules ; i++)
		{
      double xyz_temp[mol[p].num_of_atoms][3]{};
      coordinates c;

			do
			{

				for (int dim = 0; dim < s.dimensions ; dim++)
				{
					// asssign random coordinates and angles to molecule
					mol[p].com_co[i][dim] = distribution(generator) * s.box_dim ;
					mol[p].com_angles[i][dim] = distribution(generator)  * 2 * pi ;
				}

        //cerr << "COM coords and angles created for molecule # " << i << " of type " << p << endl;

				double sin_a[mol[p].num_of_molecules];  double sin_b[mol[p].num_of_molecules];  double sin_g[mol[p].num_of_molecules];
				double cos_a[mol[p].num_of_molecules];  double cos_b[mol[p].num_of_molecules];  double cos_g[mol[p].num_of_molecules];

				sin_a[i] = sin(mol[p].com_angles[i][0]);  sin_b[i] = sin(mol[p].com_angles[i][1]);  sin_g[i] = sin(mol[p].com_angles[i][2]);
				cos_a[i] = cos(mol[p].com_angles[i][0]);  cos_b[i] = cos(mol[p].com_angles[i][1]);  cos_g[i] = cos(mol[p].com_angles[i][2]);

        //cerr << "Sines and cosines calculated for molecule # " << i << " of type " << p << endl;

				for (int j = 0; j < mol[p].num_of_atoms; j++)
				{

					double coords_atoms_temp[mol[p].num_of_atoms][s.dimensions]{};
					// rotates molecules about origin using 3d rotation matrix R = Rx(g)*Ry(b)*Rz(a) ; u' = R*u

					coords_atoms_temp[j][0] = (mol[p].xyz_co[j][0] ) *  cos_a[i] * cos_b[i]  +
											              (mol[p].xyz_co[j][1] ) * (cos_a[i] * sin_b[i]  * sin_g[i] - sin_a[i] * cos_g[i]) +
																    (mol[p].xyz_co[j][2] ) * (cos_a[i] * sin_b[i]  * cos_g[i] + sin_a[i] * sin_g[i]) ;


					coords_atoms_temp[j][1] = (mol[p].xyz_co[j][0] ) *  sin_a[i] * cos_b[i]  +
					                        	(mol[p].xyz_co[j][1] ) * (sin_a[i] * sin_b[i]  * sin_g[i] + cos_a[i] * cos_g[i]) +
					                        	(mol[p].xyz_co[j][2] ) * (sin_a[i] * sin_b[i]  * cos_g[i] - cos_a[i] * sin_g[i])  ;


					coords_atoms_temp[j][2] = (mol[p].xyz_co[j][0] ) * (-1)      * sin_b[i]  +
					                        	(mol[p].xyz_co[j][1] ) * cos_b[i]  * sin_g[i]  +
					                        	(mol[p].xyz_co[j][2] ) * cos_b[i]  * cos_g[i]  ;

					for (int dim = 0; dim < s.dimensions; dim++)
					{
						xyz_temp[j][dim] =	makePeriodic( mol[p].com_co[i][dim] + coords_atoms_temp[j][dim], s.box_dim );
					}

				  generation_successful = true;

          //cerr << "molecule # " << i << " of type " << p  << " generated,rotated and placed at com "<< endl;
					//check distance between all atoms in temp and all atoms created
          //double DIST = min_dist;

          for (int k = 1; k <= coords.size() ; ++k)
          {
           	double d = distChecker(xyz_temp,coords, s, j, k);
            //cerr << "dist sqr is " << d_sqr << endl;
            if ( d < min_dist)
             {
                //DIST = d;
                generation_successful = false;
                break;
                //cerr << "Current minimum distance is " << DIST << endl;
             }
          }

          if(generation_successful == false)
          {
            break;
          }

/*           if (DIST < min_dist)
          {
            generation_successful = false;
            //cerr << "current minimum distance is lower than given criteria, generating new coordinates " << endl;
          	break;
          } */
				}
        //cerr <<  i << " molecules created with given criteria of type " << p << endl;
			}
			while (generation_successful == false);
      // now assign temp to main array
      curr_atoms = curr_atoms + mol[p].num_of_atoms ;
      cerr << curr_total_molec+1 << " molecules created! " << endl;
      //cerr << "Current atoms are " << curr_atoms << endl;
      for (int Atoms = 0; Atoms < mol[p].num_of_atoms; Atoms++)
      {
        c.x = xyz_temp[Atoms][0];
        c.y = xyz_temp[Atoms][1];
        c.z = xyz_temp[Atoms][2];
        c.elems = mol[p].elems[Atoms];
        coords.push_back(c);
      }
			curr_total_molec++;
		}
    n--;
	}
}

//-----------------------WRITE .XYZ FILE WITH FINAL COORDINATES---------------------
void writeOutputXyzFile(struct Molecule *mol, struct params &s, vector<coordinates> &coords)
{
  PROFILE_FUNCTION();
  outFile.open(s.out_file.c_str());
  outFile << s.total_atoms_in_system << endl;
  outFile << "!comment line" << endl;
    for (int i = 1; i <= s.total_atoms_in_system ; i++)
    {
      outFile << coords[i].elems << "\t" << coords[i].x << "\t" << coords[i].y << "\t" << coords[i].z << endl;
    }
  cout << "Output file written!" << endl;
  outFile.close();
}

//----------------------DEALLOCATE MEMORYYYYYYYYYYY-------------------------------
void deallocateMemory(struct Molecule *mol, struct params &s)
{
  PROFILE_FUNCTION();
  for (int i = 0; i < s.inp_mols; i++)
  {
    delete[] mol[i].elems;

    for (int j = 0; j < mol[i].num_of_atoms; j++) //array for rotated atoms (FINAL)
    {
      delete[] mol[i].xyz_co[j];
    }
    delete[] mol[i].xyz_co;

    for (int j = 0; j < mol[i].num_of_molecules; j++)
    {
      delete[] mol[i].com_co[j];
      delete[] mol[i].com_angles[j];
    }
    delete[] mol[i].com_co;
    delete[] mol[i].com_angles;
    delete[] mol[i].barycenter;
  }

  //cout << "Memory deallocated!" << endl;
}


//-----------------------------------------------------------------------
//------------------END OF ALL FUNCTION DECLARATIONS---------------------
//-----------------------------------------------------------------------


int main()
{
  Instrumentor::Get().BeginSession("Session Name");         // Begin session
  {
  params s;
  s = readParamsFile();
  Molecule mol[s.inp_mols];


  checkAndOpenFiles(s,mol);

  allocateMemory(mol, s);

  readMoleculeXyzFile(s, mol);


  if (s.inp_mols > 1)
  {
  sortMolecules(mol, s); //sorting by number of atoms in each file
  }
  //cout << "Creating " << s.tot_molec << " molecules in a simulation box of side " << s.box_dim << " units!" << endl;
  generateSystem(mol, s);
  writeOutputXyzFile(mol, s, coords);

  cout << "Initial structure created with random coordinates and random orientation!" << endl;

  deallocateMemory(mol, s);
  }

  Instrumentor::Get().EndSession();                        // End Session
  return 0;
}
//# Initial-Structure
