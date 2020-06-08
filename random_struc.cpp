/* Program to initialize system structure using random placement approach for molecules

  - the molecules are placed randomly and with a minimum distance of 6 units
  - writes output file in xyz format

    * requires single molecule File in xyz format
    * Parameters to be given:
      - desired number density
      - extents of simulation box
      - input filename and output filename
*/

/*
TODO :
 - Write help for % of solute
 - Create function to read input files from folder DONE
 - If statement for only solvent molecules
 - Set limit for loops and possibly restart with another seed
 - Generate system for larger molecules first with much larger minimum distance (10% of given distance)
 - Generate system for smaller molecules with given minimum distance requirement

*/
//  cd /mnt/c/Users/Anjali/Dropbox/Thesis/intial_structure
//	g++ random_struc.cpp -o random_struc.exe
//	./random_struc.exe 0.0002 20 5 24 test1.xyz prop_gly.xyz CNT.xyz

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <algorithm>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <map>

using namespace std;

//GLOBAL PARAMETERS
ofstream outFile;                                 	//output File
double const pi = 3.14159265358979323;

//FUNCTIONS
struct params
{
  float                 density;
  float                 percent_solute;
  float                 box_dim;
  double                r_box_dim;
  double                min_dist;
  string                out_file;
  unsigned long int     seed;
  int                   tot_molec;
  int                   dimensions;
  int                   inp_mols;
  int                   solute_mols;
  int                   solvent_mols;
  int                   total_atoms_in_system;
  vector<string>        file_name;
 };

 struct Molecule
{
  int     num_of_atoms;   // total # of atoms in file
  int     num_of_molecules;    //total molecules to create for each type
  int     total_atoms;    // molecules*atoms in each molecule
  string  *elems;        // elems - array for elements
  double  **xyz_co;      //array for storing atom coordinates
  double  **com_co;      //array for generated com coordinates
  double  **com_angles; //array for generated com angles
  double  *barycenter;
} ;

struct coordinates
{
  double x, y, z;
  string elems;
};

vector<coordinates> coords;

void add_coords(vector<coordinates> &Coords)
{
  coordinates coord;
  coord.x = 0;
  coord.y = 0;
  coord.z = 0;
  coord.elems = "0";

  Coords.push_back(coord);
}

void printHelp()
{
    cout << "Not enough arguments!" << endl;
    cout << "Try with: " << endl;
    cout << "./random_struct.exe <float desired_num_density> <float box_length> <int seed> <string output_filename.xyz> <int num_of_input_molecules> <float percentage_solute>" << endl;
    cout << "For example: " << endl;
    cout << "./random_struc.exe 0.0002 20 5 test1.xyz 35 2" << endl;
    cout << "*********** HELP ****************" << endl;
    cout << "   - Number density has units molecules/Angstrom^3. " << endl;
    cout << "   - Box length has units Angstrom." << endl;
    cout << "   - Provide a seed for random number generator (0 to 4,294,967,295)" << endl;
    cout << "   - Name of output file with .xyz extension." << endl;
    cout << "   - Percentage of solute in mixture. " << endl;
    cout << "   - To create a pure system, num_of_input_molecules = 1 and percentage_solute = 0" << endl;

}

params getAndPrintRunParameters(int argc, char *argv[])
{
  params p;
  float  f;

  p.density         =   atof(argv[1]);   // number density desired for liquid structure
  p.box_dim         =   atof(argv[2]);   //extents of the sim box
  p.seed            =   atoi(argv[3]);
  p.out_file        =   argv[4];
  p.inp_mols        =   atoi(argv[5]);
  p.percent_solute  =   atof(argv[6]);
  p.min_dist        =   3.0;             //twice the largest distance from com (approx)
  p.dimensions      =   3;


  f             =  p.density * p.box_dim*p.box_dim*p.box_dim;
  p.tot_molec   = (int)(f+0.5);    //cast to int for proper rounding and conversion
  p.r_box_dim   = 1.0/p.box_dim;


  p.solute_mols = (int)((p.tot_molec * p.percent_solute / 100) + 0.5);
  if (p.inp_mols == 1)
  {
    p.solute_mols = 0;
  }
  p.solvent_mols = p.tot_molec - p.solute_mols;

	if (argc < 6)  //parameter check
  {
    printHelp();
    exit(1);
  }

  if (p.tot_molec < 1)
  {
  cerr << "Error : Number density too low!" << endl ;
  exit(1);
  }

  cout << "Chosen parameters for creating initial structure : " << endl;
  cout << "Number density : " << p.density << " | " << "Box side : " << p.box_dim <<endl;
  cout << "Input molecule types are : " << p.inp_mols << endl;
  cout << "Molecules of solute : " << p.solute_mols << " and molecules of solvent : " << p.solvent_mols << " from a total of " << p.tot_molec << endl;
  return p;
}

void allocateMemory(struct Molecule *mol, struct params &s)
{
    // MEMORY ALLOCATION
 // cout << "Total atoms in the system are " << s.total_atoms_in_system << endl;
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

	cout << "Memory allocated!" << endl;
}


void checkAndOpenFilesInDirectory(const char *path, params &s, struct Molecule *mol)
{
  struct dirent *entry;
  struct stat filestat;
  int mols = 0;
  string line;
  s.total_atoms_in_system = 0;

  DIR *dir = opendir(path);

  if (dir == NULL)
  {
    cerr << "Specified directory is empty!" << endl;
    exit(1);
  }

  while ((entry = readdir(dir)) != NULL)
  {
    //string inf[2];
    cerr << mols << endl;
    stat(entry->d_name, &filestat);
    if (S_ISREG(filestat.st_mode))
    {
      string a = entry->d_name;
      ifstream file;
      file.open(a.c_str());
      getline(file, line); //gets num of atoms
      mol[mols].num_of_atoms = abs(atoi(line.c_str()));
      //cerr << "Num of atoms " << mol[mols].num_of_atoms << endl;

      if (!file) //unable to open file error
      {
        cerr << "Unable to open file " << a << endl;
        cerr << "Check if it exists in the specified folder!" << endl;
        exit(1); // call system to stop
      }
      (s.file_name).push_back(a);
      //ifstream f(inf[s.inp_mols].c_str(), ios::in);
      //(s.inFile).push_back(&f);
      mols++;
      file.close();
    }

  }

  closedir(dir);
  cout << "No of files " << mols << endl;

  for (int i = 0; i < s.inp_mols ; i++)
  {
      if (s.inp_mols == 1)
      {
        mol[i].num_of_molecules = s.solvent_mols;
      }
      else if ( s.inp_mols == 2)
      {
        if (i == 0)
        {
          mol[i].num_of_molecules = s.solute_mols;
        } else if (i == 1)
        {
          mol[i].num_of_molecules = s.solvent_mols;
        }
      }

      mol[i].total_atoms = mol[i].num_of_atoms * mol[i].num_of_molecules;

      cerr << "Number of molecules " << mol[i].num_of_molecules << endl;
      s.total_atoms_in_system += mol[i].total_atoms;
      cout << "Molecule " << i << " has atoms " << mol[i].total_atoms << endl;
  }
}

void readMoleculeXyzFile(struct params &s, struct Molecule *mol)
{
    //Reading single molecule file
    istringstream iss;
    string line;

    for (int j = 0; j < s.inp_mols; j++)
    {
      cerr << "Filename : " << s.file_name[j] << endl;
      ifstream inFile(s.file_name[j]);
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
    cout << "Done reading files" << endl;
}

bool compareMolecules(Molecule lhs, Molecule rhs)
{
  return lhs.num_of_atoms > rhs.num_of_atoms;
}

void sortMolecules(struct Molecule *mol, struct params &s)
{
  sort(mol, mol+s.inp_mols, compareMolecules);
  for (int i = 0 ; i < s.inp_mols; i++)
  {
    cerr << "Molecules sorted!" << endl;
    cerr << mol[i].num_of_atoms << endl;
  }
}

//function to calculate distance between coms
double distChecker(double xyz_temp[][3], vector<coordinates> coords, struct params &s, int j, int k)
{
	double DIST;
	double dist_comp[3];

    //cerr << xyz_temp[j][dim] << endl;
		dist_comp[0] = coords[k].x - xyz_temp[j][0];
		dist_comp[1] = coords[k].y - xyz_temp[j][1];
		dist_comp[2] = coords[k].z - xyz_temp[j][2];

	DIST = sqrt( (dist_comp[0]*dist_comp[0])  + (dist_comp[1]*dist_comp[1]) + (dist_comp[2]*dist_comp[2]) );

	DIST -= static_cast<int>(DIST*s.r_box_dim + 0.5)*s.box_dim ;		//effective distance with PBC
	DIST = fabs(DIST);											                        //absolute value of distance
//	cerr << "Distance is : " << DIST << endl;
	return DIST;
}

//function to return new coords of atoms after PBC implementation
double makePeriodic(double a, double length)
{
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

void findBarycenter(struct Molecule *mol, struct params &s)
{
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
    /*
    for (int dim = 0; dim < s.dimensions; dim++)
    {
      //min[i][dim] = *min_element(mol[i].xyz_co[dim], mol[i].xyz_co[dim] + mol[i].num_of_atoms);
      //cerr << mol[i].num_of_atoms << endl;
      //max[i][dim] = *max_element(mol[i].xyz_co[dim], mol[i].xyz_co[dim] + mol[i].num_of_atoms);
      //cerr << "min " << min[i][dim] << " & max " << max[i][dim] << endl;
      //mol[i].barycenter[dim] = (max[i][dim] + min[i][dim])/2.0;

      min[i][dim] = *min_element(xyz[dim], xyz[dim] + mol[i].num_of_atoms);
      //cerr << mol[i].num_of_atoms << endl;
      max[i][dim] = *max_element(xyz[dim], xyz[dim] + mol[i].num_of_atoms);
      cerr << "min " << min[i][dim] << " & max " << max[i][dim] << endl;
      mol[i].barycenter[dim] = (max[i][dim] + min[i][dim])/2.0;
    }
    */
  }

}

void generateSystem(struct Molecule *mol, struct params &s)
{
  //RNG generator : generates uniformly distributed floating point number in given range
  //  http://www.cplusplus.com/reference/random/uniform_real_distribution/

  bool generation_successful;
  cout << "Initializing random number generator with seed " << s.seed << endl;
  mt19937 generator(s.seed);
  uniform_real_distribution<double> distribution(0.0, 1.0);
  cout << "PRNG initiated!" << endl;

  findBarycenter(mol,s);
  add_coords(coords);
  int curr_atoms = 1;
  int n = s.inp_mols;
	for (int p = 0; p < s.inp_mols ; p++)
	{
		cout << "Creating " << mol[p].total_atoms << " total atoms of molecule type " << p << endl;
		double min_dist = s.min_dist + (n * 0.3) ;
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
					//check distance between all atoms in temp and mol[p].xyz_rotated
          double DIST = min_dist;

          for (int k = 1; k < coords.size()+1 ; ++k)
          {
           	double d = distChecker(xyz_temp,coords, s, j, k);
           	if ( d < DIST)
             {
                DIST = d;
                //cerr << "Current minimum distance is " << DIST << endl;
             }
          }

          if (DIST < min_dist)
          {
            generation_successful = false;
            cerr << "current minimum distance is lower than given criteria, creating new coordinates " << endl;
          	break;
          }
				} // end for loop j atoms in single molecule
        //cerr <<  i << " molecules created with given criteria of type " << p << endl;
			}
			while (generation_successful == false);
      // now assign temp to main array
      curr_atoms = curr_atoms + mol[p].num_of_atoms ;
      cerr << curr_total_molec << " molecules created! " << endl;
      //cerr << "Current atoms are " << curr_atoms << endl;
      //for (int Atoms = (curr_atoms-mol[p].num_of_atoms) ; Atoms < curr_atoms; Atoms++)
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

void writeOutputXyzFile(struct Molecule *mol, struct params &s, vector<coordinates> &coords)
{
  outFile.open(s.out_file.c_str());
  outFile << s.total_atoms_in_system << endl;
  outFile << "!comment line" << endl;
    for (int i = 1; i <= s.total_atoms_in_system ; i++)
    {
      outFile << coords[i].elems << "\t" << coords[i].x << "\t" << coords[i].y << "\t" << coords[i].z << endl;
    }
  cout << "output file written!" << endl;
  outFile.close();
}

void deallocateMemory(struct Molecule *mol, struct params &s)
{

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

  cout << "Memory deallocated!" << endl;
}

int main(int argc, char *argv[])
{
  params s;
  s = getAndPrintRunParameters(argc, argv);
  Molecule mol[s.inp_mols];
  checkAndOpenFilesInDirectory("/mnt/c/Users/Anjali/Dropbox/Thesis/intial_structure/inp_files", s,mol);
  allocateMemory(mol, s);
  readMoleculeXyzFile(s, mol);

  if (s.inp_mols > 1)
  {
  sortMolecules(mol, s); //sorting by number of atoms in each file
  }

  cout << "Creating " << s.tot_molec << " molecules in a simulation box of side " << s.box_dim << " units!" << endl;
  generateSystem(mol, s);
  writeOutputXyzFile(mol, s, coords);

  cout << "Initial structure created with random coordinates and random orientation!" << endl;

  deallocateMemory(mol, s);
  return 0;
}
//# Initial-Structure
