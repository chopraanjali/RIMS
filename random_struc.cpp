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
 - Create function to read input files from folder
 - If statement for only solvent molecules
 - Write function to determine larger molecule from input files (by no of atoms)
 - Write function to check distance between each atom of another molecule (except for its own atoms)
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
#include <filesystem>     //to read files from directory
#include <algorithm>

using namespace std;

//GLOBAL PARAMETERS
ifstream inFile1;                                  	//input file1
ifstream inFile2;                                 	//input file2
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
  string                inp_file_1;
  string                inp_file_2;
  string                out_file;
  unsigned long int     seed;
  int                   com_id;
  int                   args;
  int                   tot_molec;
  int                   dimensions;
  int                   inp_mols;
  int                   solute_mols;
  int                   solvent_mols;
  int                   total_atoms_in_system;
 };

 struct Molecule
{
  int num_of_atoms;   // total # of atoms in file
  int num_of_molecules;    //total molecules to create for each type
  int total_atoms;    // molecules*atoms in each molecule
  string *elems;        // elems - array for elements
  double **xyz_co;      //array for storing atom coordinates
  double **com_co;      //array for generated com coordinates
  double **com_angles; //array for generated com angles
  double *barycenter;
  double **xyz_rotated;
} ;

params getAndPrintRunParameters(int argc, char *argv[])
{
  params p;
  float  f;

  p.density     =   atof(argv[1]);   // number density desired for liquid structure
  p.box_dim     =   atof(argv[2]);   //extents of the sim box
  p.com_id      =   atoi(argv[3]);
  p.seed        =   atoi(argv[4]);
  p.out_file    =   argv[5];
  p.inp_file_1  =   argv[6];
  //if (argc > 7)
  //{
  p.inp_file_2  =   argv[7];
  p.percent_solute = atof(argv[8]);
  //}
  p.args        =   argc;
  p.min_dist    =   6.0;             //twice the largest distance from com (approx)
  p.dimensions  =   3;
  p.inp_mols    = argc-7;

  f             =  p.density * p.box_dim*p.box_dim*p.box_dim;
  p.tot_molec   = (int)(f+0.5);    //cast to int for proper rounding and conversion
  p.r_box_dim   = 1.0/p.box_dim;


  p.solute_mols = (int)((p.tot_molec * p.percent_solute / 100) + 0.5);
  if (p.inp_mols == 1)
  {
    p.solute_mols = 0;
  }

  p.solvent_mols = p.tot_molec - p.solute_mols ;


  cout << "Chosen parameters for creating initial structure : " << endl;
  cout << "Number density : " << p.density << " | " << "Box side : " << p.box_dim <<endl;
  cout << "Chosen com id is : " << p.com_id << endl;
  cout << "Input molecule types are : " << p.inp_mols << endl;
  cout << "Molecules of solute : " << p.solute_mols << "and molecules of solvent : " << p.solvent_mols << " from a total of " << p.tot_molec << endl;
  inFile1.open(p.inp_file_1.c_str());
  inFile2.open(p.inp_file_2.c_str());

  if (!inFile1)                                    //unable to open file error
  {
    cerr << "Unable to open file " << p.inp_file_1 << endl;
    cerr << "Check if it exists in the current folder!" << endl;
    exit(1);               // call system to stop
  }
  else if (!inFile2)
  {
    cerr << "Unable to open file " << p.inp_file_2 << endl;
    cerr << "Check if it exists in the current folder!" << endl;
    exit(1);               // call system to stop
  }

  return p;
}



void printHelp()
{
    cout << "Not enough arguments!" << endl;
    cout << "Try with: " << endl;
    cout << "./random_struct.exe <float desired_num_density> <float box_length> <int com_id> <int seed> <string output_filename.xyz> <string single_molecule_1.xyz> <string single_molecule_2.xyz> " << endl;
    cout << "For example: " << endl;
    cout << "./random_struc.exe 0.0002 20 5 24 test1.xyz prop_gly.xyz CNT.xyz" << endl;
    cout << " - Number density has units molecules/Angstrom^3. " << endl;
    cout << " - For com_id, provide line number from single_molecule.xyz of the atom at the center of the molecule." << endl;
    cout << " - Provide a seed for random number generator (0 to 4,294,967,295)" << endl;
}

void readMoleculeXyzFile(struct params &s, struct Molecule *mol, string line)
{
    //Reading single molecule file
    istringstream iss;
    getline(inFile1, line); 		//gets comment line
    getline(inFile2, line); 		//gets comment line from second file

    for (int j = 0; j < s.inp_mols; j++)
    {
      for (int i = 0; i < mol[j].num_of_atoms; i++) 	//assigns coordinates from file to array
      {
        iss.clear();

        if (j == 0) {
          getline(inFile1,line);
        } else if (j == 1) {
          getline(inFile2,line);
        }


        if(line.length() == 0)        	//empty line error
        {
          cout << "Error: empty line " << endl;
          exit(3);
        }

        iss.str(line);
        iss >> mol[j].elems[i];
        iss >> mol[j].xyz_co[i][0];
        iss >> mol[j].xyz_co[i][1];
        iss >> mol[j].xyz_co[i][2];

        //cout << mol[j].elems[i] << "|" << mol[j].xyz_co[i][0] << "|" << mol[j].xyz_co[i][1] << "|" << mol[j].xyz_co[i][2] << endl;
      }
    }

    cout << "done reading files" << endl;
    inFile1.close();
    inFile2.close();
}

void allocateMemory(struct Molecule *mol, struct params &s)
{
    // MEMORY ALLOCATION
 // cout << "Total atoms in the system are " << s.total_atoms_in_system << endl;
  for (int i = 0; i < s.inp_mols; i++)
  {

    mol[i].elems = new string[mol[i].total_atoms];
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
    mol[i].xyz_rotated = new(nothrow) double*[mol[i].total_atoms];
                    for (int j = 0; j < mol[i].total_atoms ; j++)
                    {
                      mol[i].xyz_rotated[j] = new double[s.dimensions];
                    }
  }

	cout << "Memory allocated!" << endl;
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
    cout << mol[i].num_of_atoms << endl;
  }
}

//function to calculate distance between coms
double distCheck(int i, int j, struct Molecule *mol, struct params &s, int p)
{

  double DIST;
  double dist_comp[3]{};

  for (int k = 0; k < s.dimensions; k++)
  {
   dist_comp[k] =  mol[p].com_co[j][k] - mol[p].com_co[i][k];
  }

  DIST = sqrt( (dist_comp[0]*dist_comp[0])  + (dist_comp[1]*dist_comp[1]) + (dist_comp[2]*dist_comp[2]) );

  DIST -= static_cast<int>(DIST*s.r_box_dim + 0.5)*s.box_dim ;		//effective distance with PBC
  DIST = fabs(DIST);											                        //absolute value of distance
  //cout << "Distance is : " << DIST << endl;
  return DIST;
}

void generateSystem(struct Molecule *mol, struct params &s)
{
	 //RNG generator : generates uniformly distributed floating point number in given range
    //  http://www.cplusplus.com/reference/random/uniform_real_distribution/
	//random_device device;           //for random seed generation


    bool      generation_successful;
    cout << "Initializing random number generator with seed " << s.seed << endl;

    mt19937 generator(s.seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);
    cout << "PRNG initiated!" << endl;

  int n = s.inp_mols;
  for (int p = 0; p < s.inp_mols ; p++)
  {
    cout << "Creating " << mol[p].total_atoms << " total atoms of molecule type " << p << endl;
    double min_dist = s.min_dist + (n * 0.1) ;
      //cout << "Min distance is " << min_dist << endl;
        // asssign random coordinates and angles to molecule 1 of type 1
        for (int q = 0; q < s.dimensions ; q++)
          {
            mol[p].com_co[0][q] = distribution(generator) * s.box_dim ;
            mol[p].com_angles[0][q] = distribution(generator)  * 2 * pi ;
          }

   // FLAG
    //doesn't work for when only single molecule is to be generated
   int i = 1;
   while (i < mol[p].num_of_molecules)
    {

      cout <<  "total_molecules of type " << p <<  " are " << i << endl;
      generation_successful = true;

        for (int q = 0; q < s.dimensions ; q++)
          {
            mol[p].com_co[i][q] = distribution(generator) * s.box_dim ;
            mol[p].com_angles[i][q] = distribution(generator)  * 2 * pi ;
          }

      //checking distance for all coms created and storing coordinates only when the distance is greater than minimum distance
      for (int j = 0; j < i; j++ )
      {
          double DIST;

          //DIST = distCheck(i, j, coords_com, s.r_box_dim, s.box_dim);
          DIST = distCheck(i, j, mol, s, p);

          if (DIST < min_dist)                           //sets generation_successful as false
          {
              generation_successful = false;
          }
      }

      if (generation_successful == true)               //increment molecule number (i)
          {
              i++;
          }
    }
    n--;
  }

}

//function to return new coords of atoms after PBC implementation
void makePeriodic(struct Molecule *mol, int curr_atoms, double length, int p, int i)
{
    if (mol[i].xyz_rotated[curr_atoms][p] < 0.0)
    {
		mol[i].xyz_rotated[curr_atoms][p] = mol[i].xyz_rotated[curr_atoms][p] + length;
    }
    else if (mol[i].xyz_rotated[curr_atoms][p] > length)
    {
		mol[i].xyz_rotated[curr_atoms][p] = mol[i].xyz_rotated[curr_atoms][p] - length;
    }
    else
    {
		mol[i].xyz_rotated[curr_atoms][p] = mol[i].xyz_rotated[curr_atoms][p];
	}
}

//function to rotate molecule about chosen COM
void rotateMoleculeAboutCOM(struct Molecule *mol, struct params &s)
{
    cout << "Inside the rotation function!" << endl;
    cout << "Printing function arguments " << endl;

    //barycenter stores coordinates of chosen com from single molecule file
    for (int i = 0 ; i < s.inp_mols ; i++)
    {

      for (int p = 0; p < s.dimensions ; p++)
      {
        //FLAG
          //compute actual barycenter
        mol[i].barycenter[p] = mol[i].xyz_co[s.com_id-3][p];
        //cout << "Chosen com coordinates are: " << mol[i].barycenter[p] << endl ;
      }

      for (int j = 0; j < mol[i].num_of_atoms; j++)     //translates molecule to origin based on user's chosen com_id
      {
        for (int p = 0; p < s.dimensions; p++)
        {
          mol[i].xyz_co[j][p] = mol[i].xyz_co[j][p] - mol[i].barycenter[p];
        }
      }


      int k = 0;                  //number of molecules
      int curr_atoms = 0;         //counter for total atoms
      double coords_atoms_temp[mol[i].num_of_atoms][s.dimensions]{};          // temp array
      double sin_a[mol[i].num_of_molecules];  double sin_b[mol[i].num_of_molecules];  double sin_g[mol[i].num_of_molecules];
      double cos_a[mol[i].num_of_molecules];  double cos_b[mol[i].num_of_molecules];  double cos_g[mol[i].num_of_molecules];

      while (k < mol[i].num_of_molecules)
      {
        sin_a[k] = sin(mol[i].com_angles[k][0]);  sin_b[k] = sin(mol[i].com_angles[k][1]);  sin_g[k] = sin(mol[i].com_angles[k][2]);
        cos_a[k] = cos(mol[i].com_angles[k][0]);  cos_b[k] = cos(mol[i].com_angles[k][1]);  cos_g[k] = cos(mol[i].com_angles[k][2]);

        for (int j = 0; j < mol[i].num_of_atoms; j++)
        {
          // rotates molecules about origin using 3d rotation matrix R = Rx(g)*Ry(b)*Rz(a) ; u' = R*u

          coords_atoms_temp[j][0] =     (mol[i].xyz_co[j][0] ) * cos_a[k]  * cos_b[k]  +
                                        (mol[i].xyz_co[j][1] ) * (cos_a[k] * sin_b[k]  * sin_g[k] - sin_a[k] * cos_g[k]) +
                                        (mol[i].xyz_co[j][2] ) * (cos_a[k] * sin_b[k]  * cos_g[k] + sin_a[k] * sin_g[k]) ;


          coords_atoms_temp[j][1] =     (mol[i].xyz_co[j][0] ) * sin_a[k]  * cos_b[k]  +
                                        (mol[i].xyz_co[j][1] ) * (sin_a[k] * sin_b[k]  * sin_g[k] + cos_a[k] * cos_g[k]) +
                                        (mol[i].xyz_co[j][2] ) * (sin_a[k] * sin_b[k]  * cos_g[k] - cos_a[k] * sin_g[k])  ;


          coords_atoms_temp[j][2] =     (mol[i].xyz_co[j][0] ) * (-1)      * sin_b[k]  +
                                        (mol[i].xyz_co[j][1] ) * cos_b[k]  * sin_g[k]  +
                                        (mol[i].xyz_co[j][2] ) * cos_b[k]  * cos_g[k]  ;


          for (int p = 0; p < s.dimensions; p++)
          {
            //take rotated atoms and move them to generated com locations
            mol[i].xyz_rotated[curr_atoms][p] = mol[i].com_co[k][p] + coords_atoms_temp[j][p];

            //implement periodicity and store final coordinates
            makePeriodic(mol, curr_atoms, s.box_dim, p, i);
          }

          //repeat elements after every loop for each molecule
          mol[i].elems[curr_atoms] = mol[i].elems[j];
          curr_atoms++;
        }
        cout << "total molecules are : " << k+1 << endl;
        k++;
      }
    }
}

void writeOutputXyzFile(struct Molecule *mol, struct params &s)
{
  outFile.open(s.out_file.c_str());
  outFile << s.total_atoms_in_system << endl;
  outFile << "!comment line" << endl;
  for (int j = 0; j < s.inp_mols; j++) {
    for (int i = 0; i < mol[j].total_atoms ; i++)
    {
      outFile << mol[j].elems[i] << "\t" << mol[j].xyz_rotated[i][0] << "\t" << mol[j].xyz_rotated[i][1] << "\t" << mol[j].xyz_rotated[i][2] << endl;
    }
  }
  cout << "output file written!" << endl;
  outFile.close();
}

void deallocateMemory(struct Molecule *mol, struct params &s)
{
  //DEALLOCATE MEMORY

  for (int i = 0; i < s.inp_mols; i++)
  {
    delete [] mol[i].elems;

    for (int j = 0; j < mol[i].num_of_atoms; j++)                  //array for rotated atoms (FINAL)
      {
      delete [] mol[i].xyz_co[j];
      }
      delete [] mol[i].xyz_co;

      for (int j = 0; j < mol[i].num_of_molecules; j++)
      {
      delete [] mol[i].com_co[j] ;
      delete [] mol[i].com_angles[j] ;
      }
    delete [] mol[i].com_co ;
    delete [] mol[i].com_angles ;
    delete [] mol[i].barycenter ;
      for (int j = 0; j < mol[i].total_atoms ; j++)
      {
        delete [] mol[i].xyz_rotated[j];
      }
    delete [] mol[i].xyz_rotated;
  }

  cout << "Memory deallocated!" << endl;
}

int main (int argc, char *argv[])
{

	if (argc < 7)  //parameter check
  {
    printHelp();
    return 0;
  }

  params s;
  s = getAndPrintRunParameters(argc, argv);

  Molecule mol[s.inp_mols];

  if (s.tot_molec < 1)
  {
  cout << "Error : Number density too low!" << endl ;
  return 0;
  }

  string          line;
  s.total_atoms_in_system = 0;
  for (int j = 0; j < s.inp_mols; j++)    //doing this here for mem allocation
  {
        if (j == 0) {
          getline(inFile1,line);
          mol[j].num_of_molecules = s.solvent_mols;
        } else if (j == 1) {
          getline(inFile2,line);
          mol[j].num_of_molecules = s.solute_mols;
        }

    mol[j].num_of_atoms  = abs(atoi(line.c_str()));
    mol[j].total_atoms = mol[j].num_of_atoms * mol[j].num_of_molecules;
    s.total_atoms_in_system += mol[j].total_atoms;
    cout << "Molecule " << j << " has atoms " <<  mol[j].num_of_atoms << endl;

  }

	allocateMemory(mol,s);
	readMoleculeXyzFile(s,mol,line);

  cout << "Creating " << s.tot_molec << " molecules in a simulation box of side " << s.box_dim << " units!"<< endl;


  sortMolecules(mol,s);                 //sorting by number of atoms in each file
  generateSystem(mol,s);

	rotateMoleculeAboutCOM(mol,s);
	writeOutputXyzFile(mol, s);

	cout << "Initial structure created with random coordinates and random orientation!" << endl;

  deallocateMemory(mol,s);
	return 0;
}
//# Initial-Structure
