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
 1. Read 2 files 
 2. Create new array for second file molecules
 3. Remove unnecessary arrays elements and all_elements , coords_atoms, coords_atoms_rot, etc
 4. Write function to check distance between each atom of another molecule (except for its own atoms)
 5. Write function to determine larger molecule from input files
 6. Generate system for larger molecules first with much larger minimum distance (10% of given distance)
 7. Generate system for smaller molecules with given minimum distance requirement
 
 
*/
//  cd /mnt/c/Users/Anjali/Dropbox/Thesis/intial_structure
//	g++ random_struc.cpp -o random_struc.exe
//	./random_struc.exe 0.0002 20 prop_gly.xyz CNT.xyz test1.xyz 5 24

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>

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
  float                 box_dim;
  string                inp_file_1;
  string                inp_file_2;
  string                out_file;
  int                   com_id;
  unsigned long int     seed;
  int                   args;
  double                min_dist;
  int                   tot_molec;
  double                r_box_dim;
};

params getAndPrintRunParameters(int argc, char *argv[])
{
  params p;
  float  f;
  p.density     =   atof(argv[1]);   // number density desired for liquid structure
  p.box_dim     =   atof(argv[2]);   //extents of the sim box + 1
  p.inp_file_1  =   argv[3];
  p.inp_file_2  =   argv[4];
  p.out_file    =   argv[5];
  p.com_id      =   atoi(argv[6]);
  p.seed        =   atoi(argv[7]);
  p.args        =   argc;
  p.min_dist    =  6.0;             //twice the largest distance from com (approx)

  f         =  p.density * p.box_dim*p.box_dim*p.box_dim;

  p.tot_molec   = (int)(f+0.5);    //cast to int for proper rounding and conversion
  p.r_box_dim  = 1.0/p.box_dim;



  cout << "Chosen parameters for creating initial structure : " << endl;
  cout << "Number density : " << p.density << " | " << "Box side : " << p.box_dim <<endl;
  cout << "Chosen com id is : " << p.com_id << endl;
  inFile1.open(p.inp_file_1.c_str());

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
    cout << "./random_struct.exe <float desired_num_density> <float box_length> <string single_molecule_1.xyz>  <string single_molecule_2.xyz> <string output_filename.xyz> <int com_id> <int seed>" << endl;
    cout << " - Number density has units molecules/Angstrom^3. " << endl;
    cout << " - For com_id, provide line number from single_molecule.xyz of the atom at the center of the molecule." << endl;
    cout << " - Provide a seed for random number generator (0 to 4,294,967,295)" << endl;
}

void readMoleculeXyzFile(int atoms, string line, double **&coords_atoms, string *&elements)
{
    //Reading single molecule file
    istringstream   iss;
    int             tot_atoms;

    getline(inFile1, line); 		//gets comment line

    for (int i = 0; i < atoms; i++) 	//assigns coordinates from file to array
    {
      iss.clear();
      getline(inFile1,line);

      if(line.length() == 0)        	//empty line error
      {
        cout << "Error: empty line " << endl;
        exit(3);
      }

      iss.str(line);
      iss >> elements[i];
      iss >> coords_atoms[i][0];
      iss >> coords_atoms[i][1];
      iss >> coords_atoms[i][2];

      cout << elements[i] << "|" << coords_atoms[i][0] << "|" << coords_atoms[i][1] << "|" << coords_atoms[i][2] << endl;
    }
    cout << "done reading file" << endl;
    inFile1.close();
}


void allocateMemory(struct params &s, string *&all_elements, string *&elements, double **&coords_atoms_rot, double **&coords_atoms, double **&coords_com, double **&angles_com, int tot_atoms, int atoms)
{
    // MEMORY ALLOCATION
  coords_com = new(nothrow) double*[s.tot_molec];         //xyz coordinates
  angles_com = new(nothrow) double*[s.tot_molec];         //angles
              for (int i = 0; i < s.tot_molec; i++)
              {
                coords_com[i] = new double[3];
				angles_com[i] = new double[3];
              }

  elements = new(nothrow) string[atoms];                        //elements from single molecule file
  all_elements= new(nothrow) string[tot_atoms];                //elements to write in .xyz

  coords_atoms = new(nothrow) double*[atoms];
              for (int i = 0; i < atoms ; i++)
              {
                coords_atoms[i] = new double[3];               //coordinates from input file
              }

  coords_atoms_rot = new(nothrow) double*[tot_atoms];
              for (int i = 0; i < tot_atoms; i++)                  //array for rotated atoms (FINAL)
              {
                coords_atoms_rot[i] = new double[3];
              }
	cout << "Memory allocated!" << endl;  
  }

//function to calculate distance between coms
double distCheck(int i, int j, double  **&coords_com, double r_box_dim, double box_dim)
{

  double dist_comp_x, dist_comp_y, dist_comp_z, DIST;

  dist_comp_x = (coords_com[j][0] - coords_com[i][0])  ;
  dist_comp_y = (coords_com[j][1] - coords_com[i][1])  ;
  dist_comp_z = (coords_com[j][2] - coords_com[i][2])  ;

  DIST = sqrt( (dist_comp_x*dist_comp_x)  + (dist_comp_y*dist_comp_y) + (dist_comp_z*dist_comp_z) );

  DIST -= static_cast<int>(DIST*r_box_dim + 0.5)*box_dim ;		//effective distance with PBC
  DIST = fabs(DIST);											//r=absolute value of distance
  //cout << "Distance is : " << DIST << endl;
  return DIST;
}

void generateSystem(struct params &s, double **&coords_com, double **&angles_com)
{
	 //RNG generator : generates uniformly distributed floating point number in given range
    //  http://www.cplusplus.com/reference/random/uniform_real_distribution/
	//random_device device;           //for random seed generation


    bool      generation_successful;
    cout << "Initializing random number generator with seed " << s.seed << endl;

    mt19937 generator(s.seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);

    cout << "PRNG initiated!" << endl;

   // asssign random coordinates to molecule 1
   
	coords_com[0][0] = distribution(generator) * s.box_dim ;
	coords_com[0][1] = distribution(generator) * s.box_dim ;
	coords_com[0][2] = distribution(generator) * s.box_dim ;

  cout << "coords 0: " <<coords_com[0][0] << "|" << coords_com[0][1] << "|" << coords_com[0][2] << endl;

  angles_com[0][0] = (distribution(generator)  * 2 * pi);
  angles_com[0][1] = (distribution(generator)  * 2 * pi);
  angles_com[0][2] = (distribution(generator)  * 2 * pi);

  cout << "angles 0: " << angles_com[0][0] << "|" << angles_com[0][1] << "|" << angles_com[0][2] << endl;


  int i = 1;
  while (i < s.tot_molec)
  {
	  
    cout <<  "total_molecules are " << i << endl;
    generation_successful = true;

    // assign random coordinates to all molec-1
    coords_com[i][0] = distribution(generator) * s.box_dim ;
    coords_com[i][1] = distribution(generator) * s.box_dim ;
    coords_com[i][2] = distribution(generator) * s.box_dim ;

    cout << "coords " << i << ": " <<coords_com[i][0] << "|" << coords_com[i][1] << "|" << coords_com[i][2] << endl;

    // assign random angles to all molec-1
    angles_com[i][0] = (distribution(generator)  * 2 * pi);
    angles_com[i][1] = (distribution(generator)  * 2 * pi);
    angles_com[i][2] = (distribution(generator)  * 2 * pi);

    cout << "angles  " << i << ": " << angles_com[i][0] << "|" << angles_com[i][1] << "|" << angles_com[i][2] << endl;


    //checking distance for all coms created and storing coordinates only when the distance is greater than minimum distance
    for (int j = 0; j < i; j++ )
    {
        double DIST;

        DIST = distCheck(i, j, coords_com, s.r_box_dim, s.box_dim);

        if (DIST < s.min_dist)                           //sets generation_successful as false
        {
            generation_successful = false;
        }
    }

    if (generation_successful == true)               //increment molecule number (i)
        {
            i++;
        }
  }
	
}

//function to return new coords of atoms after PBC implementation
void makePeriodic(double **&coords_atoms_rot, int curr_atoms, double length, int i)
{	
    if (coords_atoms_rot[curr_atoms][i] < 0.0)
    {	
		coords_atoms_rot[curr_atoms][i] = coords_atoms_rot[curr_atoms][i] + length;
    }
    else if (coords_atoms_rot[curr_atoms][i] > length)
    {
		coords_atoms_rot[curr_atoms][i] = coords_atoms_rot[curr_atoms][i] - length;  
    }
    else
    {
		coords_atoms_rot[curr_atoms][i] = coords_atoms_rot[curr_atoms][i];
	}
}


void rotateMoleculeAboutCOM(struct params &s, double **&coords_atoms, int atoms, double **&coords_atoms_rot, string *&all_elements, string *&elements, double **&coords_com, double **&angles_com)  //function to rotate molecule about chosen COM
{
    cout << "Inside the rotation function!" << endl;
    cout << "Printing function arguments " << endl;
    cout << "atoms  = " << atoms << endl;
    cout << "chosen com is " << s.com_id << endl;
    //dummy stores coordinates of chosen com from single molecule file
    double dummy[3]{};
    for (int i = 0; i < 3 ; i++)
      {
        cout << "dums " << dummy[i] << endl;
        dummy[i] = coords_atoms[s.com_id-3][i];
        cout << "Chosen com coordinates are: " << dummy[i] << endl ;
      }

    for (int i = 0; i < atoms; i++)     //translates molecule to origin based on user's chosen com_id
    {
      for (int j = 0; j < 3; j++)
      {
      coords_atoms[i][j] = coords_atoms[i][j] - dummy[j];
      }
    }

    int k = 0;                  //number of molecules
    int curr_atoms = 0;         //counter for total atoms
    double coords_atoms_temp[atoms][3]{};          // temp array
  //  double angles_com[s.tot_molec][3]{};
    double sin_a[s.tot_molec];  double sin_b[s.tot_molec];  double sin_g[s.tot_molec];
    double cos_a[s.tot_molec];  double cos_b[s.tot_molec];  double cos_g[s.tot_molec];

    while (k < s.tot_molec)
    {
/*
      angles_com[k][0] = (distribution(generator)  * 2 * pi);
      angles_com[k][1] = (distribution(generator)  * 2 * pi);
      angles_com[k][2] = (distribution(generator)  * 2 * pi);
*/
      sin_a[k] = sin(angles_com[k][0]);  sin_b[k] = sin(angles_com[k][1]);  sin_g[k] = sin(angles_com[k][2]);
      cos_a[k] = cos(angles_com[k][0]);  cos_b[k] = cos(angles_com[k][1]);  cos_g[k] = cos(angles_com[k][2]);

      for (int j = 0; j < atoms; j++)
      {
        // rotates molecules about origin using 3d rotation matrix R = Rx(g)*Ry(b)*Rz(a) ; u' = R*u

        coords_atoms_temp[j][0] =     (coords_atoms[j][0] ) * cos_a[k]  * cos_b[k]  +
                                      (coords_atoms[j][1] ) * (cos_a[k] * sin_b[k]  * sin_g[k] - sin_a[k] * cos_g[k]) +
                                      (coords_atoms[j][2] ) * (cos_a[k] * sin_b[k]  * cos_g[k] + sin_a[k] * sin_g[k]) ;


        coords_atoms_temp[j][1] =     (coords_atoms[j][0] ) * sin_a[k]  * cos_b[k]  +
                                      (coords_atoms[j][1] ) * (sin_a[k] * sin_b[k]  * sin_g[k] + cos_a[k] * cos_g[k]) +
                                      (coords_atoms[j][2] ) * (sin_a[k] * sin_b[k]  * cos_g[k] - cos_a[k] * sin_g[k])  ;


        coords_atoms_temp[j][2] =     (coords_atoms[j][0] ) * (-1)      * sin_b[k]  +
                                      (coords_atoms[j][1] ) * cos_b[k]  * sin_g[k]  +
                                      (coords_atoms[j][2] ) * cos_b[k]  * cos_g[k]  ;


        for (int i = 0; i < 3; i++)
        {
          //take rotated atoms and move them to generated com locations
          coords_atoms_rot[curr_atoms][i] = coords_com[k][i] + coords_atoms_temp[j][i];

          //implement periodicity and store final coordinates
          makePeriodic(coords_atoms_rot, curr_atoms, s.box_dim, i);
        }

        //repeat elements after every loop for each molecule
        all_elements[curr_atoms] = elements[j];
        curr_atoms++;
      }
        cout << "total molecules are : " << k << endl;
      k++;
    }

}

void writeOutputXyzFile(struct params &s, string *&all_elements, double **&coords_atoms_rot, int tot_atoms) //writing to .xyz file
{
  outFile.open(s.out_file.c_str());
  outFile << tot_atoms << endl;
  outFile << "!comment line" << endl;
  for (int i = 0; i < tot_atoms ; i++)
  {
    outFile << all_elements[i] << "\t" << coords_atoms_rot[i][0] << "\t" << coords_atoms_rot[i][1] << "\t" << coords_atoms_rot[i][2] << endl;
  }
  cout << "output file written!" << endl;
  outFile.close();
}


void deallocateMemory(struct params &s, string *&all_elements, string *&elements, double **&coords_atoms_rot, double **&coords_atoms, double **&coords_com, double **&angles_com, int tot_atoms, int atoms)
{
  //DEALLOCATE MEMORY

  delete [] elements;
  delete [] all_elements;

  for (int i = 0; i < tot_atoms; i++)                  //array for rotated atoms (FINAL)
  {
    delete [] coords_atoms_rot[i];
  }
  delete [] coords_atoms_rot;

  for (int i = 0; i < atoms ; i++)
  {
    delete [] coords_atoms[i];               //coordinates from input file
  }
  delete [] coords_atoms;

  for (int i = 0; i < s.tot_molec; i++)
  {
    delete [] coords_com[i];
    delete [] angles_com[i];
  }
  delete [] coords_com;
  delete [] angles_com;
}

int main (int argc, char *argv[])
{
	
	if (argc != 8)  //parameter check
  {
    printHelp();
    return 0;
  }

  params s;
  s = getAndPrintRunParameters(argc, argv);

  if (s.tot_molec < 1)
  {
  cout << "Error : Number density too low!" << endl ;
  return 0;
  }

  int             atoms, tot_atoms;
  string          line;

  getline(inFile1,line);     //gets first line with number of atoms
  atoms = abs(atoi(line.c_str()));

  tot_atoms = atoms*s.tot_molec;

  cout << "Total atoms = " << tot_atoms << " Atoms = " << atoms << endl;
    
  // MEMORY ALLOCATION
	double      **coords_com;
	double      **angles_com;
	string      *elements;
	string      *all_elements;
	double      **coords_atoms;
	double      **coords_atoms_rot;

	allocateMemory(s,all_elements,elements,coords_atoms_rot,coords_atoms,coords_com,angles_com,tot_atoms,atoms);
	readMoleculeXyzFile(atoms, line, coords_atoms, elements);

	cout << "Creating " << s.tot_molec << " molecules in a simulation box of side " << s.box_dim << " units!"<< endl;
	generateSystem(s,coords_com,angles_com);
/*
  //RNG generator : generates uniformly distributed floating point number in given range
    //  http://www.cplusplus.com/reference/random/uniform_real_distribution/
	//random_device device;           //for random seed generation


    bool      generation_successful;
    cout << "Initializing random number generator with seed " << s.seed << endl;

    mt19937 generator(s.seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);

    cout << "PRNG initiated!" << endl;

   // asssign random coordinates to molecule 1
  coords_com[0][0] = distribution(generator) * s.box_dim ;
  coords_com[0][1] = distribution(generator) * s.box_dim ;
  coords_com[0][2] = distribution(generator) * s.box_dim ;

  cout << "coords 0: " <<coords_com[0][0] << "|" << coords_com[0][1] << "|" << coords_com[0][2] << endl;

  angles_com[0][0] = (distribution(generator)  * 2 * pi);
  angles_com[0][1] = (distribution(generator)  * 2 * pi);
  angles_com[0][2] = (distribution(generator)  * 2 * pi);

  cout << "angles 0: " << angles_com[0][0] << "|" << angles_com[0][1] << "|" << angles_com[0][2] << endl;


  int i = 1;
  while (i < s.tot_molec)
  {
    cout <<  "total_molecules are " << i << endl;
    generation_successful = true;

    // assign random coordinates to all molec-1
    coords_com[i][0] = distribution(generator) * s.box_dim ;
    coords_com[i][1] = distribution(generator) * s.box_dim ;
    coords_com[i][2] = distribution(generator) * s.box_dim ;

    cout << "coords " << i << ": " <<coords_com[i][0] << "|" << coords_com[i][1] << "|" << coords_com[i][2] << endl;

    // assign random angles to all molec-1
    angles_com[i][0] = (distribution(generator)  * 2 * pi);
    angles_com[i][1] = (distribution(generator)  * 2 * pi);
    angles_com[i][2] = (distribution(generator)  * 2 * pi);

    cout << "angles  " << i << ": " << angles_com[i][0] << "|" << angles_com[i][1] << "|" << angles_com[i][2] << endl;


    //checking distance for all coms created and storing coordinates only when the distance is greater than minimum distance
    for (int j = 0; j < i; j++ )
    {
        double DIST;

        DIST = Dist_check(i, j, coords_com, s.r_box_dim, s.box_dim);

        if (DIST < s.min_dist)                           //sets generation_successful as false
        {
            generation_successful = false;
        }
    }

    if (generation_successful == true)               //increment molecule number (i)
        {
            i++;
        }
  }
*/
	//cout << "Attempting to rotate molecule. " << endl;
	
	rotateMoleculeAboutCOM(s, coords_atoms,atoms,coords_atoms_rot,all_elements,elements,coords_com,angles_com);
	//cout << "Molecule rotated about com!" << endl;

	writeOutputXyzFile(s,all_elements,coords_atoms_rot,tot_atoms);

	cout << "Initial structure created with random coordinates and random orientation!" << endl;

  //DEALLOCATE MEMORY
	
	deallocateMemory(s,all_elements,elements, coords_atoms_rot, coords_atoms, coords_com, angles_com,tot_atoms,atoms);
	cout << "Memory deallocated!" << endl;
	
	return 0;
}
//# Initial-Structure
