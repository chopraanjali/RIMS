/* Program to initialize system structure using random placement approach for molecules

  - the molecules are placed randomly and with a minimum distance of 6 units
  - writes output file in xyz format

    * requires single molecule File in xyz format
    * Parameters to be given:
      - desired number density
      - extents of simulation box
      - input filename and output filename
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>


using namespace std;

//function to return new coords of atoms after PBC implementation
double make_periodic(double x[][3], int c, double length, int b)
{

    if (x[c][b] < 0.0)
    {
      return x[c][b] + length;
    }
    else if (x[c][b] > length)
    {
      return x[c][b] - length;
    }
    else
    {
      return x[c][b];
    }
}

//function to calculate distance between coms
double Dist_check(int a, int b, double pos[][3], double r_length, double length)
{

  double dist_comp_x, dist_comp_y, dist_comp_z, DIST;

  dist_comp_x = (pos[b][0] - pos[a][0])  ;
  dist_comp_y = (pos[b][1] - pos[a][1])  ;
  dist_comp_z = (pos[b][2] - pos[a][2])  ;

  DIST = sqrt( (dist_comp_x*dist_comp_x)  + (dist_comp_y*dist_comp_y) + (dist_comp_z*dist_comp_z) );

  DIST -= static_cast<int>(DIST*r_length + 0.5)*length ;
  return DIST;
}

double TranslateToOrigin(double coords_atoms[][3], double dummy[3], int a, int b)
{
  return coords_atoms[a][b] = (coords_atoms[a][b] - dummy[b]);
}

double TranslateToCOM(int atom_count, double coords_com[][3], double coords_atoms_temp[][3], int a, int b, int c)
{

  return (coords_com[a][c] + coords_atoms_temp[b][c]);
}

double deg_to_rad(double deg) //function to convert degrees to radian
{
  double const pi = 3.14159265358979323;
  double rad = deg * (pi/180.0);
  return rad;
}

struct params
{
  float density;
  float box_dim;
  string inp_file;
  string out_file;
  int com_id;
  unsigned long int seed;
  int args;
};

params GetAndPrintRunParameters(int argc, char *argv[])
{
  params p;
  p.density     =  atof(argv[1]);   // number density desired for liquid structure
  p.box_dim     =  atof(argv[2]);   //extents of the sim box + 1
  p.inp_file    = argv[3];
  p.out_file    = argv[4];
  p.com_id      = atoi(argv[5]);
  p.seed        = atoi(argv[6]);
  p.args        = argc;
  cout << "Chosen parameters for creating intial structure : " << endl;
  cout << "Number density : " << p.density << " | " << "Box side : " << p.box_dim <<endl;

  return p;
}

void PrintHelp()
{
    cout << "Not enough arguments!" << endl;
    cout << "Try with: " << endl;
    cout << "./random_struct.exe <float desired_num_density> <float box_length> <string single_molecule.xyz> <string output_filename.xyz> <int com_id> <int seed>" << endl;
    cout << " - Number density has units molecules/Angstrom^3. " << endl;
    cout << " - For com_id, provide line number from single_molecule.xyz of the atom at the center of the molecule." << endl;
    cout << " - Provide a seed for random number generator (0 to 4,294,967,295)" << endl;
}


int main (int argc, char *argv[])
{

  if (argc != 7)  //parameter check
  {
    PrintHelp();
    return 0;
  }

  params s;
  s = GetAndPrintRunParameters(argc, argv);

  double   min_dist    =  6.0;             //twice the largest distance from com (approx)
  float      f         =  s.density * s.box_dim*s.box_dim*s.box_dim;
  int      tot_molec   = (int)(f+0.5);    //cast to int for proper rounding and conversion
  double    r_box_dim  = 1.0/s.box_dim;

  if (tot_molec < 1)
  {
  cout << "Error : Number density too low!" << endl ;
  return 0;
  }

  cout << "Creating " << tot_molec << " molecules in a simulation box of side " << s.box_dim << " units!"<< endl;

  double    coords_com[tot_molec][3]{};          //xyz coordinates
  double    angles_com[tot_molec][3]{};          //angles

  bool      flag;


  //RNG generator : generates uniformly distributed floating point number in given range
    //  http://www.cplusplus.com/reference/random/uniform_real_distribution/

  //random_device device;           //for random seed generation

  mt19937 generator(s.seed);
  uniform_real_distribution<double> distribution(0.0, 1.0);

   // asssign random coordinates to molecule 1
  coords_com[0][0] = distribution(generator) * s.box_dim ;
  coords_com[0][1] = distribution(generator) * s.box_dim ;
  coords_com[0][2] = distribution(generator) * s.box_dim ;

  angles_com[0][0] = deg_to_rad (distribution(generator)  * 360.0);
  angles_com[0][1] = deg_to_rad (distribution(generator)  * 360.0);
  angles_com[0][2] = deg_to_rad (distribution(generator)  * 360.0);


  int i = 1;
  while (i < tot_molec)
  {
    flag = true;

    // assign random coordinates to all molec-1
    coords_com[i][0] = distribution(generator) * s.box_dim ;
    coords_com[i][1] = distribution(generator) * s.box_dim ;
    coords_com[i][2] = distribution(generator) * s.box_dim ;

    // assign random angles to all molec-1
    angles_com[i][0] = deg_to_rad (distribution(generator)  * 360.0);
    angles_com[i][1] = deg_to_rad (distribution(generator)  * 360.0);
    angles_com[i][2] = deg_to_rad (distribution(generator)  * 360.0);


    //checking distance for all coms created and storing coordinates only when the distance is greater than minimum distance
    for (int j = 0; j < i; j++ )
    {
      double DIST;

      DIST = Dist_check(i, j, coords_com, r_box_dim, s.box_dim);

      if (DIST < min_dist)                           //sets flag as false
      {
        flag = false;
      }
    }

    if (flag == true)                            //increment molecule number (i)
    {
    i++;
    }

  }


  ifstream inFile;                                  //input file
  ofstream outFile;                                 //output File
  inFile.open(s.inp_file.c_str());

  if (!inFile)                                    //unable to open file error
  {
    cerr << "Unable to open file " << s.inp_file << endl;
    cerr << "Check if it exists in the current folder!" << endl;
    exit(1);               // call system to stop
  }

  //Reading single molecule file
  istringstream   iss;
  int             atoms, tot_atoms;
  string          line;

  getline(inFile,line);     //gets first line with number of atoms
  atoms = abs(atoi(line.c_str()));

  tot_atoms = atoms*tot_molec;

//  AllocateMemory(tot_molec, atoms, tot_atoms)

  //intializing arrays
  string elements[atoms];                        //elements from single molecule file
  string all_elements[tot_atoms];                //elements to write in .xyz

  double coords_atoms[atoms][3]{};               //coordinates from input file
  double coords_atoms_temp[atoms][3]{};          // temp array
  double coords_atoms_rot[tot_atoms][3]{};       //array for rotated atoms (FINAL)

  double sin_a[tot_molec];  double sin_b[tot_molec];  double sin_g[tot_molec];
  double cos_a[tot_molec];  double cos_b[tot_molec];  double cos_g[tot_molec];


  getline(inFile, line); //gets comment getline

  for (int i = 0; i < atoms; i++) //assigns coordinates from file to array
  {
    iss.clear();
    getline(inFile,line);

    if(line.length() == 0)        //empty line error
    {
      cout << "Error: empty line " << endl;
      exit(3);
    }

    iss.str(line);
    iss >> elements[i];
    iss >> coords_atoms[i][0];
    iss >> coords_atoms[i][1];
    iss >> coords_atoms[i][2];
  }

  //dummy stores coordinates of chosen com from single molecule file
  double dummy[3];
   for (int i = 0; i < 3 ; i++)
  {
  dummy[i] = coords_atoms[s.com_id-3][i];
  //cout << dummy[i] << endl ;
  }

  for (int i = 0; i < atoms; i++)     //translates molecule to origin based on user's chosen com_id
  {
    for (int j = 0; j < 3; j++)
    {
    coords_atoms[i][j] = coords_atoms[i][j] - dummy[j];
    }
  }
//creating a new function
  void RotateMoleculeAboutCOM()
  {
    int k = 0;                  //number of molecules
    int curr_atoms = 0;         //counter for total atoms
    while (k < tot_molec)
    {

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
          coords_atoms_rot[curr_atoms][i] = make_periodic(coords_atoms_rot, curr_atoms, s.box_dim, i);
        }

        //repeat elements after every loop for each molecule
        all_elements[curr_atoms] = elements[j];
        curr_atoms++;
      }

      k++;
    }

  }

  int k = 0;                  //number of molecules
  int curr_atoms = 0;         //counter for total atoms
  while (k < tot_molec)
  {

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
        coords_atoms_rot[curr_atoms][i] = make_periodic(coords_atoms_rot, curr_atoms, s.box_dim, i);
      }

      //repeat elements after every loop for each molecule
      all_elements[curr_atoms] = elements[j];
      curr_atoms++;
    }

    k++;
  }

    //writing to .xyz file
    outFile.open(s.out_file.c_str());
    outFile << tot_atoms << endl;
    outFile << "!comment line" << endl;
    for (int i = 0; i < tot_atoms ; i++)
    {
      outFile << all_elements[i] << "\t" << coords_atoms_rot[i][0] << "\t" << coords_atoms_rot[i][1] << "\t" << coords_atoms_rot[i][2] << endl;
    }


  inFile.close();
  outFile.close();

  cout << "Initial structure created with random coordinates and random orientation!" << endl;

  return 0;
}
# Initial-Structure
