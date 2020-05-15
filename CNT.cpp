/* program to gnerate CNT from graphene structure
 * unit cell with 4 atoms
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc,char *argv[])
{
  //check if parameters are enough

  if ( argc != 4 )
  {
   cout << "not enough arguments!" << endl;
   cout << "try with: " << endl;
   cout << "CNT.exe <double CNT_diameter> <double CNT_length> <string output_filename.xyz>" << endl;
   cout << "NOTE: All units are in nanometers!" << endl;
   return 0;
  }

  int nx,ny,total_atoms; 
  float a;
  double CNT_diameter, CNT_length, len_y, radius, alpha;
  double pi = 3.14159;
  string element_name, output_filename;


  //a  = 1.42 ;	//cc bond length (in Angstroms)
  a   = 0.142 ; // in nm

  CNT_diameter = atof(argv[1]);
  CNT_length = atof(argv[2]);
  output_filename = argv[3];

  element_name = "C";

  len_y = pi*CNT_diameter;                                //graphene sheet in y_dir
  
  radius = CNT_diameter/2;                                 //of CNT
  cout << "Given radius = " << radius << endl;

  nx = (int)((CNT_length / (3.0 * a)) + 0.5);              //number of repititions of unit cell in x direction
  ny = (int)((len_y / (1.732 * a)) +0.5);                  //number of repititions of unit cell in x direction

  radius = (1.732 * a * ny ) / ( 2 * pi);                  //computing again for rounding errors
  cout << "Actual radius = " << radius << endl;
  
  total_atoms = 4 * nx * ny;
  
  cout << "Run parameters: " << endl;
  cout << "Bond length = "<< a << endl;
  cout << "nx = " << nx << endl;
  cout << "ny = " << ny << endl;
  cout << "element name = "<< element_name << endl;
  cout << "file name =  "<< output_filename << endl;

  double coords[total_atoms][4];
  int curr_atoms = 0;

  for (int i = 0; i < nx; i++) {

    for (int j = 0; j < ny; j++) {

     //atom 1 (a,0,0)
      coords[curr_atoms][0] = (i*3*a)+a;
      coords[curr_atoms][1] = j*1.732*a;
      alpha                 = coords[curr_atoms][1]/radius;
      coords[curr_atoms][2] = radius*cos(alpha);
      coords[curr_atoms][3] = radius*sin(alpha);
      curr_atoms++;

     //atom 2 (2a,0,0)
      coords[curr_atoms][0] = (i*3*a)+(2*a);
      coords[curr_atoms][1] = j*1.732*a;
      alpha                 = coords[curr_atoms][1]/radius;
      coords[curr_atoms][2] = radius*cos(alpha);
      coords[curr_atoms][3] = radius*sin(alpha);
      curr_atoms++;

     //atom 3 (a/2, a*(sqrt3)/2, 0)
     //(a*0.5, 1.7*0.5*a, 0)
      coords[curr_atoms][0] = (i*3*a)+(a/2);
      coords[curr_atoms][1] = (j*1.732*a) +(1.732*a/2);
      alpha                 = coords[curr_atoms][1]/radius;
      coords[curr_atoms][2] = radius*cos(alpha);
      coords[curr_atoms][3] = radius*sin(alpha);

      curr_atoms++;

     //atom 4 (2.5*a,1.7*0.5*a,0)
      coords[curr_atoms][0] = (i*3*a)+(a*5/2);
      coords[curr_atoms][1] = (j*1.732*a) +(1.732*a/2.0);
      alpha                 = coords[curr_atoms][1]/radius;
      coords[curr_atoms][2] = radius*cos(alpha);
      coords[curr_atoms][3] = radius*sin(alpha);

      curr_atoms++;
    }

  }

  ofstream output_file;
  output_file.open(output_filename.c_str());

  output_file << curr_atoms << endl;
  output_file << "CNT from Graphene sheet" << endl;

  for (int i = 0; i < curr_atoms; i++)
  {

    output_file << element_name << "	" << coords[i][0] << "	" << coords[i][2] << "	" << coords[i][3] << endl;
  }

  output_file.close();

  return 0;
}
